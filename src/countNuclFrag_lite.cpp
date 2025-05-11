#include "myoperation.h"
#include <filesystem>
#include <getopt.h>
#include <future>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <numeric>
#include <atomic>
#include <mutex>
#include <cstdlib>
#include <cerrno>
#include <climits>
#include <sys/sysinfo.h>
#include <thread>

using namespace std;

// Get available system memory (bytes)
size_t getAvailableMemory() {
    struct sysinfo memInfo;
    if (sysinfo(&memInfo) == 0) return memInfo.freeram;
    return 0;
}

// Adjust thread count based on memory estimate
int adjustThreadNum(int desired, size_t perThreadEstimate) {
    size_t avail = getAvailableMemory();
    if (perThreadEstimate == 0) return max(1, desired);
    int maxTh = static_cast<int>(avail / perThreadEstimate);
    return max(1, min(desired, maxTh));
}

// Assign a single tag to region bins, update count
void assignTagPosLocal(int pos,
                       const vector<pair<int,int>>& regs,
                       vector<int>& tempCnt,
                       int* colSum) {
    for (size_t i = 0; i < regs.size(); ++i) {
        if (regs[i].second < pos) continue;
        if (regs[i].first > pos) break;
        tempCnt[i]++;
        colSum[i]++;
        break;
    }
}

// Simple pool counting using header-provided mapping
void assignCountPoolSimple(const map<string, vector<int>>& tagpos,
                           const count_pool& cp,
                           int* colSum,
                           int colLen) {
    for (const auto& kv : tagpos) {
        const string& chr = kv.first;
        const auto& positions = kv.second;
        if (!cp.region_id_map.count(chr)) continue;
        for (const auto& rm : cp.region_id_map.at(chr)) {
            auto coords = rm.first;
            size_t rid = rm.second;
            const auto& regs = cp.cut_reg_ve[rid];
            auto it = lower_bound(positions.begin(), positions.end(), coords.first);
            for (; it != positions.end() && *it <= coords.second; ++it) {
                static thread_local vector<int> tempCnt;
                tempCnt.assign(colLen, 0);
                assignTagPosLocal(*it, regs, tempCnt, colSum);
            }
        }
    }
}

int main(int argc, char* argv[]) {
    auto t0 = chrono::high_resolution_clock::now();

    // ----------------- parse command-line -----------------
    const char* short_opts = "f:c:o:u:d:w:p:";
    static struct option long_opts[] = {
        {"fragment",   required_argument, nullptr, 'f'},
        {"centerDir",  required_argument, nullptr, 'c'},
        {"out_prefix", required_argument, nullptr, 'o'},
        {"up",         required_argument, nullptr, 'u'},
        {"down",       required_argument, nullptr, 'd'},
        {"window",     required_argument, nullptr, 'w'},
        {"threads",    required_argument, nullptr, 'p'},
        {nullptr,0,nullptr,0}
    };
    string metaFile, regionDir, prefix = "output";
    int ups = 100, downs = 100, win = 1, threads = 1;
    int opt;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'f': metaFile = optarg; break;
            case 'c': regionDir = optarg; break;
            case 'o': prefix = optarg; break;
            case 'u': ups = stoi(optarg); break;
            case 'd': downs = stoi(optarg); break;
            case 'w': win = stoi(optarg); break;
            case 'p': threads = stoi(optarg); break;
            default: break;
        }
    }
    if (metaFile.empty() || regionDir.empty()) {
        cerr << "Usage: " << argv[0]
             << " -f <metaFile> -c <centerDir> -o <out_prefix> \n";
        return 1;
    }

    // --------------- load cell & TF lists ----------------
    vector<string> cellBeds = getPathFromMetafile(metaFile);
    vector<string> tfFiles  = getFiles(regionDir, ".txt", "tfc");
    size_t nCell = cellBeds.size();
    size_t nTF   = tfFiles.size();
    if (nCell == 0 || nTF == 0) {
        cerr << "Error: no cells or no TF files found\n";
        return 1;
    }

    // ------------- prepare result matrices ---------------
    vector<vector<int>> centerMatrix(nCell, vector<int>(nTF, 0));
    vector<vector<int>> flankMatrix (nCell, vector<int>(nTF, 0));

    // -------------- process each TF sequentially ---------------
    for (size_t tfIdx = 0; tfIdx < nTF; ++tfIdx) {
        // 1) TF-level progress
        cout << "Processing TF " << (tfIdx+1) << "/" << nTF
             << ": " << filesystem::path(tfFiles[tfIdx]).stem().string()
             << std::endl;

        // 2) build pools
        count_pool poolCenter, poolFlank;
        map<string, vector<int>> centers;
        readRegionCenter(tfFiles[tfIdx], centers);
        getTFCcpool(centers, poolCenter);
        getTFCcpool(centers, poolFlank);
        poolCenter.getregionidmap_tfc(ups, downs, win);
        poolFlank .getregionidmap_tfc(downs, ups, win);

        int lenCenter = poolCenter.counts_table.empty() ? 0
                        : poolCenter.counts_table[0].size();
        int lenFlank  = poolFlank.counts_table.empty()  ? 0
                        : poolFlank.counts_table[0].size();

        // 3) adjust threads by memory estimate
        size_t estMem = max(lenCenter, lenFlank) * sizeof(int) + 1024;
        int th = adjustThreadNum(threads, estMem);

        // 4) parallel cell processing
        atomic<size_t> idx(0);
        mutex ioM;

        auto worker = [&]() {
            map<string, vector<int>> tagpos;
            while (true) {
                size_t i = idx.fetch_add(1);
                if (i >= nCell) break;

                // read tags once per cell
                readAlignTagFile(cellBeds[i], tagpos);

                // center count
                vector<int> tmpC(lenCenter, 0);
                assignCountPoolSimple(tagpos, poolCenter, tmpC.data(), lenCenter);
                centerMatrix[i][tfIdx] = accumulate(tmpC.begin(), tmpC.end(), 0);

                // flank count
                vector<int> tmpF(lenFlank, 0);
                assignCountPoolSimple(tagpos, poolFlank, tmpF.data(), lenFlank);
                flankMatrix [i][tfIdx] = accumulate(tmpF.begin(), tmpF.end(), 0);

                // cell-level progress every 50 cells
                if ((i+1) % 50 == 0) {
                    lock_guard<mutex> lk(ioM);
                    cout << "  Cell " << (i+1) << "/" << nCell
                         << " done for TF " << (tfIdx+1) << "/" << nTF
                         << std::endl;
                }

                tagpos.clear();
            }
        };

        vector<future<void>> fut;
        for (int t = 0; t < th; ++t)
            fut.emplace_back(async(launch::async, worker));
        for (auto &f : fut) f.get();
    }

    // ---------------- write matrices to file ----------------
    string centerMatFile = prefix + ".centerNuclCountMatrix.txt";
    string flankMatFile  = prefix + ".flankNuclCountMatrix.txt";

    ofstream oc(centerMatFile), of(flankMatFile);
    // header
    oc << "Cell/TF";
    of << "Cell/TF";
    for (size_t tfIdx = 0; tfIdx < nTF; ++tfIdx) {
        string name = filesystem::path(tfFiles[tfIdx]).stem().string();
        oc << "\t" << name;
        of << "\t" << name;
    }
    oc << "\n";
    of  << "\n";

    // rows
    for (size_t i = 0; i < nCell; ++i) {
        string cellName = filesystem::path(cellBeds[i]).filename().string();
        oc << cellName;
        of << cellName;
        for (size_t tfIdx = 0; tfIdx < nTF; ++tfIdx) {
            oc << "\t" << centerMatrix[i][tfIdx];
            of << "\t" << flankMatrix [i][tfIdx];
        }
        oc << "\n";
        of  << "\n";
    }

    if (remove("tmp2.txt") != 0 && errno != ENOENT) {
        cerr << "Warning: Failed to remove tmp2.txt" << endl;
    }
    if (remove("files.txt") != 0 && errno != ENOENT) {
        cerr << "Warning: Failed to remove files.txt" << endl;
    }
    if (remove("tfc_files.txt") != 0 && errno != ENOENT) {
        cerr << "Warning: Failed to remove tfc_files.txt" << endl;
    }

    auto t1 = chrono::high_resolution_clock::now();
    cout << "Total time: "
         << chrono::duration<double>(t1 - t0).count()
         << " s\n";
    return 0;
}
