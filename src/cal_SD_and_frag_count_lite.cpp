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
#include <iomanip>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <cstdlib>
#include <cerrno>
#include <sys/sysinfo.h>
#include <set>

using namespace std;

size_t getAvailableMemory() {
    struct sysinfo memInfo;
    if (sysinfo(&memInfo) == 0) {
        return memInfo.freeram;
    }
    return 0;
}

int adjustThreadNum(int desiredThreads, size_t perThreadEstimate) {
    size_t availMem = getAvailableMemory();
    if (perThreadEstimate == 0) {
        return max(1, desiredThreads);
    }
    int maxThreads = (int)(availMem / perThreadEstimate);
    return max(1, min(desiredThreads, maxThreads));
}



void readAlignTagFile_multiThread(const string &infile, map<string, vector<int>> &tagpos) {
    ifstream inf(infile);
    if (!inf.good()) {
        cerr << "file " << infile << " not found " << endl;
        exit(1);
    }
    while (!inf.eof()) {
        string line;
        getline(inf, line);
        if (line.empty()) {
            break;
        }
        vector<string> parseditem = parseString(line);
        if (parseditem.size() != 6) {
            cerr << "error line: " << line << endl;
            cerr << "By default, you need to use BED6 format with the 6th column indicating strand." << endl;
            exit(1);
        }
        string chr  = parseditem[0];
        int start   = stoi(parseditem[1]);
        int end     = stoi(parseditem[2]);
        char strand = parseditem[5][0];
        if (strand != '+' && strand != '-') {
            cerr << "error strand " << strand << endl;
            exit(1);
        }
        int shiftpos = start + (end - start) / 2;
        tagpos[chr].push_back(shiftpos);
    }
    inf.close();

    for (auto &p : tagpos) {
        sort(p.second.begin(), p.second.end());
    }
}


void assignTagPos_multiThread(int tagpos, 
                              const vector<pair<int, int>> &cut_reg, 
                              vector<int> &counts, 
                              int* colSum, 
                              bool* allHitArr, 
                              double &allHits, 
                              bool* centerHitArr, 
                              double &centerHits, 
                              int sdStartIdx, 
                              int sdEndIdx, 
                              int mymaxReadNum) 
{
    for (size_t i = 0; i < cut_reg.size(); ++i) {
        if (cut_reg[i].second < tagpos) {
            continue;
        }
        if (cut_reg[i].first > tagpos) {
            break;
        }
        counts[i] += 1;
        colSum[i] += 1;

        if (counts[i] > mymaxReadNum) {
            fill(counts.begin(), counts.end(), 0);
            fill(colSum, colSum + cut_reg.size(), 0);
            break;
        }

        if (!allHitArr[i]) {
            allHitArr[i] = true;
            allHits += 1.0;
            if (i >= (size_t)sdStartIdx && i < (size_t)sdEndIdx) {
                centerHitArr[i] = true;
                centerHits += 1.0;
            }
        }
        break;
    }
}


void assignCountPool_multiThread(const map<string, vector<int>> &tagpos, 
                                 const count_pool &cpool, 
                                 int* colSum, 
                                 bool* allHitArr, 
                                 double &allHits,
                                 bool* centerHitArr, 
                                 double &centerHits,
                                 int sdStartIdx, 
                                 int sdEndIdx, 
                                 int maxReadNum, 
                                 int colLen) 
{
    for (auto &kv : tagpos) {
        const string &chr     = kv.first;
        const auto   &posVec  = kv.second;
        if (cpool.region_id_map.count(chr) == 0) {
            continue;
        }
        // region_id_map[chr]: map<pair<int,int>, size_t> 
        // cut_reg_ve[id]     : vector<pair<int,int>>
        
        for (const auto &ri : cpool.region_id_map.at(chr)) {
            const auto &rStartEnd = ri.first;  // pair<int,int>
            size_t regionID       = ri.second;

            const auto &cutRegVec = cpool.cut_reg_ve[regionID];

            int rStart = rStartEnd.first;
            int rEnd   = rStartEnd.second;

            auto lowerIt = std::lower_bound(posVec.begin(), posVec.end(), rStart);
            if (lowerIt == posVec.end()) {
                break;
            }
            for (auto it = lowerIt; it != posVec.end(); ++it) {
                int curPos = *it;
                if (curPos > rEnd) {
                    break;
                }
                // if in [rStart, rEnd], then assignTagPos_multiThread
                static thread_local vector<int> tempCounts; 
                tempCounts.assign(colLen, 0);

                assignTagPos_multiThread(curPos, 
                                         cutRegVec, 
                                         tempCounts,
                                         colSum,
                                         allHitArr, 
                                         allHits,
                                         centerHitArr, 
                                         centerHits,
                                         sdStartIdx, 
                                         sdEndIdx, 
                                         maxReadNum);
            }
        }
    }
}


double get_summitDist_multiThread(int* pileup, int colLen, int motifShift = 200, int nuclShift = 73, 
        			int posScoreOutterFlank = 100, int posScoreInnerFlank = 15, int centerShift = 20) {
    int centerIdx = colLen / 2;
    int psStartIndx = centerIdx - motifShift - nuclShift;
    int psEndIndx = centerIdx + motifShift + nuclShift + 1;

    if (psStartIndx < 0 || psEndIndx > colLen) {
        cout << "Error: index out of range" << endl;
        return -1;
    }

    vector<double> posScore(colLen, 0.0);
    for (int k = psStartIndx; k < psEndIndx; k++) {
        double tmpMidpointSum_innerBin = 0.0;
        double tmpMidpointSum_outterBin = 0.0;
        for (int j = k - posScoreInnerFlank; j <= k + posScoreInnerFlank; j++) {
            if(j < 0 || j >= colLen) continue;
            tmpMidpointSum_innerBin += pileup[j];
        }
        for (int j = k - posScoreOutterFlank; j <= k + posScoreOutterFlank; j++) {
            if(j < 0 || j >= colLen) continue;
            tmpMidpointSum_outterBin += pileup[j];
        }
        if (tmpMidpointSum_outterBin == 0)
            posScore[k] = 0;
        else 
            posScore[k] = tmpMidpointSum_innerBin / tmpMidpointSum_outterBin;
    }
    
    int gwStartIdx = centerIdx - motifShift;
    int gwEndIdx = centerIdx + motifShift + 1;
    
    // calculate gwScore and make sure that the index is in the range
    vector<double> gwScore(colLen, 0.0);
    for (int k = gwStartIdx; k < gwEndIdx; k++) {
        double tmp_k_gwScore = 0.0;
        for (int j = -nuclShift; j <= nuclShift; j++) {
            int idx = k + j;
            if (idx < 0 || idx >= colLen) continue;
            double gw = exp(-pow(j / 20.0, 2) / 2);
            tmp_k_gwScore += posScore[idx] * gw;
        }
        gwScore[k] = tmp_k_gwScore;
    }
    gwScore.assign(gwScore.begin() + gwStartIdx, gwScore.begin() + gwEndIdx);

    int gwSize = gwScore.size();  // expectation: (centerIdx+motifShift+1) - (centerIdx-motifShift) = 2*motifShift+1
    int gwCenterIdx = gwSize / 2;
    vector<double> uCwScore(gwScore.begin(), gwScore.begin() + gwCenterIdx);
    vector<double> dCwScore(gwScore.begin() + gwCenterIdx, gwScore.end());
    
    reverse(uCwScore.begin(), uCwScore.end());
    int umidx_flip = max_element(uCwScore.begin(), uCwScore.end()) - uCwScore.begin();
    int umax_center_shift = 0 - umidx_flip - 1;
    int umidx = uCwScore.size() - umidx_flip - 1;
    
    int dmidx = max_element(dCwScore.begin(), dCwScore.end()) - dCwScore.begin();
    int dmax_center_shift = dmidx;
    
    vector<int> slicedPileup(pileup + gwStartIdx, pileup + gwEndIdx);
    int slicedSize = slicedPileup.size();
    int slicedCenterIndex = slicedSize / 2;
    vector<double> my_array(slicedPileup.begin(), slicedPileup.end());

    // dynamic center shift
    int upstream_start = max(0, slicedCenterIndex - centerShift);
    int upstream_end = slicedCenterIndex;  // [upstream_start, slicedCenterIndex)
    int index_upstream = distance(my_array.begin(), min_element(my_array.begin() + upstream_start, my_array.begin() + upstream_end));
    
    int downstream_start = slicedCenterIndex;
    int downstream_end = min(slicedSize, slicedCenterIndex + centerShift);
    int index_downstream = distance(my_array.begin(), min_element(my_array.begin() + downstream_start, my_array.begin() + downstream_end));
    
    if (my_array[index_upstream] < my_array[index_downstream])
        slicedCenterIndex = index_upstream;
    else
        slicedCenterIndex = index_downstream;
    
    int uMidPointCountSum = 0;
    double uShiftSum = 0;
    for (int u = umidx; u < slicedCenterIndex; u++) {
        int tmpCenterShift = slicedCenterIndex - u;
        int tmpMidPointCount = slicedPileup[u];
        uMidPointCountSum += tmpMidPointCount;
        uShiftSum += tmpMidPointCount * tmpCenterShift;
    }
    double uShiftMean = (uMidPointCountSum == 0) ? -1 * umax_center_shift : uShiftSum / uMidPointCountSum;
    
    int dMidPointCountSum = 0;
    double dShiftSum = 0;
    for (int d = slicedCenterIndex + 1; d < dmidx + slicedCenterIndex && d < slicedSize; d++) {
        int tmpCenterShift = d - slicedCenterIndex;
        int tmpMidPointCount = slicedPileup[d];
        dMidPointCountSum += tmpMidPointCount;
        dShiftSum += tmpMidPointCount * tmpCenterShift;
    }
    double dShiftMean = (dMidPointCountSum == 0) ? dmax_center_shift : dShiftSum / dMidPointCountSum;
    
    double summitDist = uShiftMean + dShiftMean;
    return summitDist;
}



tuple<int, int, double> processOneCellOneTF(const string &cellBed,
                        const count_pool &cpool,
                        int colLen,
                        int ups, 
                        int downs, 
                        int win,
                        int motifShift,
                        double covThre, 
                        int sdReionReadHitThre,
                        int maxReadNum,
                        int center_range = 100,
                        int flank_range  = 200)
{
    map<string, vector<int>> localTagposmap;
    readAlignTagFile_multiThread(cellBed, localTagposmap);

    vector<int> colSum(colLen, 0);
    bool *allHitArr    = new bool[colLen]();
    bool *centerHitArr = new bool[colLen]();
    double allHits     = 0.0;
    double centerHits  = 0.0;

    int centerIdx = colLen / 2;
    int startIdx  = max(0, centerIdx - motifShift);
    int endIdx    = min(colLen, centerIdx + motifShift + 1);

    assignCountPool_multiThread(localTagposmap, cpool,
                                colSum.data(),
                                allHitArr, allHits,
                                centerHitArr, centerHits,
                                startIdx, endIdx,
                                maxReadNum, colLen);

    double nonZeroPct    = allHits / colLen;
    double tmpSummitDist = 0.0;
    if (nonZeroPct >= covThre && centerHits > sdReionReadHitThre) {
        tmpSummitDist = get_summitDist_multiThread(colSum.data(), colLen, motifShift);
    }

    int tmpCenterMatSum = 0;
    for (int i = max(0, centerIdx - center_range); i < min(colLen, centerIdx + center_range + 1); ++i) {
        tmpCenterMatSum += colSum[i];
    }

    int tmpFlankMatSum = 0;
    for (int i = max(0, centerIdx - flank_range); i < min(colLen, centerIdx + flank_range + 1); ++i) {
        tmpFlankMatSum += colSum[i];
    }

    delete[] allHitArr;
    delete[] centerHitArr;
    localTagposmap.clear();

    return make_tuple(tmpCenterMatSum, tmpFlankMatSum, tmpSummitDist);
}



int main(int argc, char* const argv[]) {
    auto start = chrono::high_resolution_clock::now();

    const char* short_options = "f:c:o:u:d:w:m:t:n:r:k:p:";
    static struct option long_options[] = {
        {"fragment",             required_argument, nullptr, 'f'},
        {"center",               required_argument, nullptr, 'c'},
        {"out_prefix",           required_argument, nullptr, 'o'},
        {"up",                   required_argument, nullptr, 'u'},
        {"down",                 required_argument, nullptr, 'd'},
        {"window_size",          required_argument, nullptr, 'w'},
        {"motif_shift",          required_argument, nullptr, 'm'},
        {"coverage_threshold",   required_argument, nullptr, 't'},
        {"region_hit_threshold", required_argument, nullptr, 'n'},
        {"center_range",         required_argument, nullptr, 'r'},
        {"flank_range",          required_argument, nullptr, 'k'},
        {"threads",              required_argument, nullptr, 'p'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    string metaFile, regionPath, prefix = "output";
    int ups = 400, downs = 400, win = 1, motifShift = 200, sdReionReadHitThre = 1, myCenterRange = 100, myFlankRange = 200;
    double covThre = 0.05;
    int threadNum  = 1;

    while ((opt = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1) {
        switch (opt) {
            case 'f':
                metaFile = optarg;
                break;
            case 'c':
                regionPath = optarg;
                break;
            case 'o':
                prefix = optarg;
                break;
            case 'u':
                ups = stoi(optarg);
                break;
            case 'd':
                downs = stoi(optarg);
                break;
            case 'w':
                win = stoi(optarg);
                break;
            case 'm':
                motifShift = stoi(optarg);
                break;
            case 't':
                covThre = stod(optarg);
                break;
            case 'n':
                sdReionReadHitThre = stoi(optarg);
                break;
            case 'r':
                myCenterRange = stoi(optarg);
                break;
            case 'k':
                myFlankRange = stoi(optarg);
                break;
            case 'p':
                threadNum = stoi(optarg);
                break;
            case '?':
                cerr << "Unknown option or missing argument." << endl;
                return 1;
            default:
                break;
        }
    }

    if (metaFile.empty() || regionPath.empty()) {
        cerr << "Error: Missing required options." << endl;
        cerr << "Usage: " << argv[0] << " <options>" << endl;
        cerr << "-f:        A meta file with single-cell BED file paths (1st column) + labels (2nd col)\n"
             << "-c:        TFBS_center directory or file\n"
             << "-o:        out_prefix\n"
             << "-u:        center upstream, default: 400\n"
             << "-d:        center downstream, default: 400\n"
             << "-w:        window_size, default: 1\n"
             << "-m:        motif_shift_size, default: 200\n"
             << "-t:        coverage_threshold, default: 0.05\n"
             << "-n:        region_hit_threshold, default: 1\n"
             << "-r:        center range (upstream and downstream) used in fragment count calculation, default: 100\n"
             << "-k:        flank range (upstream and downstream) used in fragment count calculation, default: 200\n"
             << "-p:        number of threads, default: 1\n";
        return 1;
    }

    if (threadNum < 1) {
        cerr << "Error: Number of threads (-p) must be >= 1." << endl;
        return 1;
    }

    std::filesystem::path parentDir = std::filesystem::path(prefix).parent_path();
    if (!parentDir.empty() && !std::filesystem::exists(parentDir)) {
        std::filesystem::create_directories(parentDir);
    }

    vector<string> cellBedFiles = getPathFromMetafile(metaFile);
    size_t totalCellNum = cellBedFiles.size();

    vector<string> TFCfiles = getFiles(regionPath, ".txt", "tfc");
    size_t totalTFC = TFCfiles.size();

    vector<vector<double>> summitDist_mat(totalCellNum, vector<double>(totalTFC, 0.0));
    vector<vector<int>> center_countSum_mat(totalCellNum, vector<int>(totalTFC, 0));
    vector<vector<int>> flank_countSum_mat(totalCellNum, vector<int>(totalTFC, 0));

    if (TFCfiles.empty()) {
        cerr << "No TFC files found in " << regionPath << " . Exit." << endl;
        return 1;
    }

    int maxReadNum = 5;

    size_t estimatedMemoryPerThread = 1024;

    // for each TF center
    for (size_t tfIdx = 0; tfIdx < totalTFC; tfIdx++) {
        string TFC = TFCfiles[tfIdx];
        size_t pos1 = TFC.find_last_of('/');
        string filename = (pos1 == string::npos) ? TFC : TFC.substr(pos1 + 1);
        size_t pos2 = filename.find_last_of('.');
        string centerName = (pos2 == string::npos) ? filename : filename.substr(0, pos2);

        cout << "# Processing TF center: " << centerName << endl;

        ifstream myfile(TFC);
        if (!myfile.good() || myfile.peek() == ifstream::traits_type::eof()) {
            cerr << "### WARNING: TFC file " << TFC << " is empty or not good. Skipping." << endl;
            continue;
        }

        count_pool cpool;
        {
            map<string, vector<int>> tfcenters;
            readRegionCenter(TFC, tfcenters);
            getTFCcpool(tfcenters, cpool);
            cpool.getregionidmap_tfc(ups, downs, win);
        }


        int colLen = 0;
        if (!cpool.counts_table.empty()) {
            colLen = (int)cpool.counts_table[0].size();
        } else {
            cerr << "### WARNING: cpool.counts_table is empty for " << centerName << ". Skipping." << endl;
            continue;
        }

        // uodate approximate memory usage for each thread
        estimatedMemoryPerThread = (size_t)colLen * (sizeof(int) + 2*sizeof(bool)) + 1024;
        threadNum = adjustThreadNum(threadNum, estimatedMemoryPerThread);
        cout << "Using " << threadNum << " threads for this TF center." << endl;

        // for all the cell bed files
        atomic<size_t> nextIndex(0);
        mutex coutMutex;

        auto worker = [&]() {
            while (true) {
                size_t i = nextIndex.fetch_add(1);
                if (i >= totalCellNum) break;

                auto [centerMatSum_val, flankMatSum_val, summitDist_val] = processOneCellOneTF(cellBedFiles[i], cpool, colLen,ups, downs, win,
                                                                                                motifShift,covThre, sdReionReadHitThre, maxReadNum,
                                                                                                myCenterRange, myFlankRange);
                summitDist_mat[i][tfIdx] = summitDist_val;
                center_countSum_mat[i][tfIdx] = centerMatSum_val;
                flank_countSum_mat[i][tfIdx] = flankMatSum_val;

                {
                    lock_guard<mutex> lock(coutMutex);
                    static size_t processed = 0;
                    processed++;
                    if (processed % 50 == 0) {
                        cout << "Processed " << processed << " / " << totalCellNum 
                             << " cells for TF center: " << centerName << endl;
                    }
                }
            }
        };

        vector<future<void>> futures;
        for (int th = 0; th < threadNum; th++) {
            futures.push_back(std::async(std::launch::async, worker));
        }
        for (auto &f : futures) {
            f.get();
        }
    }



    // output center counts
    string centerfile = prefix + ".centerNuclCountSum.txt";
    ofstream outf_center(centerfile);
    if (!outf_center.good()) {
        cerr << "Cannot open output file: " << centerfile << endl;
        return 1;
    }

    for (size_t tfIdx = 0; tfIdx < TFCfiles.size(); tfIdx++) {
        const string &TFC = TFCfiles[tfIdx];
        size_t pos1 = TFC.find_last_of('/');
        string filename = (pos1 == string::npos) ? TFC : TFC.substr(pos1 + 1);
        size_t pos2 = filename.find_last_of('.');
        string centerName = (pos2 == string::npos) ? filename : filename.substr(0, pos2);

        outf_center << centerName;
        if (tfIdx < TFCfiles.size() - 1) {
            outf_center << "\t";
        }
    }
    outf_center << "\n";

    for (size_t i = 0; i < totalCellNum; i++) {
        std::filesystem::path p(cellBedFiles[i]);
        string cellName = p.filename().string();
        outf_center << cellName << "\t";

        for (size_t tfIdx = 0; tfIdx < TFCfiles.size(); tfIdx++) {
            outf_center << center_countSum_mat[i][tfIdx];
            if (tfIdx < TFCfiles.size() - 1) {
                outf_center << "\t";
            }
        }
        outf_center << "\n";
    }

    outf_center.close();


    // output flank counts
    string flankfile = prefix + ".flankNuclCountSum.txt";
    ofstream outf_flank(flankfile);
    if (!outf_flank.good()) {
        cerr << "Cannot open output file: " << flankfile << endl;
        return 1;
    }

    for (size_t tfIdx = 0; tfIdx < TFCfiles.size(); tfIdx++) {
        const string &TFC = TFCfiles[tfIdx];
        size_t pos1 = TFC.find_last_of('/');
        string filename = (pos1 == string::npos) ? TFC : TFC.substr(pos1 + 1);
        size_t pos2 = filename.find_last_of('.');
        string centerName = (pos2 == string::npos) ? filename : filename.substr(0, pos2);

        outf_flank << centerName;
        if (tfIdx < TFCfiles.size() - 1) {
            outf_flank << "\t";
        }
    }
    outf_flank << "\n";

    for (size_t i = 0; i < totalCellNum; i++) {
        std::filesystem::path p(cellBedFiles[i]);
        string cellName = p.filename().string();
        outf_flank << cellName << "\t";

        for (size_t tfIdx = 0; tfIdx < TFCfiles.size(); tfIdx++) {
            outf_flank << flank_countSum_mat[i][tfIdx];
            if (tfIdx < TFCfiles.size() - 1) {
                outf_flank << "\t";
            }
        }
        outf_flank << "\n";
    }

    outf_flank.close();
    

    // output summit distance matrix
    string distfile = prefix + ".summitDist.txt";
    ofstream outf(distfile);
    if (!outf.good()) {
        cerr << "Cannot open output file: " << distfile << endl;
        return 1;
    }

    for (size_t tfIdx = 0; tfIdx < TFCfiles.size(); tfIdx++) {
        const string &TFC = TFCfiles[tfIdx];
        size_t pos1 = TFC.find_last_of('/');
        string filename = (pos1 == string::npos) ? TFC : TFC.substr(pos1 + 1);
        size_t pos2 = filename.find_last_of('.');
        string centerName = (pos2 == string::npos) ? filename : filename.substr(0, pos2);

        outf << centerName;
        if (tfIdx < TFCfiles.size() - 1) {
            outf << "\t";
        }
    }
    outf << "\n";

    for (size_t i = 0; i < totalCellNum; i++) {
        std::filesystem::path p(cellBedFiles[i]);
        string cellName = p.filename().string();
        outf << cellName << "\t";

        for (size_t tfIdx = 0; tfIdx < TFCfiles.size(); tfIdx++) {
            outf << summitDist_mat[i][tfIdx];
            if (tfIdx < TFCfiles.size() - 1) {
                outf << "\t";
            }
        }
        outf << "\n";
    }

    outf.close();

    if (remove("tmp2.txt") != 0 && errno != ENOENT) {
        cerr << "Warning: Failed to remove tmp2.txt" << endl;
    }
    if (remove("files.txt") != 0 && errno != ENOENT) {
        cerr << "Warning: Failed to remove files.txt" << endl;
    }
    if (remove("tfc_files.txt") != 0 && errno != ENOENT) {
        cerr << "Warning: Failed to remove tfc_files.txt" << endl;
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_seconds = end - start;
    cout << "Total execution time: " << elapsed_seconds.count() << " seconds" << endl;

    return 0;
}
