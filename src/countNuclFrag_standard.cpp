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
    int maxThreads = availMem / perThreadEstimate;
    return std::max(1, std::min(desiredThreads, maxThreads));
}

// ============================

void assignTagPos_occuScore_multiThread(int tagpos, 
                                        vector<pair<int, int>> &cut_reg, 
                                        int* colSum) {

    for (size_t i = 0; i < cut_reg.size(); ++i) {
        if (cut_reg[i].second < tagpos) {  // region's upper limit is less than tag's position. skip it...
            continue;
        }
        if (cut_reg[i].first > tagpos) {   // region's lower limit is beyond tag's position. finish scanning...
            break;
        }
        colSum[i] += 1;
        break; 
    }

}


void assignCountPool_occuScore_multiThread(const map<string, vector<int>> &tagpos, 
                                           const count_pool &cpool, 
                                           int* colSum){

    int L = 100000;
    map<string, map<int, set<pair<int, int>>>> chr_index_region;
    for (const auto &ite : cpool.region_strand_map) {
        string chr = ite.first;
        for (const auto &si : ite.second) {
            int start = si.first.first;
            int end = si.first.second;
            int index1 = start / L;
            int index2 = end / L;
            chr_index_region[chr][index1].insert(si.first);
            if (index2 > index1) {
                for (int i = index1 + 1; i <= index2; ++i) {
                    chr_index_region[chr][i].insert(si.first);
                }
            }
        }
    }
    
    for (const auto &ite : tagpos) {
        string chr = ite.first;
        for (const auto &pos : ite.second) {
            int index = pos / L;
            if (chr_index_region.find(chr) != chr_index_region.end() &&
                chr_index_region[chr].find(index) != chr_index_region[chr].end()) {
                for (auto ci = chr_index_region[chr][index].begin(); ci != chr_index_region[chr][index].end(); ++ci) {
                    if (ci->second < pos)
                        continue;
                    if (ci->first > pos)
                        break;
                    size_t id = cpool.region_id_map.at(chr).at(*ci);
                    // copy for memory saving
                    vector<pair<int, int>> localCutRegions = cpool.cut_reg_ve[id];
                    assignTagPos_occuScore_multiThread(pos, localCutRegions, colSum);
				}
			}

		}
	}
}




void readAlignTagFile_multiThread(string infile, map<string, vector<int>> &tagpos){
    ifstream inf(infile.data());
    if (!inf.good()){
        cout << "file " << infile << " not found " << endl; 
        exit(1);
    }
    while(!inf.eof()){
        string line;
        getline(inf, line);
        if (line.empty()){
            break;
        }
        vector<string> parseditem = parseString(line);
        if (parseditem.size() != 6){
            cout << "error line: " << line << endl; 
            cout << "By default, you need to use BED6 format with the 6th column indicating strand." << endl;
            exit(1);
        }
        string chr = parseditem[0];
        int start = atoi(parseditem[1].c_str());
        int end = atoi(parseditem[2].c_str());
        char strand = parseditem[5][0];
        if (strand != '+' && strand != '-'){
            cout << "error strand " << strand << endl;
            exit(1);
        }
        int shiftpos = start + (end - start) / 2;
        tagpos[chr].push_back(shiftpos);
    }
    inf.close();
    for(auto &p : tagpos) {
        sort(p.second.begin(), p.second.end());
    }
}




pair<vector<int>, vector<int>> processCellBedFile_readCount(const string& cellBed,
                                                            const vector<count_pool>& cpool_center_vec,
                                                            const vector<count_pool>& cpool_flank_vec,
                                                            int colLen_center,
                                                            int colLen_flank) {

    map<string, vector<int>> localTagposmap;
    vector<int> centerMatSum_vec, flankMatSum_vec;
    
    cout << "# Processing bed file: " << cellBed << endl;
    auto eachCellStart = chrono::high_resolution_clock::now();
    
    readAlignTagFile_multiThread(cellBed, localTagposmap);
    
    for (size_t k = 0; k < cpool_center_vec.size(); k++) {

        const count_pool &tmp_cpool_center = cpool_center_vec[k];
        const count_pool &tmp_cpool_flank = cpool_flank_vec[k];
        
        vector<int> colSum_center(colLen_center, 0);
        vector<int> colSum_flank(colLen_flank, 0);
        
        assignCountPool_occuScore_multiThread(localTagposmap, tmp_cpool_center, colSum_center.data());
        assignCountPool_occuScore_multiThread(localTagposmap, tmp_cpool_flank, colSum_flank.data());

        int centerMatSum = accumulate(colSum_center.begin(), colSum_center.end(), 0);
        int flankMatSum = accumulate(colSum_flank.begin(), colSum_flank.end(), 0);

        centerMatSum_vec.push_back(centerMatSum);
        flankMatSum_vec.push_back(flankMatSum);
    }
    
    localTagposmap.clear();
    
    auto eachCellEnd = chrono::high_resolution_clock::now();
    chrono::duration<double> eachCellDuration = eachCellEnd - eachCellStart;
    cout << "current cell processing duration: " << eachCellDuration.count() << " seconds" << endl;

    pair<vector<int>, vector<int>> result = make_pair(centerMatSum_vec, flankMatSum_vec);
    return result;
}

// ============================
int main(int argc, char* const argv[]){
    
    auto start = chrono::high_resolution_clock::now();
    
    // ################ set options ################
	// const char* short_options = "f:c:o:u:d:l:r:w:t:p:";
    const char* short_options = "f:c:o:u:d:l:r:w:p:";
    static struct option long_options[] = {
        {"fragment", required_argument, nullptr, 'f'},
        {"center", required_argument, nullptr, 'c'},
        {"out_prefix", required_argument, nullptr, 'o'},
        {"center_up", required_argument, nullptr, 'u'},
        {"center_down", required_argument, nullptr, 'd'},
        {"flank_up", required_argument, nullptr, 'l'},
        {"flank_down", required_argument, nullptr, 'r'},
        {"window_size", required_argument, nullptr, 'w'},
        // {"out_type", required_argument, nullptr, 't'},
        {"thread", required_argument, nullptr, 'p'},
        {nullptr, 0, nullptr, 0}
    };

    int opt, threadNum = 1;
    // parameters
	string metaFile, regionPath, prefix = "output";
    int center_ups = 100, center_downs = 100, flank_ups = 200, flank_downs = 200, win = 1, out_type = 1;

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
				center_ups = std::stoi(optarg);
                break;
            case 'd':
				center_downs = std::stoi(optarg);
                break;
            case 'l':
				flank_ups = std::stoi(optarg);
                break;
            case 'r':
				flank_downs = std::stoi(optarg);
                break;
            case 'w':
                win = std::stoi(optarg);
                break;
            // case 't':
            //     out_type = std::stoi(optarg);
            //     break;
            case 'p':
                threadNum = std::stoi(optarg);
                break;
            case '?':
                // Handle unknown option or missing argument
                std::cerr << "Unknown option or missing argument." << std::endl;
				return 1;
                break;
            default:
                // Handle other options if needed
                break;
        }
    }

	// Check if required options are provided
    if (metaFile.empty() || regionPath.empty()) {
        cerr << "Error: Missing required options." << std::endl;
        cerr << "Usage: " << argv[0] << "<options>" << endl;
        cerr << "-f:        A meta file with the first column containing single-cell BED file paths," << endl;
        cerr << "           and the second column containing labels" << endl;
        cerr << "-c:        TFBS_center" << endl;
        cerr << "-o:        out_prefix" << endl;
        cerr << "-u:        center upstream, default: 100" << endl;
        cerr << "-d:        center downstream, default: 100" << endl;
        cerr << "-l:        flank upstream, default: 200" << endl;
        cerr << "-r:        flank downstream, default: 200" << endl;
        cerr << "-w:        window_size, default: 1" << endl;
        // cerr << "-t:        output type" << endl;
        // cerr << "           1: output nucleosome read fraction" << endl;
        // cerr << "           2: output 2 matrices of the sum of nuclosome read number" << endl;
        // cerr << "           3: only output raw count tables around centers" << endl;
        // cerr << "           4: only the bool value for each region around centers" << endl;
        // cerr << "           5: the hit of read count of each bin on each region around centers" << endl;
        cerr << "-p:        thread, default: 1" << endl;
        return 1;
    }

    // Create parent directory if it doesn't exist
    std::filesystem::path parentDir = std::filesystem::path(prefix).parent_path();
    if (!std::filesystem::exists(parentDir)) {
        std::filesystem::create_directories(parentDir);
    }

	// ################ set options ################
    
    vector<string> cellBedFiles = getPathFromMetafile(metaFile);
    vector<string> TFCfiles = getFiles(regionPath, ".txt", "tfc");
    
    size_t totalCellNum = cellBedFiles.size();
    size_t totalCenterNum = TFCfiles.size();

    vector<vector<int>> centerCount_mat(totalCellNum), flankCount_mat(totalCellNum);
    
    vector<count_pool> cpool_center_vec, cpool_flank_vec;
    vector<string> centerName_vec;

    cpool_center_vec.reserve(TFCfiles.size());
    cpool_flank_vec.reserve(TFCfiles.size());
    centerName_vec.reserve(TFCfiles.size());
    
    cout << "# Reading TF center: " << endl;
    for (const auto& TFC : TFCfiles) {
        count_pool cpool_center, cpool_flank;
        string tmpTFC = TFC;
        cout << "Processing TF center: " << tmpTFC << endl;
        
        size_t pos1 = tmpTFC.find_last_of('/');
        string filename = (pos1 == string::npos) ? tmpTFC : tmpTFC.substr(pos1 + 1);
        size_t pos2 = filename.find_last_of('.');
        string centerName = (pos2 == string::npos) ? filename : filename.substr(0, pos2);
        if(centerName.empty()){
            cout << "### WARNING: Invalid TF center name from file " << tmpTFC << ". Skipping." << endl;
            continue;
        }
        centerName_vec.push_back(centerName);
        
        ifstream myfile(tmpTFC);
        if (!myfile.good() || myfile.peek() == ifstream::traits_type::eof()){
            cout << "### WARNING: File " << tmpTFC << " is empty which will cause an error!" << endl;
            cout << "Skipping this file..." << endl;
            centerName_vec.pop_back();
        } else {
            map<string, vector<int>> tfcenters;
            readRegionCenter(tmpTFC, tfcenters);
            // (1) for cpool_center
            getTFCcpool(tfcenters, cpool_center);
            cpool_center.getregionidmap_tfc(center_ups, center_downs, win);
            cpool_center_vec.push_back(cpool_center);
            // (2) for cpool_flank
            getTFCcpool(tfcenters, cpool_flank);
            cpool_flank.getregionidmap_tfc(flank_ups, flank_downs, win);
            cpool_flank_vec.push_back(cpool_flank);
        }
    }
    
    size_t validCenterNum = centerName_vec.size();
    
    int centerNum = 0;
    int colLen_center = 0, colLen_flank = 0;
    count_pool advanceCpool;
    if (!cpool_center_vec.empty() && !cpool_flank_vec.empty()){
        // for cpool_center_vec
        advanceCpool = cpool_center_vec[0];
        advanceCpool.getregionidmap_tfc(center_ups, center_downs, win);
        centerNum = advanceCpool.counts_table.size();
        colLen_center = centerNum > 0 ? advanceCpool.counts_table[0].size() : 0;
        // for cpool_flank_vec
        advanceCpool = cpool_flank_vec[0];
        advanceCpool.getregionidmap_tfc(flank_ups, flank_downs, win);
        centerNum = advanceCpool.counts_table.size();
        colLen_flank = centerNum > 0 ? advanceCpool.counts_table[0].size() : 0;
    }
    
    size_t estimatedMemoryPerThread = (colLen_center + colLen_flank) * (sizeof(int) + 2 * sizeof(bool)) + 1024;
    threadNum = adjustThreadNum(threadNum, estimatedMemoryPerThread);
    cout << "Using " << threadNum << " threads." << endl;
    
    if (threadNum == 1) {
        double progress = 1.0;
        for (size_t i = 0; i < totalCellNum; i++){
            pair<vector<int>, vector<int>> res_pair = processCellBedFile_readCount(cellBedFiles[i], cpool_center_vec, cpool_flank_vec, colLen_center, colLen_flank);
            centerCount_mat[i] = res_pair.first;
            flankCount_mat[i] = res_pair.second;
            cout << setprecision(4) << (progress / totalCellNum) * 100 << "%" << endl;
            progress += 1.0;
        }
    } else {
        atomic<size_t> nextIndex(0);
        mutex coutMutex;
        
        auto worker = [&]() {
            while (true) {
                size_t i = nextIndex.fetch_add(1);
                if (i >= totalCellNum)
                    break;
                pair<vector<int>, vector<int>> res_pair = processCellBedFile_readCount(cellBedFiles[i], cpool_center_vec, cpool_flank_vec, colLen_center, colLen_flank);
                centerCount_mat[i] = res_pair.first;
                flankCount_mat[i] = res_pair.second;
                {
                    lock_guard<mutex> lock(coutMutex);
                    // cout << setprecision(4) << ((i + 1) / double(totalCellNum)) * 100 << "%" << endl;
                }
            }
        };
        
        vector<future<void>> futures;
        for (int j = 0; j < threadNum; j++) {
            futures.push_back(async(launch::async, worker));
        }
        for (auto &f : futures) {
            f.get();
        }
    }
    
    string centerMatSumFile = prefix + ".centerNuclCountSum.txt";
    string flankMatSumFile = prefix + ".flankNuclCountSum.txt";

    ofstream outf_center(centerMatSumFile.data());
    ofstream outf_flank(flankMatSumFile.data());

    for (int i = 0; i < totalCenterNum; i++) {
        outf_center << centerName_vec[i] << "\t";
    }
    outf_center << "\n";

    for (int i = 0; i < totalCenterNum; i++) {
        outf_flank << centerName_vec[i] << "\t";
    }
    outf_flank << "\n";


    if (!centerCount_mat.empty()){
        for (int i = 0; i < centerCount_mat.size(); i++){
            std::filesystem::path p(cellBedFiles[i]);
            std::string tmp_basename = p.filename().string();
            outf_center << tmp_basename << "\t";
            for (int j = 0; j < centerCount_mat[i].size(); j++){
                outf_center << centerCount_mat[i][j] << "\t";
            }
            outf_center << "\n";
        }
    } 
    else{
        std::cout << "center sum matrix is empty, nothing to output." << endl;
    }

    if (!flankCount_mat.empty()){
        for (int i = 0; i < flankCount_mat.size(); i++){
            std::filesystem::path p(cellBedFiles[i]);
            std::string tmp_basename = p.filename().string();
            outf_flank << tmp_basename << "\t";
            for (int j = 0; j < flankCount_mat[i].size(); j++){
                outf_flank << flankCount_mat[i][j] << "\t";
            }
            outf_flank << "\n";
        }
    } 
    else{
        std::cout << "flank sum matrix is empty, nothing to output." << endl;
    }


    outf_center.close();
    outf_flank.close();
    
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
