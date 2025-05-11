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
void assignTagPos_multiThread(int tagpos, 
                  vector<pair<int, int>> &cut_reg, 
                  vector<int> &counts, 
                  int* colSum, 
                  bool* allHitArr, 
                  double &allHits, 
                  bool* centerHitArr, 
                  double &centerHits, 
                  int sdStartIdx, 
                  int sdEndIdx, 
                  int mymaxReadNum) {
    for (size_t i = 0; i < cut_reg.size(); ++i) {
        if (cut_reg[i].second < tagpos) {  // region's upper limit is less than tag's position. skip it...
            continue;
        }
        if (cut_reg[i].first > tagpos) {   // region's lower limit is beyond tag's position. finish scanning...
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
            allHits += 1;
            if (i >= sdStartIdx && i < sdEndIdx) {
                centerHitArr[i] = true;
                centerHits += 1;
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
                     int colLen) {
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
                    vector<pair<int, int>> localCutRegions = cpool.cut_reg_ve[id];
                    vector<int> tempCounts(colLen, 0);
                    assignTagPos_multiThread(pos, localCutRegions, tempCounts, colSum, allHitArr, allHits, 
                                            centerHitArr, centerHits, sdStartIdx, sdEndIdx, maxReadNum);
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



vector<double> processCellBedFile(const string& cellBed,
                                  const vector<count_pool>& cpool_vec,
                                  int colLen,
                                  int ups,
                                  int downs,
                                  int win,
                                  int motifShift,
                                  double covThre,
                                  int sdReionReadHitThre,
                                  int maxReadNum) {
    map<string, vector<int>> localTagposmap;
    vector<double> localSummitDist_vec;
    
    cout << "# Processing bed file: " << cellBed << endl;
    auto eachCellStart = chrono::high_resolution_clock::now();
    
    readAlignTagFile_multiThread(cellBed, localTagposmap);
    
    for (size_t k = 0; k < cpool_vec.size(); k++) {
        const count_pool &curCpool = cpool_vec[k];
        
        vector<int> colSum(colLen, 0);
        bool* allHitArr = new bool[colLen]();      // allocate bool array and init to false
        double allHits = 0.0;
        int centerIdx = colLen / 2;
        int startIdx = max(0, centerIdx - motifShift);
        int endIdx = min(colLen, centerIdx + motifShift + 1);
        bool* centerHitArr = new bool[colLen](); 
        double centerHits = 0.0;
        
        assignCountPool_multiThread(localTagposmap, curCpool, colSum.data(), allHitArr, allHits,
                                    centerHitArr, centerHits, startIdx, endIdx, maxReadNum, colLen);
        
        double nonZeroPct = allHits / colLen;
        double tmpSummitDist = 0.0;
        if (nonZeroPct >= covThre && centerHits > sdReionReadHitThre) {
            tmpSummitDist = get_summitDist_multiThread(colSum.data(), colLen, motifShift);
        }
        localSummitDist_vec.push_back(tmpSummitDist);
        
        delete[] allHitArr;
        delete[] centerHitArr;
    }
    
    localTagposmap.clear();
    
    auto eachCellEnd = chrono::high_resolution_clock::now();
    chrono::duration<double> eachCellDuration = eachCellEnd - eachCellStart;
    cout << "current cell processing duration: " << eachCellDuration.count() << " seconds" << endl;
    
    return localSummitDist_vec;
}


int main(int argc, char* const argv[]){
    auto start = chrono::high_resolution_clock::now();
    
    const char* short_options = "f:c:o:u:d:w:m:t:n:p:"; 
    static struct option long_options[] = {
        {"fragment", required_argument, nullptr, 'f'},
        {"center", required_argument, nullptr, 'c'},
        {"out_prefix", required_argument, nullptr, 'o'},
        {"up", required_argument, nullptr, 'u'},
        {"down", required_argument, nullptr, 'd'},
        {"window_size", required_argument, nullptr, 'w'},
        {"motif_shift", required_argument, nullptr, 'm'},
        {"coverage_threshold", required_argument, nullptr, 't'},
        {"region_hit_threshold", required_argument, nullptr, 'n'},
        {"threads", required_argument, nullptr, 'p'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    string metaFile, regionPath, prefix = "output";
    int ups = 400, downs = 400, win = 1, motifShift = 200, sdReionReadHitThre = 1;
    double covThre = 0.05;
    int threadNum = 1; 

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
                ups = std::stoi(optarg);
                break;
            case 'd':
                downs = std::stoi(optarg);
                break;
            case 'w':
                win = std::stoi(optarg);
                break;
            case 'm':
                motifShift = std::stoi(optarg);
                break;
            case 't':
                covThre = std::stod(optarg);
                break;
            case 'n':
                sdReionReadHitThre = std::stoi(optarg);
                break;
            case 'p':
                threadNum = std::stoi(optarg);
                break;
            case '?':
                cerr << "Unknown option or missing argument." << endl;
                return 1;
            default:
                break;
        }
    }
    
    if (metaFile.empty() || regionPath.empty()){
        cerr << "Error: Missing required options." << endl;
        cerr << "Usage: " << argv[0] << " <options>" << endl;
        cerr << "-f:        A meta file with the first column containing single-cell BED file paths," << endl;
        cerr << "           and the second column containing labels" << endl;
        cerr << "-c:        TFBS_center" << endl;
        cerr << "-o:        out_prefix" << endl;
        cerr << "-u:        center upstream, default: 400" << endl;
        cerr << "-d:        center downstream, default: 400" << endl;
        cerr << "-w:        window_size, default: 1" << endl;
        cerr << "-m:        motif_shift_size, default: 200" << endl;
        cerr << "-t:        coverage_threshold, default: 0.05" << endl;
        cerr << "-n:        region_hit_threshold, default: 1" << endl;
        cerr << "-p:        number of threads, default: 1 (single-threaded)" << endl;
        return 1;
    }
    
    if (threadNum < 1) {
        cerr << "Error: Number of threads (-p) must be at least 1." << endl;
        return 1;
    }
    
    std::filesystem::path parentDir = std::filesystem::path(prefix).parent_path();
    if (!parentDir.empty() && !std::filesystem::exists(parentDir)){
        std::filesystem::create_directories(parentDir);
    }
    
    vector<string> cellBedFiles = getPathFromMetafile(metaFile);
    vector<string> TFCfiles = getFiles(regionPath, "", "tfc");
    
    size_t totalCellNum = cellBedFiles.size();
    vector<vector<double>> summitDist_mat(totalCellNum);
    
    vector<count_pool> cpool_vec;
    vector<string> centerName_vec;
    cpool_vec.reserve(TFCfiles.size());
    centerName_vec.reserve(TFCfiles.size());
    
    cout << "# Reading TF center: " << endl;
    for (const auto& TFC : TFCfiles) {
        count_pool cpool;
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
            getTFCcpool(tfcenters, cpool);
            cpool.getregionidmap_tfc(ups, downs, win);
            cpool_vec.push_back(std::move(cpool));
        }
    }
    
    size_t validCenterNum = centerName_vec.size();
    
    int centerNum = 0;
    int colLen = 0;
    if (!cpool_vec.empty()){
        count_pool advanceCpool = cpool_vec[0];
        advanceCpool.getregionidmap_tfc(ups, downs, win);
        centerNum = advanceCpool.counts_table.size();
        colLen = centerNum > 0 ? advanceCpool.counts_table[0].size() : 0;
    }
    
    size_t estimatedMemoryPerThread = colLen * (sizeof(int) + 2 * sizeof(bool)) + 1024;
    threadNum = adjustThreadNum(threadNum, estimatedMemoryPerThread);
    cout << "Using " << threadNum << " threads." << endl;
    
    int maxReadNum = 5;
    
    if (threadNum == 1) {
        double progress = 1.0;
        for (size_t i = 0; i < totalCellNum; i++){
            summitDist_mat[i] = processCellBedFile(cellBedFiles[i], cpool_vec,
                                                     colLen, ups, downs, win, motifShift,
                                                     covThre, sdReionReadHitThre, maxReadNum);
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
                summitDist_mat[i] = processCellBedFile(cellBedFiles[i], cpool_vec,
                                                         colLen, ups, downs, win, motifShift,
                                                         covThre, sdReionReadHitThre, maxReadNum);
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
    
    string distfile = prefix + ".summitDist.txt";
    ofstream outf(distfile);
    for (size_t i = 0; i < validCenterNum; i++){
        outf << centerName_vec[i] << "\t";
    }
    outf << "\n";
    
    if (!summitDist_mat.empty()){
        for (size_t i = 0; i < summitDist_mat.size(); i++){
            std::filesystem::path p(cellBedFiles[i]);
            string tmp_basename = p.filename().string();
            outf << tmp_basename << "\t";
            for (size_t j = 0; j < summitDist_mat[i].size(); j++){
                outf << summitDist_mat[i][j] << "\t";
            }
            outf << "\n";
        }
    } else {
        cout << "summit distance matrix is empty, nothing to output." << endl;
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
