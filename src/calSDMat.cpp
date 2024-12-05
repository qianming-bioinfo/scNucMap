#include "myoperation.h"
#include <filesystem>
#include <getopt.h>

int main(int argc, char* const argv[]){

    auto start = chrono::high_resolution_clock::now();
    
    // ################ set options ################
	const char* short_options = "f:c:o:u:d:w:m:t:n:";
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
        {nullptr, 0, nullptr, 0}
    };

    int opt, threadNum = 1;
    // parameters
	string metaFile, regionPath, prefix = "output";
    int ups = 400, downs = 400, win = 1, motifShift = 200, sdReionReadHitThre = 1;
    double covThre = 0.05;

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
        return 1;
    }

    // Create parent directory if it doesn't exist
    std::filesystem::path parentDir = std::filesystem::path(prefix).parent_path();
    if (!std::filesystem::exists(parentDir)) {
        std::filesystem::create_directories(parentDir);
    }

	// ################ set options ################

    // directory as the input
    // vector<string> cellBedFiles = getFiles(metaFile, ".bed", "cell"); // get cell bed files' path set
    vector<string> cellBedFiles = getPathFromMetafile(metaFile);
    vector<string> TFCfiles = getFiles(regionPath, ".txt", "tfc"); // get TF center files' path set
    
    double i = 1.0;
    size_t totalCellNum = cellBedFiles.size();
    size_t totalCenterNum = TFCfiles.size();

    map<string, vector<int>> tagposmap;
    vector<double> summitDist_vec;
    vector<vector<double>> summitDist_mat; // matrix of totalCellNum rows and totalCenterNum columns


    // allocate memory
    summitDist_mat.reserve(totalCellNum);

    // only read TFCfiles once
    // count_pool cpool;
    vector<count_pool> cpool_vec;
    vector<string> centerName_vec;

    cpool_vec.reserve(totalCenterNum);
    centerName_vec.reserve(totalCenterNum);

    std::cout << "# Reading TF center: " << endl;


    for (const auto& TFC : TFCfiles){
		
        count_pool cpool;
        // count_pool cpool;

        string tmpTFC = TFC;
        std::cout << "Processing TF center: " << tmpTFC << endl;

        size_t pos1 = tmpTFC.find_last_of('/'); // find the position of the last '/'
        size_t pos2 = tmpTFC.find_last_of('.'); // find the position of the last '.'
        string centerName = tmpTFC.substr(pos1 + 1, pos2 - pos1 - 1);
        centerName_vec.push_back(centerName);

        ifstream myfilefile(tmpTFC);

        if (isFileEmpty(myfilefile)) {
            std::cout << "### WARNING: File " << tmpTFC << " is empty which will cause an error!" << endl;
            std::cout << "skip this file..." << endl;
            centerName_vec.pop_back();
            // when handling empty files
        }
        else {
            // map<string, vector<int>> tfcenters;
            map<string, vector<int>> tfcenters;
            readRegionCenter(tmpTFC, tfcenters);
            getTFCcpool(tfcenters, cpool);
            cpool.getregionidmap_tfc(ups, downs, win);
            cpool_vec.push_back(cpool);
        }

    }


    // get centerNum and colLen in advance
    int centerNum = 0;
    int colLen = 0;
    if (!cpool_vec.empty()){
        count_pool advanceCpool = cpool_vec[0];
        advanceCpool.getregionidmap_tfc(ups, downs, win);
        centerNum = advanceCpool.counts_table.size();
        colLen = centerNum > 0 ? advanceCpool.counts_table[0].size() : 0;
    }



    for (const auto& cellBed : cellBedFiles) {

        std::cout << "# Processing bed file: " << cellBed << endl;
        auto eachCellStart = chrono::high_resolution_clock::now();

        readAlignTagFile(cellBed, tagposmap);
        // auto readAlignTagFile_end = chrono::high_resolution_clock::now();
        // std::chrono::duration<double> readAlignTagFile_Duration = readAlignTagFile_end - eachCellStart;
        // std::cout << "reading file takes: "  << readAlignTagFile_Duration.count() << " seconds" << endl;
        
        for (int k = 0; k < cpool_vec.size(); k++) { // loop through each file 
            
            count_pool tmpcpool = cpool_vec[k];
                    
            // ####################################################################################
            // ######################################## QC ########################################
            
            int colSum[colLen] = {0};

            bool allHitArr[colLen] ={};
            double allHits = 0.0;
            int centerIdx = colLen / 2;
            int startIdx = max(0, centerIdx - motifShift);
            int endIdx = min(colLen, centerIdx + motifShift + 1);
            bool centerHitArr[colLen] ={};
            double centerHits = 0.0;

            // ######################################## QC ########################################
            // ####################################################################################

            assignCountPool(tagposmap, tmpcpool, colSum, allHitArr, allHits, centerHitArr, centerHits, startIdx, endIdx);

            double nonZeroPct = allHits / colLen;
            double tmpSummitDist;

            if (nonZeroPct >= covThre && centerHits > sdReionReadHitThre){
                tmpSummitDist = get_summitDist(colSum, colLen);
            }
            else {
                tmpSummitDist = 0.0;
            }

            // if (!tmpSummitDist > 0){
            //     tmpSummitDist = 0.0;
            //     std::cout << "error when calculating summit distance of:" << endl;
            //     std::cout << "cell: " << cellBed << endl;
            //     std::cout << "and" << endl;
            //     std::cout << "center: " << centerName_vec[k] << endl;
            //     std::cout << "the result has been recorded as 0." << endl;
            // }

            // std::cout << tmpSummitDist << endl;

            summitDist_vec.push_back(tmpSummitDist);
            
        }

        summitDist_mat.push_back(summitDist_vec);

        // map<string, vector<int>>().swap(tagposmap);
        // vector<double>().swap(summitDist_vec);
        tagposmap.clear();
        summitDist_vec.clear();

        auto eachCellEnd = chrono::high_resolution_clock::now();
        std::chrono::duration<double> eachCellDuration = eachCellEnd - eachCellStart;
        std::cout << "current cell processing duration: "  << eachCellDuration.count() << " seconds" << endl;

        std::cout << setprecision(4) << (i/totalCellNum)*100 << "%" << endl;
        i += 1;

    }


    string distfile = prefix + ".summitDist.txt";
    ofstream outf(distfile.data());

    for (int i = 0; i < totalCenterNum; i++) {
        outf << centerName_vec[i] << "\t";
    }
    outf << "\n";


    if (!summitDist_mat.empty()){
        for (int i = 0; i < summitDist_mat.size(); i++){
            std::filesystem::path p(cellBedFiles[i]);
            std::string tmp_basename = p.filename().string();
            outf << tmp_basename << "\t";
            for (int j = 0; j < summitDist_mat[i].size(); j++){
                outf << summitDist_mat[i][j] << "\t";
            }
            outf << "\n";
        }
    } 
    else{
        std::cout << "summit distance matrix is empty, nothing to output." << endl;
    }


    outf.close();

    remove("tmp2.txt");
    remove("files.txt");


    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Total execution time: "  << elapsed_seconds.count() << " seconds" << endl;
    
    return 0;

}
