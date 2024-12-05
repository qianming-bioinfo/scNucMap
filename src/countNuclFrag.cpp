#include "myoperation.h"
#include <filesystem>
#include <getopt.h>



// similar to "noParall_procSDMat.cpp": 
// the only difference is the metric calculating part which replaces summitDist with nuclOccuScore

int main(int argc, char* const argv[]){
    
    auto start = chrono::high_resolution_clock::now();
    
    // ################ set options ################
	// const char* short_options = "f:c:o:u:d:l:r:w:t:p:";
    const char* short_options = "f:c:o:u:d:l:r:w:t:";
    static struct option long_options[] = {
        {"fragment", required_argument, nullptr, 'f'},
        {"center", required_argument, nullptr, 'c'},
        {"out_prefix", required_argument, nullptr, 'o'},
        {"center_up", required_argument, nullptr, 'u'},
        {"center_down", required_argument, nullptr, 'd'},
        {"flank_up", required_argument, nullptr, 'l'},
        {"flank_down", required_argument, nullptr, 'r'},
        {"window_size", required_argument, nullptr, 'w'},
        {"out_type", required_argument, nullptr, 't'},
        // {"thread", required_argument, nullptr, 'p'},
        {nullptr, 0, nullptr, 0}
    };

    int opt, threadNum = 1;
    // parameters
	string metaFile, regionPath, prefix = "output";
    int center_ups = 100, center_downs = 100, flank_ups = 200, flnak_downs = 200, win = 1, out_type = 1;

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
				flnak_downs = std::stoi(optarg);
                break;
            case 'w':
                win = std::stoi(optarg);
                break;
            case 't':
                out_type = std::stoi(optarg);
                break;
            // case 'p':
            //     threadNum = std::stoi(optarg);
            //     break;
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
        cerr << "-t:        output type" << endl;
        cerr << "           1: output nucleosome read fraction" << endl;
        cerr << "           2: output 2 matrices of the sum of nuclosome read number" << endl;
        cerr << "           3: only output raw count tables around centers" << endl;
        cerr << "           4: only the bool value for each region around centers" << endl;
        cerr << "           5: the hit of read count of each bin on each region around centers" << endl;
        // cerr << "-p:        thread" << endl;
        return 1;
    }

    // Create parent directory if it doesn't exist
    std::filesystem::path parentDir = std::filesystem::path(prefix).parent_path();
    if (!std::filesystem::exists(parentDir)) {
        std::filesystem::create_directories(parentDir);
    }

	// ################ set options ################


    // ## declaration of vital variables and objects not mentioned above
    vector<string> cellBedFiles = getPathFromMetafile(metaFile);
    vector<string> TFCfiles = getFiles(regionPath, ".txt", "tfc");
    
    double i = 1.0;
    size_t totalCellNum = cellBedFiles.size();
    size_t totalCenterNum = TFCfiles.size();

    // map<string, vector<int>> tagposmap;
    vector<double> centerMatSum_vec, flankMatSum_vec, nuclFraction_vec;
    vector<vector<double>> centerMatSum_mat, flankMatSum_mat, nuclFraction_mat; // matrix of totalCellNum rows and totalCenterNum columns

    // allocate memory
    nuclFraction_mat.reserve(totalCellNum);

    // only read TFCfiles once
    // count_pool cpool;
    vector<count_pool> cpool_center_vec, cpool_flank_vec;
    vector<string> centerName_vec;

    cpool_center_vec.reserve(totalCenterNum);
    cpool_flank_vec.reserve(totalCenterNum);
    centerName_vec.reserve(totalCenterNum);

    std::cout << "# Reading center: " << endl;

    // ## 1. generate count_pool vector
    for (const auto& TFC : TFCfiles) {
		
        count_pool cpool_center, cpool_flank; // flank means from the near region around the center is included
        // count_pool cpool;

        string tmpTFC = TFC;
        std::cout << "Processing center: " << tmpTFC << endl;

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
            map<string, vector<int>> tfcenters;
            readRegionCenter(tmpTFC, tfcenters);
            // (1) for cpool_center
            getTFCcpool(tfcenters, cpool_center);
            cpool_center.getregionidmap_tfc(center_ups, center_downs, win);
            cpool_center_vec.push_back(cpool_center);

            // (2) for cpool_flank
            getTFCcpool(tfcenters, cpool_flank);
            cpool_flank.getregionidmap_tfc(flank_ups, flnak_downs, win);
            cpool_flank_vec.push_back(cpool_flank);
        }

    }


    // ## 2. get the pooled read raw number in taotal around center of +/- 100 bp (Count center) for each cell
    // & ## 3. get the pooled read raw number in taotal around center of +/- 100 bp (Count flanking) for each cell
    // & ## 4. calculate the nucleosome fraction of nucleosomes around center (avoid 0 as denominator)

    if (out_type == 5) { // output the sum of read count in each region

        map<string, vector<int>> tagposmap;
        for (const auto& cellBed : cellBedFiles) {

            std::filesystem::path p(cellBed);
            std::string tmp_cell_basename = p.filename().string();
            string tmp_cell_dir = prefix + tmp_cell_basename + "/";
            std::filesystem::create_directories(tmp_cell_dir);

            std::cout << "# Processing bed file: " << cellBed << endl;
            auto eachCellStart = chrono::high_resolution_clock::now();

            readAlignTagFile(cellBed, tagposmap);
            
            for (int k = 0; k < cpool_center_vec.size(); k++) { // loop through each file 
                
                count_pool tmp_cpool_center = cpool_center_vec[k];

                assignCountPool_occuScore(tagposmap, tmp_cpool_center);

                vector<vector<int>> tmp_center_rawcount = tmp_cpool_center.counts_table;

                // sum of hit count for each region
                vector<int> tmp_hit_sum_vec;
                for (const auto& row : tmp_center_rawcount) {
                    int tmp_sum = std::accumulate(row.begin(), row.end(), 0, [](int acc, int val) {
                        return acc + (val > 0);
                    });
                    tmp_hit_sum_vec.push_back(tmp_sum);
                }

                // separated directory for each cell
                string tmp_table_prefix = tmp_cell_basename + "_" + centerName_vec[k];
                string tmp_save_pre = tmp_cell_dir + tmp_table_prefix;
                writeVec(tmp_save_pre, tmp_hit_sum_vec);

            }

            auto eachCellEnd = chrono::high_resolution_clock::now();
            std::chrono::duration<double> eachCellDuration = eachCellEnd - eachCellStart;
            std::cout << "current cell processing duration: "  << eachCellDuration.count() << " seconds" << endl;

            std::cout << setprecision(4) << (i/totalCellNum)*100 << "%" << endl;
            i += 1;

        }

    }


    if (out_type == 4) { // only out put the bool value of each region which stands for the presence of nucleosome

        map<string, vector<int>> tagposmap;
        for (const auto& cellBed : cellBedFiles) {

            std::filesystem::path p(cellBed);
            std::string tmp_cell_basename = p.filename().string();
            string tmp_cell_dir = prefix + tmp_cell_basename + "/";
            std::filesystem::create_directories(tmp_cell_dir);

            std::cout << "# Processing bed file: " << cellBed << endl;
            auto eachCellStart = chrono::high_resolution_clock::now();

            readAlignTagFile(cellBed, tagposmap);

            for (int k = 0; k < cpool_center_vec.size(); k++) { 
                count_pool tmp_cpool_center = cpool_center_vec[k];

                assignCountPool_occuScore(tagposmap, tmp_cpool_center);

                // 处理并生成布尔值并转为int的矩阵
                vector<vector<int>> tmp_center_rawcount = tmp_cpool_center.counts_table;
                vector<vector<int>> tmp_center_bool(tmp_center_rawcount.size(), vector<int>(tmp_center_rawcount[0].size()));

                // 使用std::transform来优化行的遍历，并将布尔值直接转换为int
                std::transform(tmp_center_rawcount.begin(), tmp_center_rawcount.end(), tmp_center_bool.begin(),
                    [](const vector<int>& row) {
                        return std::vector<int>(1, std::accumulate(row.begin(), row.end(), 0) > 0 ? 1 : 0);
                    });

                // 输出转换为int的布尔值表格
                string tmp_table_prefix = tmp_cell_basename + "_" + centerName_vec[k];
                string tmp_save_pre = tmp_cell_dir + tmp_table_prefix;
                writeRawCountTable(tmp_save_pre, tmp_center_bool);  // tmp_center_bool是int类型的矩阵
            }
        }

    }


    if (out_type == 3 && threadNum == 1){
        
        map<string, vector<int>> tagposmap;
        for (const auto& cellBed : cellBedFiles) {

            std::filesystem::path p(cellBed);
            std::string tmp_cell_basename = p.filename().string();
            string tmp_cell_dir = prefix + tmp_cell_basename + "/";
            std::filesystem::create_directories(tmp_cell_dir);

            std::cout << "# Processing bed file: " << cellBed << endl;
            auto eachCellStart = chrono::high_resolution_clock::now();

            readAlignTagFile(cellBed, tagposmap);
            
            for (int k = 0; k < cpool_center_vec.size(); k++) { // loop through each file 
                
                count_pool tmp_cpool_center = cpool_center_vec[k];
                // count_pool tmp_cpool_flank = cpool_flank_vec[k];

                assignCountPool_occuScore(tagposmap, tmp_cpool_center);
                // assignCountPool_occuScore(tagposmap, tmp_cpool_flank);

                // # output raw count table for hypergeometric test
                vector<vector<int>> tmp_center_rawcount;
                tmp_center_rawcount = tmp_cpool_center.counts_table;
                // separated directory for each cell
                string tmp_table_prefix = tmp_cell_basename + "_" + centerName_vec[k];
                string tmp_save_pre = tmp_cell_dir + tmp_table_prefix;
                writeRawCountTable(tmp_save_pre, tmp_center_rawcount);

            }

            auto eachCellEnd = chrono::high_resolution_clock::now();
            std::chrono::duration<double> eachCellDuration = eachCellEnd - eachCellStart;
            std::cout << "current cell processing duration: "  << eachCellDuration.count() << " seconds" << endl;

            std::cout << setprecision(4) << (i/totalCellNum)*100 << "%" << endl;
            i += 1;

        }
    }


    // muilti-threading
    if (out_type == 3 && threadNum > 1) {

        // get the number of threads supported by the hardware
        unsigned int max_threads = std::thread::hardware_concurrency();

        if(threadNum > max_threads) {
            cout << "## warning: the number of threads is larger than the number of threads supported by the hardware." << endl;
            cout << "## thread number will be set as: " << max_threads << endl;
            threadNum = std::min(max_threads, static_cast<unsigned int>(threadNum));
        }
        else {
            cout << "## thread number: " << threadNum << endl;
        }


        ThreadPool pool(threadNum);
        vector<future<void>> futures;

        size_t processedFiles = 0;

        for (const auto& cellBed : cellBedFiles) {
            futures.push_back(pool.enqueue(processCell, cellBed, std::cref(cpool_center_vec), prefix, centerName_vec));
        }

        // Wait for all threads to finish
        size_t totalFiles = futures.size();
        for (size_t i = 0; i < totalFiles; ++i) {
            futures[i].get();
            cout << "Processed file " << (i + 1) << " of " << totalFiles << " (" << setprecision(4) << static_cast<int>((static_cast<double>(i + 1) / totalFiles) * 100) << "% complete)" << endl;
        }

        // // Wait for all threads to finish
        // for (auto& future : futures) {
        //     future.get();
        // }
    }


    else{ // raw count tables will not be output

        map<string, vector<int>> tagposmap;
        for (const auto& cellBed : cellBedFiles) {
            
            std::cout << "# Processing bed file: " << cellBed << endl;
            auto eachCellStart = chrono::high_resolution_clock::now();

            readAlignTagFile(cellBed, tagposmap);
            
            for (int k = 0; k < cpool_center_vec.size(); k++) { // loop through each file 
                
                count_pool tmp_cpool_center = cpool_center_vec[k];
                count_pool tmp_cpool_flank = cpool_flank_vec[k];

                assignCountPool_occuScore(tagposmap, tmp_cpool_center);
                assignCountPool_occuScore(tagposmap, tmp_cpool_flank);

                double centerMatSum = calculateCpoolMatSum(tmp_cpool_center);
                double flankMatSum = calculateCpoolMatSum(tmp_cpool_flank);
                double tmpNuclFraction = centerMatSum / flankMatSum;
                // double tmpNuclFraction = nuclFraction(tmp_cpool_center, tmp_cpool_flank);

                centerMatSum_vec.push_back(centerMatSum);
                flankMatSum_vec.push_back(flankMatSum);
                nuclFraction_vec.push_back(tmpNuclFraction);
                
            }

            centerMatSum_mat.push_back(centerMatSum_vec);
            flankMatSum_mat.push_back(flankMatSum_vec);
            nuclFraction_mat.push_back(nuclFraction_vec);

            // map<string, vector<int>>().swap(tagposmap);
            // vector<double>().swap(nuclFraction_vec);
            tagposmap.clear();
            centerMatSum_vec.clear();
            flankMatSum_vec.clear();
            nuclFraction_vec.clear();

            auto eachCellEnd = chrono::high_resolution_clock::now();
            std::chrono::duration<double> eachCellDuration = eachCellEnd - eachCellStart;
            std::cout << "current cell processing duration: "  << eachCellDuration.count() << " seconds" << endl;

            std::cout << setprecision(4) << (i/totalCellNum)*100 << "%" << endl;
            i += 1;

        }
    }

    // 5. compare with the background for each cell as the final score
    // ************* some code *************

    // ## output
    if (out_type == 1){

        string nuclFracFile = prefix + ".nuclFraction.txt";
        ofstream outf(nuclFracFile.data());

        for (int i = 0; i < totalCenterNum; i++) {
            outf << centerName_vec[i] << "\t";
        }
        outf << "\n";


        if (!nuclFraction_mat.empty()){
            for (int i = 0; i < nuclFraction_mat.size(); i++){
                std::filesystem::path p(cellBedFiles[i]);
                std::string tmp_basename = p.filename().string();
                outf << tmp_basename << "\t";
                for (int j = 0; j < nuclFraction_mat[i].size(); j++){
                    outf << nuclFraction_mat[i][j] << "\t";
                }
                outf << "\n";
            }
        } 
        else{
            std::cout << "nucleosome fraction matrix is empty, nothing to output." << endl;
        }


        outf.close();

        remove("tmp2.txt");
        remove("files.txt");

        auto end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "Total execution time: "  << elapsed_seconds.count() << " seconds" << endl;
        
    }
    else if (out_type == 2){
        
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


        if (!centerMatSum_mat.empty()){
            for (int i = 0; i < centerMatSum_mat.size(); i++){
                std::filesystem::path p(cellBedFiles[i]);
                std::string tmp_basename = p.filename().string();
                outf_center << tmp_basename << "\t";
                for (int j = 0; j < centerMatSum_mat[i].size(); j++){
                    outf_center << centerMatSum_mat[i][j] << "\t";
                }
                outf_center << "\n";
            }
        } 
        else{
            std::cout << "center sum matrix is empty, nothing to output." << endl;
        }

        if (!flankMatSum_mat.empty()){
            for (int i = 0; i < flankMatSum_mat.size(); i++){
                std::filesystem::path p(cellBedFiles[i]);
                std::string tmp_basename = p.filename().string();
                outf_flank << tmp_basename << "\t";
                for (int j = 0; j < flankMatSum_mat[i].size(); j++){
                    outf_flank << flankMatSum_mat[i][j] << "\t";
                }
                outf_flank << "\n";
            }
        } 
        else{
            std::cout << "flank sum matrix is empty, nothing to output." << endl;
        }


        outf_center.close();
        outf_flank.close();

        remove("tmp2.txt");
        remove("files.txt");

        auto end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "Total execution time: "  << elapsed_seconds.count() << " seconds" << endl;

    }
    else if (out_type == 3){

        auto end = chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "Total execution time: "  << elapsed_seconds.count() << " seconds" << endl;

    }
    else{

        std::cout << "Error: unknown output type." << endl;
    
    }

    return 0;

}



    

