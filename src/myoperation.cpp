#include "myoperation.h"



// ####### system #######

vector<string> getFiles(const string& folder, const string& suffix, string fileListName){

    vector<string> files; 
    string command = "ls " + folder + "/*" + suffix + " > " + fileListName + "_files.txt";
    system(command.c_str());
    
    string fileList = fileListName + "_files.txt";
    ifstream fin(fileList);

    if (fin.is_open()){
        string filename; 
        while(getline(fin, filename)){
            files.push_back(filename);
        }
        fin.close(); 
    }
    else {
        cout << "Error: cannot open files matching with input format" << endl;
    }

    return files;

}



bool isFileEmpty(ifstream& file) {

    return file.peek() == ifstream::traits_type::eof();
    
}



vector<string> getPathFromMetafile(const string& filename) {
    vector<string> pathCol;
    ifstream file(filename);

    if (!file.is_open()) {
        cerr << "Could not open meta file!" << endl;
        return pathCol;
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string firstColumn;

        if (getline(ss, firstColumn, '\t')) {
            pathCol.push_back(firstColumn);
        }
    }

    file.close(); 
    return pathCol;
}



vector<vector<string>> splitVector(const vector<string>& vec, int numParts){
    
    vector<vector<string>> result;
    int totalSize = vec.size();
    int baseSize = totalSize / numParts;
    int remainder = totalSize % numParts;

    auto it = vec.begin();
    for (int i = 0; i < numParts; ++i) {
        int currentSize = baseSize + (i < remainder ? 1 : 0);
        vector<string> part(it, it + currentSize);
        result.push_back(part);
        it += currentSize;
    }

    return result;
}

// ####### system #######



// ####### read file #######

vector<string> parseString(string & instr){

    vector<string> bedVec;
    string s = "";
    for (size_t i = 0; i < instr.size(); ++i) {
        if (instr[i] == '\t' || instr[i] == ' ') {
            if (!s.empty()) {
                bedVec.push_back(s);
                s = "";
            }
        } else {
            s += instr[i];
        }
    }
    if (!s.empty()) {
        bedVec.push_back(s);
    }
    return bedVec;

}




void readAlignTagFile(string infile, map<string, vector<int>> &tagpos){

    // std::cout << "start readAlignTagFile" << infile << endl;

	ifstream inf(infile.data());
	if (!inf.good())
	{
		cout << "file " << infile << " not found " << endl; 
        exit(1);
	}

	// cout << "read file " << infile << endl;
	
	while(!inf.eof())
	{
		string line;
		getline(inf, line);
		if (line.empty()){
			break;
		}
				
		// vector<string> parseditem = parse_string(line);
        vector<string> parseditem = parseString(line);

		if (parseditem.size() != 6)
		{
			cout << "error line: " << line << endl; 
			cout << "By default, you need to use BED6 format with the 6th column indicating strand." << endl;
			// cout << "Otherwise, set -f to 0 if you use BED3 format" << endl;
			exit(1);
		}

		string chr = parseditem[0];
		int start = atoi(parseditem[1].c_str());
		int end = atoi(parseditem[2].c_str());
		char strand = parseditem[5][0];
		
        if (strand != '+' && strand != '-')
		{
			cout << "error strand " << strand << endl;
			exit(1);
		}

		int shiftpos = start + (end-start)/2;
		tagpos[chr].push_back(shiftpos);
	}
	
	inf.close();

    // std::cout << "end readAlignTagFile" << infile << endl;

}




void readRegionCenter(string &infile, map<string, vector<int>> &regioncenter){

	ifstream inf(infile.data());

	if (!inf.good()){
		cout << "file " << infile << " not found!" << endl;
		exit(1);
	}

	// cout << "read file " << infile << endl;
	map<string, set<int>> chr_pos_map;
	string line;

	while(!inf.eof() )
	{
		getline(inf, line);
		if (line.empty()){
			break;
		}
		
		vector<string> parseditem = parseString(line);

		string chr = parseditem[0];
		int c = atoi(parseditem[1].c_str());
		chr_pos_map[chr].insert(c);
	}
	
	inf.close();
	
	regioncenter.clear();
	for (map<string, set<int>>::iterator ite = chr_pos_map.begin(); ite != chr_pos_map.end(); ++ite)
	{
		string chr = ite->first;
		for (set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			regioncenter[chr].push_back(*si);
		}
	}
}

// ####### read file #######





// ####### count pool #######

string int2str(int i){
	
	string Res;
	ostringstream convert;
	convert << i;
	Res = convert.str();
	return Res;

}




void getTFCcpool(map<string, vector<int>> &tfcenter, count_pool &cpool){

	for (map<string, vector<int>>::iterator ite = tfcenter.begin(); ite != tfcenter.end(); ++ite){ // on every chromosome
		for (vector<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si){ // for each center

			string sid = ite->first+"_"+int2str(*si); // like "chr1_3333"
			cpool.id_ve.push_back(sid);

			vector<int> ini_ve;
			cpool.counts_table.push_back(ini_ve); // 1 center, 1 line
			
            vector<double> ini_rve;
			cpool.tss_strand_map[ite->first][*si] = '+';
			cpool.tss_id_map[ite->first][*si] = cpool.id_ve.size()-1; // map<string, map<int, size_t>> tss_id_map
		
        }
	}

}


// ####### count pool #######





// ####### in function "assignCountPool_tss" #######

void count_pool::getregionidmap_tfc(int ups, int downs, int win){

	cut_reg_ve.clear(); // defined in class count_pool | vector<vector<pair<int, int>>> cut_reg_ve;

	for (size_t i = 0; i < id_ve.size(); ++i){
		vector<pair<int, int>> nullvec; // 1 center, 1 element
		cut_reg_ve.push_back(nullvec);
	}
	
	for (map<string, map<int, char>>::iterator ite = tss_strand_map.begin(); ite != tss_strand_map.end(); ++ite){ // tss_strand_map was defined in class count_pool

		string chr = ite->first;
		for (map<int, char>::iterator si = ite->second.begin(); si != ite->second.end(); ++si){
			
			size_t id = tss_id_map[ite->first][si->first]; // map<string, map<int, size_t>> tss_id_map
			int tss = si->first;
			int start = 0;
			int end = 0;
			if (si->second == '+'){
				start = tss - ups;
				end = tss + downs -1;
			}
            else if (si->second == '-'){
				start = tss - downs + 1;
				end = tss + ups;
			}
            else {
				cout << "error strand in count_pool::getreagionidmap_tss " << si->second << endl; exit(1);
			}

			vector<pair<int, int>> cut_reg; 
			vector<int> ct;
			int sp = start;
			int ep = end;

			region_id_map[chr][make_pair(sp, ep)] = id; // map<string, map<pair<int, int>, size_t>> region_id_map
			region_strand_map[chr][make_pair(sp, ep)] = si->second; // map<string, map<pair<int, int>, char>> region_strand_map

			while(start < end){
				int tend = start + win - 1;
				if (tend > end){
				    tend = end;
                }
				cut_reg.push_back(make_pair(start, tend));
				
				ct.push_back(0); // optimize: resize(len) or push_back? avoid segment fault first then efficiency
				start = tend + 1;
			}

			cut_reg_ve[id] = cut_reg;
			counts_table[id] = ct;

		}
	}

}




// overlap between cut region and tag
void assignTagPos(int tagpos, vector<pair<int, int>> &cut_reg, vector<int> &counts, int* colSum, bool* allHitArr, double &allHits, 
					bool* centerHitArr, double &centerHits, int sdStartIdx, int sdEndIdx, int mymaxReadNum){

	
	for (size_t i = 0; i < cut_reg.size(); ++i){
		if (cut_reg[i].second < tagpos){  // region's upper limit is less than tag's position. skip it... next region window
			continue;
		}
		if (cut_reg[i].first > tagpos ){ // region's lowwer limit is beyond tag's position. finish scanning...
			break;
		}
		// if not the situation of being outside of the region
		counts[i] += 1; // there is a tag at at position i
		colSum[i] += 1;

		// ## QC ##
		// for maximum single-cell read number at a specific position
		if (counts[i] > mymaxReadNum){
            fill(counts.begin(), counts.end(), 0);
			break;
        }

		// for coverage filtering
		if (!allHitArr[i]){
			allHitArr[i] = true;
			allHits += 1;

			// for center ± motifShift region
			if (i >= sdStartIdx && i < sdEndIdx){
				centerHitArr[i] = true;
				centerHits += 1;
			}

		}
	
		break;  // finish searching
	}	
	
}



void assignCountPool(map<string, vector<int>> &tagpos, count_pool &cpool, int* colSum, bool* allHitArr, double &allHits,
					bool* centerHitArr, double &centerHits, int sdStartIdx, int sdEndIdx, int maxReadNum){

	int L = 100000; // based on chromosome size
	map<string, map<int, set<pair<int, int>>>> chr_index_region;

	for (map<string, map<pair<int, int>, char>>::iterator ite = cpool.region_strand_map.begin();
		ite != cpool.region_strand_map.end(); ++ite){  // map<string, map<pair<int, int>, char>> region_strand_map. looks like: <chr1, <(1,100), +>>>

		string chr = ite->first;
		for (map<pair<int, int>, char>::iterator si = ite->second.begin(); si != ite->second.end(); ++si){ // for object like <(1,100), +>
			int start = si->first.first;
			int end = si->first.second;
			int index1 = start / L;
			int index2 = end / L;
			chr_index_region[chr][index1].insert(si->first); // like <chr1, <idx1, set{(1,100)...}>>

			if (index2 > index1){
				for (int i = index1+1; i <= index2; ++i){
					chr_index_region[chr][i].insert(si->first);
                }
            }
		}

	} 

	
	for (map<string, vector<int>>::iterator ite = tagpos.begin(); ite != tagpos.end(); ++ite){
		
        string chr = ite->first;
		for (vector<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si){
            
			int index = *si / L;
			if (chr_index_region[chr].find(index) != chr_index_region[chr].end()){
                for (set<pair<int, int>>::iterator ci = chr_index_region[chr][index].begin(); ci != chr_index_region[chr][index].end(); ++ci){
					if (ci->second < *si)
						continue;
					if (ci->first > *si)
						break;
					size_t id = cpool.region_id_map[chr][*ci];
					assignTagPos(*si, cpool.cut_reg_ve[id], cpool.counts_table[id], colSum, allHitArr, allHits, 
								centerHitArr, centerHits, sdStartIdx, sdEndIdx, maxReadNum);
				}
			}

		}
	}
}




// ####### summit distance #######
// vector<double> get_posScore_upDown(int* pileup, int scoreRange, int nuclShift, int posScoreOutterFlank, int posScoreInnerFlank){

//     int colLen = pileup.size();
//     int centerIdx = colLen / 2;
    
//     int psStartIndx = centerIdx - scoreRange - nuclShift;
//     int psEndIndx = centerIdx + scoreRange + nuclShift + 1;


//     // Check the index range before accessing the vector elements
//     if (psStartIndx < 0 || psEndIndx > colLen){
//         // cout << psStartIndx << endl;
//         // cout << psEndIndx << endl;
//         cout << "Error: index out of range" << endl;
//         return -1;
//     }
    
//     // Initialize and resize the vector before using it
//     vector<double> posScore(colLen, 0.0);
//     for (int k = psStartIndx; k < psEndIndx; k++){
//         double tmpMidpointSum_innerBin = 0.0;
//         double tmpMidpointSum_outterBin = 0.0;
//         for (int j = k - posScoreInnerFlank; j <= k + posScoreInnerFlank; j++){
//             tmpMidpointSum_innerBin += pileup[j];
//         }
//         for (int j = k - posScoreOutterFlank; j <= k + posScoreOutterFlank; j++){
//             tmpMidpointSum_outterBin += pileup[j];
//         }
//         if (tmpMidpointSum_outterBin == 0){
//             posScore[k] = 0;
//         }
//         else {
//             posScore[k] = tmpMidpointSum_innerBin / tmpMidpointSum_outterBin;
//         }
//     }
    
//     return posScore;

// }



double get_summitDist(int* pileup, int colLen, int motifShift, int nuclShift, int posScoreOutterFlank, int posScoreInnerFlank, int centerShift){

    int centerIdx = colLen / 2;
    
    int psStartIndx = centerIdx - motifShift - nuclShift;
    int psEndIndx = centerIdx + motifShift + nuclShift + 1;


    // Check the index range before accessing the vector elements
    if (psStartIndx < 0 || psEndIndx > colLen){
        // cout << psStartIndx << endl;
        // cout << psEndIndx << endl;
        cout << "Error: index out of range" << endl;
        return -1;
    }
    
    // Initialize and resize the vector before using it
    vector<double> posScore(colLen, 0.0);
    for (int k = psStartIndx; k < psEndIndx; k++){
        double tmpMidpointSum_innerBin = 0.0;
        double tmpMidpointSum_outterBin = 0.0;
        for (int j = k - posScoreInnerFlank; j <= k + posScoreInnerFlank; j++){
            tmpMidpointSum_innerBin += pileup[j];
        }
        for (int j = k - posScoreOutterFlank; j <= k + posScoreOutterFlank; j++){
            tmpMidpointSum_outterBin += pileup[j];
        }
        if (tmpMidpointSum_outterBin == 0){
            posScore[k] = 0;
        }
        else {
            posScore[k] = tmpMidpointSum_innerBin / tmpMidpointSum_outterBin;
        }
    }
    
    int gwStartIdx = centerIdx - motifShift;
    int gwEndIdx = centerIdx + motifShift + 1;
    
    // Initialize and resize the vector before using it
    vector<double> gwScore(colLen, 0.0);
    for (int k = gwStartIdx; k < gwEndIdx; k++){
        double tmp_k_gwScore = 0.0;
        for (int j = -nuclShift; j <= nuclShift; j++){
            double gw = exp(-pow(j / 20.0, 2) / 2);
            tmp_k_gwScore += posScore[k + j] * gw;
        }
        gwScore[k] = tmp_k_gwScore;
        // cout << gwScore[k] << endl;
    }

    gwScore.assign(gwScore.begin() + gwStartIdx, gwScore.begin() + gwEndIdx);

    int gwCenterIdx = gwScore.size() / 2;
    vector<double> uCwScore(gwScore.begin(), gwScore.begin() + gwCenterIdx);
    vector<double> dCwScore(gwScore.begin() + gwCenterIdx, gwScore.end());
    
    reverse(uCwScore.begin(), uCwScore.end());
    int umidx_flip = max_element(uCwScore.begin(), uCwScore.end()) - uCwScore.begin();
    int umax_center_shift = 0 - umidx_flip - 1;
    int umidx = uCwScore.size() - umidx_flip - 1;
    
    int dmidx = max_element(dCwScore.begin(), dCwScore.end()) - dCwScore.begin();
    int dmax_center_shift = dmidx;
    
    vector<int> slicedPileup(pileup + gwStartIdx, pileup + gwEndIdx); // for array pileup

    int slicedCenterIndex = slicedPileup.size() / 2;
    vector<double> my_array(slicedPileup.begin(), slicedPileup.end());


    // dynamic center
    int index_upstream = distance(my_array.begin(), min_element(my_array.begin() + slicedCenterIndex - centerShift, my_array.begin() + slicedCenterIndex));
    int index_downstream = distance(my_array.begin(), min_element(my_array.begin() + slicedCenterIndex, my_array.begin() + slicedCenterIndex + centerShift));
    if (my_array[index_upstream] < my_array[index_downstream]){
        slicedCenterIndex = index_upstream;
    }
    else {
        slicedCenterIndex = index_downstream;
    }
    

    // 1. from umidx to center
    int uMidPointCountSum = 0;
    double uShiftSum = 0;
    for (int u = umidx; u < slicedCenterIndex; u++){
        int tmpCenterShift = slicedCenterIndex - u;
        int tmpMidPointCount = slicedPileup[u];
        uMidPointCountSum += tmpMidPointCount;
        uShiftSum += tmpMidPointCount * tmpCenterShift;
    }
    double uShiftMean = (uMidPointCountSum == 0) ? -1 * umax_center_shift : uShiftSum / uMidPointCountSum;
    
    // 2. from center to dmidx
    int dMidPointCountSum = 0;
    double dShiftSum = 0;
    // for (int d = slicedCenterIndex + 1; d < dmidx + slicedCenterIndex + 1; d++){
    for (int d = slicedCenterIndex + 1; d < dmidx + slicedCenterIndex; d++){
        int tmpCenterShift = d - slicedCenterIndex;
        int tmpMidPointCount = slicedPileup[d];
        dMidPointCountSum += tmpMidPointCount;
        dShiftSum += tmpMidPointCount * tmpCenterShift;
    }
    double dShiftMean = (dMidPointCountSum == 0) ? dmax_center_shift : dShiftSum / dMidPointCountSum;
    
    double summitDist = uShiftMean + dShiftMean;

    return summitDist;

}



double get_summitDist_vec(vector<int>& pileup, int motifShift, int nuclShift, int posScoreOutterFlank, int posScoreInnerFlank, int centerShift){

    int centerIdx = pileup.size() / 2;
    
    int psStartIndx = centerIdx - motifShift - nuclShift;
    int psEndIndx = centerIdx + motifShift + nuclShift;

    // cout << "pileup.size(): " << pileup.size() << endl; 
    // cout << "psStartIndx: " << psStartIndx << endl;
    // cout << "psEndIndx: " << psEndIndx << endl;

    // Check the index range before accessing the vector elements
    if (psStartIndx < 0 || psEndIndx > pileup.size()) {
        // cout << psStartIndx << endl;
        // cout << psEndIndx << endl;
        cout << "Error: index out of range" << endl;
        return -1;
    }
    
    // Initialize and resize the vector before using it
    vector<double> posScore(pileup.size(), 0.0);
    for (int k = psStartIndx; k < psEndIndx; k++) {
        double tmpMidpointSum_innerBin = 0.0;
        double tmpMidpointSum_outterBin = 0.0;
        for (int j = k - posScoreInnerFlank; j <= k + posScoreInnerFlank; j++) {
            tmpMidpointSum_innerBin += pileup[j];
        }
        for (int j = k - posScoreOutterFlank; j <= k + posScoreOutterFlank; j++) {
            tmpMidpointSum_outterBin += pileup[j];
        }
        if (tmpMidpointSum_outterBin == 0) {
            posScore[k] = 0;
        }
        else {
            posScore[k] = tmpMidpointSum_innerBin / tmpMidpointSum_outterBin;
        }
    }
    
    int gwStartIdx = centerIdx - motifShift;
    int gwEndIdx = centerIdx + motifShift + 1;
    
    // Initialize and resize the vector before using it
    vector<double> gwScore(pileup.size(), 0.0);
    for (int k = gwStartIdx; k < gwEndIdx; k++) {
        double tmp_k_gwScore = 0.0;
        for (int j = -nuclShift; j <= nuclShift; j++) {
            double gw = exp(-pow(j / 20.0, 2) / 2);
            tmp_k_gwScore += posScore[k + j] * gw;
        }
        gwScore[k] = tmp_k_gwScore;
        // cout << gwScore[k] << endl;
    }

    gwScore.assign(gwScore.begin() + gwStartIdx, gwScore.begin() + gwEndIdx);

    int gwCenterIdx = gwScore.size() / 2;
    vector<double> uCwScore(gwScore.begin(), gwScore.begin() + gwCenterIdx);
    vector<double> dCwScore(gwScore.begin() + gwCenterIdx, gwScore.end());
    
    reverse(uCwScore.begin(), uCwScore.end());
    int umidx_flip = max_element(uCwScore.begin(), uCwScore.end()) - uCwScore.begin();
    int umax_center_shift = 0 - umidx_flip - 1;
    int umidx = uCwScore.size() - umidx_flip - 1;
    
    int dmidx = max_element(dCwScore.begin(), dCwScore.end()) - dCwScore.begin();
    int dmax_center_shift = dmidx;
    
    vector<double> slicedPileup(pileup.begin() + gwStartIdx, pileup.begin() + gwEndIdx);

    int slicedCenterIndex = slicedPileup.size() / 2;
    vector<double> my_array(slicedPileup.begin(), slicedPileup.end());


    // dynamic center
    int index_upstream = distance(my_array.begin(), min_element(my_array.begin() + slicedCenterIndex - centerShift, my_array.begin() + slicedCenterIndex));
    int index_downstream = distance(my_array.begin(), min_element(my_array.begin() + slicedCenterIndex, my_array.begin() + slicedCenterIndex + centerShift));
    if (my_array[index_upstream] < my_array[index_downstream]){
        slicedCenterIndex = index_upstream;
    }
    else {
        slicedCenterIndex = index_downstream;
    }
    

    // 1. from umidx to center
    int uMidPointCountSum = 0;
    double uShiftSum = 0;
    for (int u = umidx; u < slicedCenterIndex; u++) {
        int tmpCenterShift = slicedCenterIndex - u;
        int tmpMidPointCount = slicedPileup[u];
        uMidPointCountSum += tmpMidPointCount;
        uShiftSum += tmpMidPointCount * tmpCenterShift;
    }
    double uShiftMean = (uMidPointCountSum == 0) ? -1 * umax_center_shift : uShiftSum / uMidPointCountSum;
    
    // 2. from center to dmidx
    int dMidPointCountSum = 0;
    double dShiftSum = 0;
    for (int d = slicedCenterIndex + 1; d < dmidx + slicedCenterIndex + 1; d++) {
        int tmpCenterShift = d - slicedCenterIndex;
        int tmpMidPointCount = slicedPileup[d];
        dMidPointCountSum += tmpMidPointCount;
        dShiftSum += tmpMidPointCount * tmpCenterShift;
    }
    double dShiftMean = (dMidPointCountSum == 0) ? dmax_center_shift : dShiftSum / dMidPointCountSum;
    
    double summitDist = uShiftMean + dShiftMean;

    return summitDist;
}



double get_summitDist_new(int* pileup, int colLen, int motifShift, int nuclShift, int posScoreOutterFlank, int posScoreInnerFlank, int centerShift){

    int centerIdx = colLen / 2;
    int psStartIndx = centerIdx - motifShift - nuclShift;
    int psEndIndx = centerIdx + motifShift + nuclShift + 1;

    if (psStartIndx < 0 || psEndIndx > colLen) {
        std::cout << "Error: index out of range" << std::endl;
        return -1;
    }

    vector<double> posScore(colLen, 0.0);
    for (int k = psStartIndx; k < psEndIndx; ++k) {
        double tmpMidpointSum_innerBin = std::accumulate(pileup + k - posScoreInnerFlank, pileup + k + posScoreInnerFlank + 1, 0.0);
        double tmpMidpointSum_outterBin = std::accumulate(pileup + k - posScoreOutterFlank, pileup + k + posScoreOutterFlank + 1, 0.0);

        posScore[k] = (tmpMidpointSum_outterBin == 0) ? 0 : tmpMidpointSum_innerBin / tmpMidpointSum_outterBin;
    }

    int gwStartIdx = centerIdx - motifShift;
    int gwEndIdx = centerIdx + motifShift + 1;

    vector<double> gwScore(colLen, 0.0);
    for (int k = gwStartIdx; k < gwEndIdx; ++k) {
        double tmp_k_gwScore = 0.0;
        for (int j = -nuclShift; j <= nuclShift; ++j) {
            double gw = exp(-pow(j / 20.0, 2) / 2);
            tmp_k_gwScore += posScore[k + j] * gw;
        }
        gwScore[k] = tmp_k_gwScore;
    }

    gwScore.assign(gwScore.begin() + gwStartIdx, gwScore.begin() + gwEndIdx);

    int gwCenterIdx = gwScore.size() / 2;
    vector<double> uCwScore(gwScore.begin(), gwScore.begin() + gwCenterIdx);
    vector<double> dCwScore(gwScore.begin() + gwCenterIdx, gwScore.end());

    std::reverse(uCwScore.begin(), uCwScore.end());
    int umidx_flip = std::max_element(uCwScore.begin(), uCwScore.end()) - uCwScore.begin();
    int umax_center_shift = 0 - umidx_flip - 1;
    int umidx = uCwScore.size() - umidx_flip - 1;

    int dmidx = std::max_element(dCwScore.begin(), dCwScore.end()) - dCwScore.begin();
    int dmax_center_shift = dmidx;

    vector<int> slicedPileup(pileup + gwStartIdx, pileup + gwEndIdx);
    int slicedCenterIndex = slicedPileup.size() / 2;

    vector<double> my_array(slicedPileup.begin(), slicedPileup.end());

    int index_upstream = std::distance(my_array.begin(), std::min_element(my_array.begin() + slicedCenterIndex - centerShift, my_array.begin() + slicedCenterIndex));
    int index_downstream = std::distance(my_array.begin(), std::min_element(my_array.begin() + slicedCenterIndex, my_array.begin() + slicedCenterIndex + centerShift));
    slicedCenterIndex = (my_array[index_upstream] < my_array[index_downstream]) ? index_upstream : index_downstream;

    int uMidPointCountSum = 0;
    double uShiftSum = 0;
    for (int u = umidx; u < slicedCenterIndex; ++u) {
        int tmpCenterShift = slicedCenterIndex - u;
        int tmpMidPointCount = slicedPileup[u];
        uMidPointCountSum += tmpMidPointCount;
        uShiftSum += tmpMidPointCount * tmpCenterShift;
    }
    double uShiftMean = (uMidPointCountSum == 0) ? -1 * umax_center_shift : uShiftSum / uMidPointCountSum;

    int dMidPointCountSum = 0;
    double dShiftSum = 0;
    for (int d = slicedCenterIndex + 1; d < dmidx + slicedCenterIndex + 1; ++d) {
        int tmpCenterShift = d - slicedCenterIndex;
        int tmpMidPointCount = slicedPileup[d];
        dMidPointCountSum += tmpMidPointCount;
        dShiftSum += tmpMidPointCount * tmpCenterShift;
    }
    double dShiftMean = (dMidPointCountSum == 0) ? dmax_center_shift : dShiftSum / dMidPointCountSum;

    return uShiftMean + dShiftMean;
}

// ####### summit distance #######





// ####### for nucl. occupancy score #######

void assignTagPos_occuScore(int tagpos, vector<pair<int, int>> &cut_reg, vector<int> &counts){

	
	for (size_t i = 0; i < cut_reg.size(); ++i){
		if (cut_reg[i].second < tagpos){  // region's upper limit is less than tag's position. skip it... next region window
			continue;
		}
		if (cut_reg[i].first > tagpos ){ // region's lowwer limit is beyond tag's position. finish scanning...
			break;
		}
		// if not the situation of being outside of the region
		counts[i] += 1; // there is a tag at at position i
	
		break;
	}	
	
}



void assignCountPool_occuScore(map<string, vector<int>> &tagpos, count_pool &cpool){

	int L = 100000; // based on chromosome size
	map<string, map<int, set<pair<int, int>>>> chr_index_region;

	for (map<string, map<pair<int, int>, char>>::iterator ite = cpool.region_strand_map.begin();
		ite != cpool.region_strand_map.end(); ++ite){  // map<string, map<pair<int, int>, char>> region_strand_map. looks like: <chr1, <(1,100), +>>>

		string chr = ite->first;
		for (map<pair<int, int>, char>::iterator si = ite->second.begin(); si != ite->second.end(); ++si){ // for object like <(1,100), +>
			int start = si->first.first;
			int end = si->first.second;
			int index1 = start / L;
			int index2 = end / L;
			chr_index_region[chr][index1].insert(si->first); // like <chr1, <idx1, set{(1,100)...}>>

			if (index2 > index1){
				for (int i = index1+1; i <= index2; ++i){
					chr_index_region[chr][i].insert(si->first);
                }
            }
		}

	} 

	
	for (map<string, vector<int>>::iterator ite = tagpos.begin(); ite != tagpos.end(); ++ite){
		
        string chr = ite->first;
		for (vector<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si){
            
			int index = *si / L;
			if (chr_index_region[chr].find(index) != chr_index_region[chr].end()){
                for (set<pair<int, int>>::iterator ci = chr_index_region[chr][index].begin(); ci != chr_index_region[chr][index].end(); ++ci){
					if (ci->second < *si)
						continue;
					if (ci->first > *si)
						break;
					size_t id = cpool.region_id_map[chr][*ci];
					assignTagPos_occuScore(*si, cpool.cut_reg_ve[id], cpool.counts_table[id]);
				}
			}

		}
	}
}



double calculateCpoolMatSum(const count_pool& mypool){
    double sum = 0.0;
    for (const auto &row : mypool.counts_table) {
        for (const auto &element : row) {
            sum += element;
        }
    }
    return sum;
}



double nuclFraction(const count_pool& centerCpool, const count_pool& flankCpool){
    double centerSum = calculateCpoolMatSum(centerCpool);
    double flankSum = calculateCpoolMatSum(flankCpool);
    return centerSum / flankSum;
}



void writeRawCountTable(const string& prefix, const vector<vector<int>>& rawCountTable){

    string countTableFile = prefix + ".rawCountTable.txt";
    std::ofstream outf(countTableFile.data());

    if (!rawCountTable.empty()){
        for (int i = 0; i < rawCountTable.size(); i++){
            for (int j = 0; j < rawCountTable[i].size(); j++){
                outf << rawCountTable[i][j] << "\t";
            }
            outf << "\n";
        }
    } 
    else {
        std::cout << "prefix: " << prefix << " | nucleosome fraction matrix is empty, nothing to output." << std::endl;
    }
    outf.close();

}



void writeVec(const string& prefix, const vector<int>& myvec){

    string countTableFile = prefix + ".rawCountTable.txt";
    std::ofstream outf(countTableFile.data());

    if (!myvec.empty()){
        for (int i = 0; i < myvec.size(); i++){
            outf << myvec[i] << "\n";
        }
    } 
    else {
        std::cout << "prefix: " << prefix << " | nucleosome fraction matrix is empty, nothing to output." << std::endl;
    }
    outf.close();

}



// process each cell | for multi-threading
void processCell(const string& cellBed, const vector<count_pool>& cpool_center_vec, 
                 const string& prefix, const vector<string>& centerName_vec) {

    std::filesystem::path p(cellBed);
    std::string tmp_cell_basename = p.filename().string();
    string tmp_cell_dir = prefix + tmp_cell_basename + "/";
    std::filesystem::create_directories(tmp_cell_dir);
    
    thread_local map<string, vector<int>> tagposmap;
    readAlignTagFile(cellBed, tagposmap);
    
    for (int k = 0; k < cpool_center_vec.size(); k++) {
        count_pool tmp_cpool_center = cpool_center_vec[k];
        
        assignCountPool_occuScore(tagposmap, tmp_cpool_center);
        
        vector<vector<int>> tmp_center_rawcount = tmp_cpool_center.counts_table;
        string tmp_table_prefix = tmp_cell_basename + "_" + centerName_vec[k];
        string tmp_save_pre = tmp_cell_dir + tmp_table_prefix;
        writeRawCountTable(tmp_save_pre, tmp_center_rawcount);
    }

}

// ####### for nucl. occupancy score #######
