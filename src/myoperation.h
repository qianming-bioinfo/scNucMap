#ifndef MYOPERATION_H
#define MYOPERATION_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <mutex>
#include <chrono>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <filesystem>
#include <future>
#include <queue>

using namespace std;



// ####### read file #######
vector<string> getFiles(const string& folder, const string& suffix, string fileListName = "my");
bool isFileEmpty(ifstream& file);
vector<string> getPathFromMetafile(const string& filename);
vector<vector<string>> splitVector(const vector<string>& vec, int numParts);

vector<string> parseString(string& instr);
void readAlignTagFile(string infile, map<string, vector<int>> &tagpos);
void readRegionCenter(string &infile, map<string, vector<int>> &regioncenter);




// ####### count pool #######
class count_pool{

public:
	
	vector<string> id_ve;
	vector<vector<int>> counts_table;
	map<string, map<int, char>> tss_strand_map;
	map<string, map<int, size_t>> tss_id_map;
	
	map<string, map<pair<int, int>, size_t>> region_id_map;
	map<string, map<pair<int, int>, char>> region_strand_map;
	vector<vector<pair<int, int>>> cut_reg_ve;

	void getregionidmap_tfc(int ups, int downs, int win);
	
};




// ####### multi-threading #######
class ThreadPool {
public:
    ThreadPool(size_t numThreads) : stop(false) {
        for (size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back(
                [this] {
                    for (;;) {
                        std::function<void()> task;

                        {
                            std::unique_lock<std::mutex> lock(this->queue_mutex);
                            this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                            if (this->stop && this->tasks.empty())
                                return;
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                        }

                        task();
                    }
                }
            );
        }
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<void> {
        using return_type = void; // 这里指定返回类型为 void

        auto task = std::make_shared<std::packaged_task<return_type()>>(std::bind(std::forward<F>(f), std::forward<Args>(args)...));

        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            // 不允许在池停止后加入任务
            if (stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            tasks.emplace([task] { (*task)(); });
        }
        condition.notify_one();
        return task->get_future(); // 返回对应的 future
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread& worker : workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};




// ####### processing #######
string int2str(int i);
void getTFCcpool(map<string, vector<int>> &tfcenter, count_pool &cpool);
void assignTagPos(int tagpos, vector<pair<int, int>> &cut_reg, vector<int> &counts, int* colSum, bool* allHitArr, double &allHits, 
					bool* centerHitArr, double &centerHits, int sdStartIdx, int sdEndIdx, int mymaxReadNum = 3);
void assignCountPool(map<string, vector<int>> &tagpos, count_pool &cpool, int* colSum, bool* allHitArr, double &allHits,
					bool* centerHitArr, double &centerHits, int sdStartIdx, int sdEndIdx, int maxReadNum = 3);


// ####### summit distance #######
// vector<double> get_posScore_upDown(int* pileup, int scoreRange, int nuclShift = 73, 
//                                 int posScoreOutterFlank = 100, int posScoreInnerFlank = 15);
double get_summitDist(int* pileup, int colLen, int motifShift = 200, int nuclShift = 73, 
        			int posScoreOutterFlank = 100, int posScoreInnerFlank = 15, int centerShift = 20);
double get_summitDist_vec(vector<int>& pileup, int motifShift = 200,  int nuclShift = 73, 
                        int posScoreOutterFlank = 100, int posScoreInnerFlank = 15, int centerShift = 20);
double get_summitDist_new(int* pileup, int colLen, int motifShift = 200, int nuclShift = 73, 
                      int posScoreOutterFlank = 100, int posScoreInnerFlank = 15, int centerShift = 20);


// ####### for nucl. occupancy score #######
void assignTagPos_occuScore(int tagpos, vector<pair<int, int>> &cut_reg, vector<int> &counts);
void assignCountPool_occuScore(map<string, vector<int>> &tagpos, count_pool &cpool);
double calculateCpoolMatSum(const count_pool& mypool);
double nuclFraction(const count_pool& centerCpool, const count_pool& flankCpool);
void writeRawCountTable(const string& prefix, const vector<vector<int>>& rawCountTable);
void writeVec(const string& prefix, const vector<int>& myvec);
void processCell(const string& cellBed, const vector<count_pool>& cpool_center_vec, const string& prefix, const vector<string>& centerName_vec);

#endif