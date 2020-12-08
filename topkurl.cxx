/* Assumptions:
 * - There are no excessively long urls.
 *   As the map uses only size to determine whether memory exceeds...
 * - The file should not be more than a hundred times larger than the memory.
 *   As the max intermediate files num is hard coded to 1024
 */
// Author: zhuo1ang

#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

//! \brief Determine the number of files based on the given size
//! \param file_size input file size in bytes
//! \param mem_size memory size in bytes
//! \return div by 10 to make sure that no single file exceeds the memory
#define FileNum() ((file_size) / ((mem_size) / 10))

using namespace std;

const int max_file = 1024;
const int max_thread = 128;

char *file_path;
bool clean_up = true;

size_t file_size = (size_t)100 * 1e9;
size_t mem_size = (size_t)1 * 1e9;
size_t k = 100;

unsigned int thread_num = 8;  // in my machine

typedef std::pair<std::string, int> url_cnt_pair;

struct url_pair_comp {
    //! \brief the overload of pair, used by priority_queue to find the max 100 target in each small file
    bool operator()(url_cnt_pair const &a, url_cnt_pair const &b) const {
        if (a.second < b.second)
            return true;
        else if (a.second > b.second)
            return false;
        else
            return a.first.size() >= b.first.size();
    }
};

auto Hash = std::hash<std::string>{};
ofstream hash_sorted_files[max_file];
std::map<string, size_t> mp;
std::priority_queue<url_cnt_pair, std::vector<url_cnt_pair>, url_pair_comp> pq;

std::thread th[max_thread];
std::mutex hash_files_mutex[max_file];

inline void OpenHashFileStreams() {
    for (int i = 0; i < FileNum(); i++) {
        hash_sorted_files[i].open(to_string(i));
    }
}
inline void CloseHashFileStreams() {
    for (int i = 0; i < FileNum(); i++) {
        hash_sorted_files[i].close();  // = fstream(to_string(i), ios::in | ios::out);
    }
}

//! \brief Generate small files which content strings are sorted by hash code. It
//! \param input_file_path the origin input file (~100G)
//! \param from index of the starting byte
//! \param to one plus, the index of the last byte
void SortByHashValue(string input_file_path, int from, int to, int idx) {
    ifstream urls;
    urls.open(input_file_path);
    urls.seekg(from);

    if (from != 0) {
        string _;
        getline(urls, _);
    }

    from = urls.tellg();

    while (!urls.eof() && from <= to) {
        string tmp;
        getline(urls, tmp);
        hash_files_mutex[Hash(tmp) % FileNum()].lock();
        hash_sorted_files[Hash(tmp) % FileNum()] << tmp << endl;
        hash_files_mutex[Hash(tmp) % FileNum()].unlock();
        from = urls.tellg();
    }
}

//! \brief Get down to each file and calc the top-K url of each file.
//         As it requires a hashmap to hold the data, these procedures could not run in parallel.
void CalcPriorityQueue() {
    for (int i = 0; i < FileNum(); i++) {
        std::ifstream file;
        file.open(to_string(i));
        while (!file.eof()) {
            std::string tmp;
            file >> tmp;
            if (tmp.size() == 0)
                continue;
            if (!mp.count(tmp))
                mp[tmp] = 1;
            else
                mp[tmp]++;
        }
        for (auto [key, value] : mp) {
            pq.push(make_pair(key, value));
        }

        vector<url_cnt_pair> tmpv;
        for (int j = 0; j < k && !pq.empty(); j++) {
            tmpv.push_back(pq.top());
            pq.pop();
        }
        pq = priority_queue<url_cnt_pair, vector<url_cnt_pair>, url_pair_comp>();  // reset it
        for (auto p : tmpv)
            pq.push(p);

        mp.clear();
    }
}

//! \param filename file path in fs
//! \return file size in bytes
std::ifstream::pos_type GetFileSize(const char *filename) {
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

//! \brief Print help messages
void PrintHelpMsg() {
    cerr << "Usage: topurl [options]... [FILE]" << endl;
    cerr << "Count the most frequent URLs" << endl;
    cerr << endl;
    cerr << "  -k <number>\t\tcount top <number> URLs, default to 100" << endl;
    cerr << "  -f <file_size>\tthe approximate size of the file, in bytes, default to 100G" << endl;
    cerr << "  -m <mem_size>\t\tthe approximate size of the memory, in bytes, default to 1G" << endl;
    cerr << "  -j <thread_num>\tallow thread_num jobs at once" << endl;
    cerr << "  -c\t\t\tclean up temporary files" << endl;
    cerr << "  -h\t\t\tprint this help message" << endl;
    exit(1);
}

int main(int argc, char *argv[]) {
    char c;
    while ((c = getopt(argc, argv, "chk:m:j:")) != -1) {
        switch (c) {
            case 'h':
                PrintHelpMsg();
                break;
            case 'c':
                clean_up = false;
                break;
            case 'k':
                k = atoi(optarg);
                break;
            case 'm':
                mem_size = atoi(optarg);
                break;
            case 'j':
                thread_num = atoi(optarg);
                break;
            case '?':
                if (optopt == 'k')
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint(optopt))
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);

                cerr << endl;
                PrintHelpMsg();
                exit(1);
            default:
                abort();
        }
    }

    if (optind < argc - 1) {
        cerr << "Warning: Ignoring extra files." << endl;
    } else if (optind == argc) {  // there's no chance that optind > argc
        cerr << "Error: No input file, aborted." << endl;
        cerr << endl;
        PrintHelpMsg();
        exit(1);
    }
    file_path = argv[optind];
    file_size = GetFileSize(file_path);

    OpenHashFileStreams();

    for (int i = 0; i < thread_num; i++) {
        // SortByHashValue(file_path);
        th[i] = std::thread(SortByHashValue,
                            file_path,
                            file_size / thread_num * i,
                            i + 1 != thread_num ? file_size / thread_num * (i + 1) : file_size,
                            i);
    }
    for (int i = 0; i < thread_num; i++) {
        th[i].join();
    }
    CloseHashFileStreams();
    CalcPriorityQueue();

    while (!pq.empty()) {
        std::cout << pq.top().first << " " << pq.top().second << endl;
        pq.pop();
    }

    if (clean_up) {
        std::string tmp_files;
        for (int i = 0; i < max_file; i++) {
            tmp_files += to_string(i) + " ";
        }
        system((string("rm ") + tmp_files + " 2> /dev/null").c_str());
    }

    return 0;
}