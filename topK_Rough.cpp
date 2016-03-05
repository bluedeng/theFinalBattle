//2016-02-24, dtz
//just for take the answer of top k
//basical method, no optimization, rough one

#include "Time.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <queue>

#include <stdlib.h>
#include <string.h>

using namespace std;

struct stringInfo{
	string s_string;
	unsigned s_length;

	stringInfo(string s, unsigned length) {
		s_string = s;
		s_length = length;
	};
};

struct queue_entry{
	unsigned m_sid;
	unsigned m_dist;

	queue_entry(unsigned sid, unsigned dist) {
		m_sid = sid;
		m_dist = dist;
	};

	bool operator < (const queue_entry &other) const {
		return m_dist < other.m_dist;
	};
};

vector<stringInfo> dataset;
vector<stringInfo> queryset;
vector<int> tk;
int topk;
unsigned gl_maxLen = 0;
bool output = false;

unsigned processed;
int threshold;
int** edcost;
priority_queue<queue_entry> m_queue;

map<unsigned, unsigned> candidate;
unsigned candidate_average;
map<unsigned, unsigned> max_ed;
double max_ed_average;
vector<double> time_all;
double time_all_average;

//usage instruction
void usage() {
	cout << "************************************************" << endl;
	cout << "------------------------------------------------" << endl;
	cout << "Usage: simple_index [OPTION] ${INPUTDB} ${QUERY}" << endl;
//	cout << "-k val     set the value for top-k heap         " << endl;
	cout << "-o         output the result set of sequence ids" << endl;
//	cout << "      we will process the topK for 5 only       " << endl;
	cout << " we will process the topK for 5, 10, 20, 30, 40 " << endl;
	cout << "------------------------------------------------" << endl;
	cout << "************************************************" << endl;
}

//get the input parameters: output
void parseOptions(int argc, char* argv[]) {
	if (argc > 3)
		for (int i = 1; i < argc - 2; i++) {
			/*
			if (strncmp(argv[i], "-k", 2) == 0) {
				tk.push_back(atoi(argv[i + 1]));
				i++;
			}
			*/
			if (strncmp(argv[i], "-o", 2) == 0)
				output = true;
		}
}

//read the input dataset file
void read(string filename, vector<stringInfo>& data) {
	ifstream _fs;
	_fs.open(filename.c_str(), ios::in);
	if (_fs.fail()) {
		cerr << "Error: Failed to open the data file: " << filename << endl;
		cerr << "Please check the file name and restart this program" << endl;
		cerr << "The program will be terminated!" << endl;
		exit(-1);
	}

	string line;
	while (!_fs.eof()) {
		getline(_fs, line);

		if (line.empty())
			continue;

		for (unsigned i = 0; i < line.length(); i++)
			if (line[i] >= 'A' && line[i] <= 'Z')
				line[i] += 32;
		
		unsigned l = line.length();

		gl_maxLen = gl_maxLen > l ? gl_maxLen : l;

		data.push_back(stringInfo(line, l));
	}

	_fs.close();
}

//change int to string
string changeI2S(int i) {
	stringstream ss;
	string temp;
	ss << i;
	ss >> temp;
	ss.clear();
	return temp;
}

//get the edit distance between s1 and s2, if the intermediate result is large than t, return it.
int getRealEditDistance(const stringInfo& sI1, const stringInfo& sI2, const int& t) {
	const size_t len1 = sI1.s_length, len2 = sI2.s_length;
	const string s1 = sI1.s_string, s2 = sI2.s_string;
	int mincost;

	edcost[0][0] = 0;
	for (unsigned int i = 1; i <= len1; i++) edcost[i][0] = i;
	for (unsigned int i = 1; i <= len2; i++) edcost[0][i] = i;

	for (unsigned int i = 1; i <= len1; i++) {
		mincost = len2;
		for (unsigned int j = 1; j <= len2; j++) {
			edcost[i][j] = min(min(edcost[i - 1][j] + 1, edcost[i][j - 1] + 1),
				edcost[i - 1][j - 1] + (s1[i - 1] == s2[j - 1] ? 0 : 1));
			mincost = (mincost <= edcost[i][j] ? mincost : edcost[i][j]);
		}

		if (mincost > t)
			return mincost;
	}

	return edcost[len1][len2];
}

//initial the threshold for query 
//the length of query temporary as the t-value in ed function
void init_threshold(const stringInfo& query) {
	processed = 0;
	int ed;

	for (int i = 0; i < topk; i++) {
		ed = getRealEditDistance(dataset[i], query, query.s_length);
		processed++;
		m_queue.push(queue_entry((unsigned)i, (unsigned)ed));
	}

	threshold = m_queue.top().m_dist;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		usage();
		exit(-1);
	}

	parseOptions(argc, argv);

	//read the data file
	read(argv[argc - 2], dataset);
	if (dataset.size() == 0) {
		cerr << "Error: Failed to open the data file: " << argv[argc - 2] << endl;
		cerr << "Please check the file name and restart this program" << endl;
		cerr << "The program will be terminated!" << endl;
		exit(-1);
	}

	//read the query file
	read(argv[argc - 1], queryset);
	if (queryset.size() == 0) {
		cerr << "Error: Failed to open the data file: " << argv[argc - 1] << endl;
		cerr << "Please check the file name and restart this program" << endl;
		cerr << "The program will be terminated!" << endl;
		exit(-1);
	}

	//build the edcost matrix, prepare for the calculation
	edcost = new int *[gl_maxLen + 1];
	for (unsigned i = 0; i < gl_maxLen + 1; i++)
		edcost[i] = new int[gl_maxLen + 1];

	string stemp = argv[argc - 1];
	ofstream fout((stemp + "_Rough").c_str());
	fout << "K,  Average_max_ed, Average_candidate_number,  Average_time" << endl;

	//push values into topk array
//	tk.push_back(5);
	tk.push_back(5); tk.push_back(10); tk.push_back(20); tk.push_back(30); tk.push_back(40);
	for (unsigned j = 0; j < tk.size(); j++) {
		topk = tk[j];
		int ed = 0;

		ofstream fdetail((stemp + "_Rough_" + changeI2S(topk)).c_str());

		candidate.clear();
		max_ed.clear();
		time_all.clear();

		for (unsigned i = 0; i < queryset.size(); i++) {
			cout << "processing the " << i << "th" << endl;

			class TSINGHUA_CLIPSE_UTIL::TimeRecorder time;
			
			init_threshold(queryset[i]);

			for (unsigned u = topk; u < dataset.size(); u++) {
				if (queryset[i].s_length > dataset[u].s_length + threshold)
					continue;
				if (dataset[u].s_length > queryset[i].s_length + threshold)
					break;

				ed = getRealEditDistance(dataset[u], queryset[i], threshold);
				processed++;
				if (ed >= threshold)
					continue;

				m_queue.pop();
				m_queue.push(queue_entry((unsigned)u, (unsigned)ed));
				threshold = m_queue.top().m_dist;
			}
			time.check();

			candidate[processed] += 1;
			time_all.push_back(time.diffTime(0, 1));
			max_ed[m_queue.top().m_dist] += 1;

			if (output) {
				cout << "Query for: " << queryset[i].s_string << endl;
				fdetail << "Query for: " << queryset[i].s_string << endl;
				while (!m_queue.empty()) {
					const queue_entry& entry = m_queue.top();
					cout << entry.m_sid << " " << dataset[entry.m_sid].s_string << " " << entry.m_dist << endl;
					fdetail << dataset[entry.m_sid].s_string << " " << entry.m_dist << endl;
					m_queue.pop();
				}
			}
		}
		fdetail.close();

		//statics
		map<unsigned, unsigned>::iterator iter;
		unsigned temp1 = 0;
		double temp2 = 0.0;
		unsigned len = queryset.size();

		for (iter = candidate.begin(); iter != candidate.end(); iter++)
			temp1 += iter->first * iter->second;
		candidate_average = temp1 / len;

		temp1 = 0;
		for (iter = max_ed.begin(); iter != max_ed.end(); iter++)
			temp1 += iter->first * iter->second;
		max_ed_average = (double)temp1 / len;

		for (unsigned i = 0; i < time_all.size(); i++)
			temp2 += time_all[i];
		time_all_average = temp2 / len;

		//output
		fout << topk << "  " << max_ed_average << "  " << candidate_average << "  " << time_all_average << endl;
		cout << "Over for topk:    " << topk << endl;
	}
	fout.close();
}