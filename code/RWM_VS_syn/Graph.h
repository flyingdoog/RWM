#pragma once
#ifndef _HEAD_GRAPH
#define _HEAD_GRAPH
#define tr1ns std
#include<unordered_map>
#include<string>
#include"tools.h"
#include<assert.h>
using namespace std;
struct sparsevec {
	//mxIndex ->int
	typedef std::unordered_map<unsigned int, double> map_type;
	map_type map;
	/** Get an element and provide a default value when it doesn't exist
	* This command does not insert the element into the vector
	*/
	double get(unsigned int index, double default_value = 0.0) {
		map_type::iterator it = map.find(index);
		if (it == map.end()) {
			return default_value;
		}
		else {
			return it->second;
		}
	}

	/** Compute the sum of all the elements
	* Implements compensated summation
	*/
	double sum() {
		double s = 0.;
		for (map_type::iterator it = map.begin(), itend = map.end(); it != itend; ++it) {
			s += it->second;
		}
		return s;
	}

	/** Compute the max of the element values
	* This operation returns the first element if the vector is empty.
	*/
	unsigned int max_index() {
		unsigned int index = 0;
		double maxval = std::numeric_limits<double>::min();
		for (map_type::iterator it = map.begin(), itend = map.end(); it != itend; ++it) {
			if (it->second>maxval) { maxval = it->second; index = it->first; }
		}
		return index;
	}
};

struct mtuple {
	unsigned int row, col;
	double w;
	bool operator < (const mtuple& str) const
	{
		if (row < str.row)
			return true;
		else if (str.row < row)
			return false;
		else if (col < str.col)
			return true;
		else
			return false;
	}
};
typedef mtuple myEdge;

struct sparserow {
	unsigned int nodeSize, edgeSize;
	bool weighted;
	bool directed;
	unsigned int *indexs;//indexs
	unsigned int *nbs;//nbs
	double *ws;
	double *volumns;
public:
	void check();
	sparserow() {
		nodeSize = edgeSize = 0;
		weighted = false;
		directed = false;
		indexs = NULL;
		nbs = NULL;
	}

	uint degree(uint id) const{
		//assert(id < nodeSize+1);
		if (id > nodeSize + 1)
			return 0;
		return indexs[id + 1] - indexs[id];
	}

	double vol(uint id) const {
		//assert(id < nodeSize+1);
		if (id > nodeSize + 1)
			return 0;
		return volumns[id];
	}


};

void createGraph(const std::string, sparserow &);
void createGraph(ifstream &fin, sparserow & graph, char endchar);
void createGraph(ifstream &fin, sparserow & graph, char endchar, bool);
void createGraph(const sparserow & from, sparserow & to, bool reverse);
void createGraph(ifstream &fin, sparserow & graph, char endchar, double, set<myEdge> &);
double extractPart(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double theta, uint seed, bool *visited);
void calConductance(const sparserow &graph,vector<double> &conductances,const vector<int_double_pair_sorted_by_score> tep);
double calConductance(const sparserow &graph, const set<uint> &cmty);
double calConductance(const sparserow &graph, const vector<uint> &cmty);
void calConductance(const sparserow &graph, vector<double> &conductances, const vector<int_double_pair_sorted_by_score> tep, const sparserow& inter);
struct local_pagerank_stats {
	double conductance;
	double volume;
	double support;
	double steps;
	double eps;
	double cut;
};

/*
res = alpha*P*x
*/
void multi(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double alpha);
void multi(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double alpha, bool);
void amulti(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double alpha, double theta, uint seed, bool *visited, bool self);
void addArray(map<uint, double>&, const map<uint, double>&, double);

#endif // !_HEAD_GRAPH

