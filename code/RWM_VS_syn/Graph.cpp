#include "stdafx.h"
#include"Graph.h"
#include<fstream>
#include<string>
#include"tools.h"
#include<algorithm>
using namespace std;
bool debug_mtuple = false;
int mtplcnt = 0;
extern uint topk;
extern int min_cmty_size, max_cmty_size;

void createGraph(const string path, sparserow & graph) {
	uint hashroot = 39569;
	set<uint> sedges;
	sedges.insert(0);
	ifstream fin(path);
	string line;
	vector<mtuple> edges;
	int maxid = -1;
	while (!fin.eof()) {
		getline(fin, line);
		if (line.length() < 2 || line.c_str()[0] == '#')
			continue;
		vector<string> ss;
		split(ss, line, '\t');
		int rid = atoi(ss[0].c_str());
		int cid = atoi(ss[1].c_str());
		if (cid > maxid) {
			maxid = cid;
		}

		if (rid > maxid) {
			maxid = rid;
		}

		double w = 1.0;
		if (ss.size() >= 3) {
			w = atof(ss[2].c_str());
		}
		mtuple ele;
		ele.row = rid;
		ele.col = cid;
		ele.w = w;
		uint hashc = rid*hashroot + cid;
		if (sedges.find(hashc) == sedges.end()) {
			edges.push_back(ele);
			sedges.insert(hashc);
		}
		if (!graph.directed) {
			ele.row = cid;
			ele.col = rid;

			hashc = cid*hashroot + rid;
			if (sedges.find(hashc) == sedges.end()) {
				edges.push_back(ele);
				sedges.insert(hashc);
			}
		}

	}

	sort(edges.begin(), edges.end());

	fin.close();

	graph.nodeSize = maxid;
	graph.edgeSize = edges.size();
	graph.nbs = (unsigned int*)malloc(sizeof(unsigned int)*graph.edgeSize);
	graph.indexs = (unsigned int*)malloc(sizeof(unsigned int)*(graph.nodeSize + 2));
	//graph.ws = (double*)malloc(sizeof(double)*graph.edgeSize);

	cout << graph.nodeSize << endl;


	graph.indexs[0] = 0;
	int old = 0;
	for (size_t i = 0; i < graph.edgeSize; i++) {
		graph.nbs[i] = edges[i].col;
		//graph.ws[i] = edges[i].w;
		if (edges[i].row != old) {
			for (unsigned int next = old + 1; next <= edges[i].row; next++) {
				graph.indexs[next] = i;
			}
		}
		old = edges[i].row;
	}
	for (unsigned int next = old + 1; next <= maxid + 1; next++) {
		graph.indexs[next] = graph.edgeSize;
	}
}
void createGraph(ifstream &fin, sparserow & graph, char endchar) {
	createGraph(fin, graph, endchar, false);
}

/*
* Create graph from file stream.
* skip the lines startswith '#'
* end with the lines startswith '-'
*/
void createGraph(ifstream &fin, sparserow & graph, char endchar, bool self) {
	uint hashroot = 39569;
	set<uint> sedges;
	string line = "";
	vector<mtuple> edges;
	int maxid = -1;
	string title;
	graph.weighted = false;
	while (!fin.eof()) {
		getline(fin, line);
		if (line.c_str()[0] == endchar) {
			if (DEBUG_CREATE_GRAPH)
				cout << "begin to read " << line << endl;
			title = line;
			break;
		}
		if (line.length() < 2 || line.c_str()[0] == '#')
			continue;
		vector<string> ss;
		split(ss, line, '\t');
		if (ss.size() <2) {
			cout << "skip :" << line << endl;
			continue;
		}
		int rid = atoi(ss[0].c_str());
		int cid = atoi(ss[1].c_str());
		if (rid > maxid) {
			maxid = rid;
		}

		double w = 1.0;
		if (ss.size() >= 3) {
			w = atof(ss[2].c_str());
			graph.weighted = true;
		}
		mtuple ele;
		ele.row = rid;
		ele.col = cid;
		ele.w = w;
		uint hashc = rid*hashroot + cid;
		if (sedges.find(hashc) == sedges.end()) {
			edges.push_back(ele);
			sedges.insert(hashc);
		}
		if (!graph.directed) {
			ele.row = cid;
			ele.col = rid;

			if (cid > maxid) {
				maxid = cid;
			}

			hashc = cid*hashroot + rid;
			if (sedges.find(hashc) == sedges.end()) {
				edges.push_back(ele);
				sedges.insert(hashc);
			}
		}

	}

	if (self) {
		for (uint i = 0; i <= maxid; i++) {
			mtuple ele;
			ele.row = i;
			ele.col = i;
			ele.w = 1;
			edges.push_back(ele);
		}
	}
	sort(edges.begin(), edges.end());
	if (edges.size() == 0) {
		cout << title << "\tedge size zero" << endl;
		return;
	}
	graph.nodeSize = maxid;
	cout << graph.nodeSize << endl;
	graph.edgeSize = edges.size();
	graph.nbs = new unsigned int[graph.edgeSize];
	graph.indexs = new unsigned int[(graph.nodeSize + 2)];
	if (graph.weighted) {
		graph.ws = new double[graph.edgeSize];
		graph.volumns = new double[graph.nodeSize + 1];
	}
	else {
		graph.ws = NULL;
		graph.volumns = NULL;
	}

	if (DEBUG_CREATE_GRAPH)
		cout << "graph.nodeSize\t" << graph.nodeSize << endl;

	graph.indexs[0] = 0;
	int old = 0;
	if (!graph.weighted) {
		for (unsigned int i = 0; i < graph.edgeSize; i++) {
			graph.nbs[i] = edges[i].col;
			if (edges[i].row != old) {
				for (unsigned int next = old + 1; next <= edges[i].row; next++) {
					graph.indexs[next] = i;
				}
			}
			old = edges[i].row;
		}
	}
	else {
		double volumn = 0;
		for (unsigned int i = 0; i < graph.edgeSize; i++) {
			graph.nbs[i] = edges[i].col;
			graph.ws[i] = edges[i].w;

			if (edges[i].row != old) {
				graph.volumns[old] = volumn;
				volumn = 0;
				for (unsigned int next = old + 1; next <= edges[i].row; next++) {
					graph.indexs[next] = i;
					graph.volumns[next] = 0;
				}
			}
			volumn += edges[i].w;
			old = edges[i].row;
		}
		graph.volumns[maxid] = volumn;
	}

	for (unsigned int next = old + 1; next <= maxid + 1; next++) {
		graph.indexs[next] = graph.edgeSize;
	}

	if (DEBUG_CREATE_GRAPH) {
		cout << "done\t" << title << "\r";
	}
}

void createGraph(const sparserow & from, sparserow & to, bool reverse) {
	if (!reverse) {
		cout << "maybe wrong here" << endl;
		return;
	}
	vector<mtuple> edges;
	to.directed = from.directed;
	to.edgeSize = from.edgeSize;
	to.weighted = from.weighted;
	uint maxid = 0;
	for (uint i = 0; i <= from.nodeSize; i++) {
		uint begin = from.indexs[i];
		uint end = from.indexs[i + 1];
		for (uint j = begin; j < end; j++) {
			uint id2 = from.nbs[j];
			mtuple ele;
			ele.row = id2;
			ele.col = i;
			if (from.weighted) {
				double weight = from.ws[j];
				ele.w = weight;
			}

			if (ele.row > maxid)
				maxid = ele.row;
			edges.push_back(ele);
		}
	}
	sort(edges.begin(), edges.end());
	if (edges.size() == 0) {
		cout << "\tedge size zero" << endl;
		return;
	}

	to.nodeSize = maxid;
	to.nbs = new uint[to.edgeSize];
	to.indexs = new uint[to.nodeSize + 2];
	//cout << "node size " << to.nodeSize << endl;
	if (to.weighted) {
		to.ws = new double[to.edgeSize];
		to.volumns = new double[to.nodeSize + 1];
	}
	else {
		to.ws = NULL;
		to.volumns = NULL;
	}

	to.indexs[0] = 0;
	int old = 0;
	if (!to.weighted) {
		for (unsigned int i = 0; i < to.edgeSize; i++) {
			to.nbs[i] = edges[i].col;
			if (edges[i].row != old) {
				for (unsigned int next = old + 1; next <= edges[i].row; next++) {
					to.indexs[next] = i;
				}
			}
			old = edges[i].row;
		}
	}
	else {
		double volumn = 0;
		for (unsigned int i = 0; i < to.edgeSize; i++) {
			to.nbs[i] = edges[i].col;
			to.ws[i] = edges[i].w;

			if (edges[i].row != old) {
				to.volumns[old] = volumn;
				volumn = 0;
				for (unsigned int next = old + 1; next <= edges[i].row; next++) {
					to.indexs[next] = i;
					to.volumns[next] = 0;
				}
			}
			volumn += edges[i].w;
			old = edges[i].row;
		}
		to.volumns[maxid] = volumn;
	}

	for (unsigned int next = old + 1; next <= maxid + 1; next++) {
		to.indexs[next] = to.edgeSize;
	}
	//to.indexs[to.nodeSize] = to.edgeSize;
}

void sparserow::check() {
	if (!directed) {
		for (uint i = 0; i < edgeSize; i++) {
			if (nbs[i] > nodeSize) {
				cout << "graph check Error\t" << i << "\t" << nbs[i] << endl;
				system("pause");
			}
		}
	}

	for (uint i = 0; i < nodeSize; i++) {
		if (indexs[i] > edgeSize || indexs[i]>indexs[i + 1]) {
			cout << "graph indexs check Error\t" << i << "\t" << indexs[i] << "\tsize "<< edgeSize<<"\ti+1\t"<< indexs[i + 1]<<endl;
			exit(1);
		}
	}
	if (weighted) {
		for (uint i = 0; i < nodeSize; i++) {
			if (indexs[i + 1] > indexs[i] && volumns[i] == 0) {
				cout << "weighted graph check Error\t" << i << "\t" << indexs[i] << "\t" << indexs[i + 1] << endl;
				system("pause");
			}
		}
	}
	return;
}

/*
* res = res + alpha*x;
*/
void addArray(map<uint, double>& res, const map<uint, double>& x, double alpha) {
	map<uint, double>::const_iterator mite = x.begin();
	map<uint, double>::iterator index;

	for (; mite != x.end(); mite++) {
		index = res.find(mite->first);
		if (index == res.end()) {
			res.insert(map<int, float>::value_type(mite->first, alpha*mite->second));
		}
		else
			index->second += alpha*mite->second;
	}
}

/*
* input: x and alpha
* output:alpha*P*x- through pointer
* calculate alpha*P*x
*/
void multi(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double alpha, bool self) {

	uint *nbs = P.nbs;
	double * volumns = P.volumns;
	double *ws = P.ws;
	uint *indexs = P.indexs;
	res.clear();
	if (alpha < 1e-6)
		return;
	if (P.weighted) {
		for (map<uint, double>::const_iterator mite = x.begin(); mite != x.end(); mite++) {
			if (mite->first > P.nodeSize)
				continue;
			int cnodeid = mite->first;
			if (indexs[cnodeid] == indexs[cnodeid + 1]) {
				if (self)
					res[cnodeid] += alpha*mite->second;
			}
			else
				for (int i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
					uint nb = nbs[i];
					if (res.find(nb) == res.end()) {
						res[nb] = alpha*mite->second*ws[i] / volumns[cnodeid];
					}
					else
						res[nb] += alpha*mite->second*ws[i] / volumns[cnodeid];
				}
		}
		return;
	}


	//unweighted case
	for (map<uint, double>::const_iterator mite = x.begin(); mite != x.end(); mite++) {
		uint cnodeid = mite->first;
		if (cnodeid > P.nodeSize) {
			continue;
		}

		if (indexs[cnodeid] == indexs[cnodeid + 1]) {
			if (self)
				res[cnodeid] += alpha*mite->second;
		}
		else
			for (uint i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
				uint size = indexs[cnodeid + 1] - indexs[cnodeid];
				uint nb = nbs[i];
				if (res.find(nb) == res.end()) {
					res[nb] = alpha*mite->second / size;
				}
				else
					res[nb] += alpha*mite->second / size;
			}
	}
}


/*
* input: x and alpha
* output:alpha*P*x- through pointer
* calculate alpha*P*x
*/
void multi(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double alpha) {
	double * volumns = P.volumns;
	uint *nbs = P.nbs;
	double *ws = P.ws;
	uint *indexs = P.indexs;
	res.clear();
	if (alpha < 1e-6)
		return;
	if (P.weighted) {
		for (map<uint, double>::const_iterator mite = x.begin(); mite != x.end(); mite++) {
			if (mite->first > P.nodeSize)
				continue;
			int cnodeid = mite->first;
			if (indexs[cnodeid] == indexs[cnodeid + 1]) {

				res[cnodeid] += alpha*mite->second;
			}
			else
				for (int i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
					uint nb = nbs[i];
					if (res.find(nb) == res.end()) {
						res[nb] = alpha*mite->second*ws[i] / volumns[cnodeid];
					}
					else
						res[nb] += alpha*mite->second*ws[i] / volumns[cnodeid];
				}
		}
		return;
	}


	//unweighted case
	for (map<uint, double>::const_iterator mite = x.begin(); mite != x.end(); mite++) {
		uint cnodeid = mite->first;
		if (cnodeid > P.nodeSize) {
			continue;
		}

		if (indexs[cnodeid] == indexs[cnodeid + 1]) {
			res[cnodeid] += alpha*mite->second;
		}
		else
			for (uint i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
				uint size = indexs[cnodeid + 1] - indexs[cnodeid];
				uint nb = nbs[i];
				if (res.find(nb) == res.end()) {
					res[nb] = alpha*mite->second / size;
				}
				else
					res[nb] += alpha*mite->second / size;
			}
	}
}


/*
* input: x and alpha, theta
* x = x0+x1
* output:alpha*(P*x0+x1)- through pointer
*/

void amulti(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double alpha, double theta, uint seed, bool *visited, bool self) {
	uint *nbs = P.nbs;
	double *ws = P.ws;
	double sum = 0;
	double * volumns = P.volumns;
	uint *indexs = P.indexs;
	res.clear();
	if (alpha < 1e-6)
		return;
	queue<uint> que;
	memset(visited, 0, sizeof(bool)*(P.nodeSize + 1));
	que.push(seed);
	visited[seed] = true;
	if (P.weighted) {
		while (sum < theta && !que.empty()) {
			uint cnodeid = que.front();
			que.pop();
			if (cnodeid > P.nodeSize) {
				continue;
			}
			map<uint, double>::const_iterator mite = x.find(cnodeid);
			if (mite == x.end())
				continue;
			if (indexs[cnodeid] == indexs[cnodeid + 1]) {
				if (self)
					res[cnodeid] += alpha*mite->second;
			}
			else
				for (uint i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
					uint size = indexs[cnodeid + 1] - indexs[cnodeid];
					uint nb = nbs[i];
					if (nb<P.nodeSize + 1 && !visited[nb]) {
						que.push(nb);
						visited[nb] = true;
					}
					if (res.find(nb) == res.end()) {
						res[nb] = alpha*mite->second*ws[i] / volumns[cnodeid];
					}
					else
						res[nb] += alpha*mite->second*ws[i] / volumns[cnodeid];
				}
			sum += mite->second;
		}
		res[seed] += 1 - sum;
		return;
	}

	//unweighted case
	while (sum < theta && !que.empty()) {
		uint cnodeid = que.front();
		que.pop();
		if (cnodeid > P.nodeSize) {
			//cout << "ERROR Approximate Calculation" << endl;
			continue;
		}
		map<uint, double>::const_iterator mite = x.find(cnodeid);
		if (mite == x.end())
			continue;
		if (indexs[cnodeid] == indexs[cnodeid + 1]) {
			if (self)
				res[cnodeid] += alpha*mite->second;
		}
		else
			for (uint i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
				uint size = indexs[cnodeid + 1] - indexs[cnodeid];
				uint nb = nbs[i];
				if (nb<P.nodeSize + 1 && !visited[nb]) {
					que.push(nb);
					visited[nb] = true;
				}
				if (res.find(nb) == res.end()) {
					res[nb] = alpha*mite->second / size;
				}
				else
					res[nb] += alpha*mite->second / size;
			}
		sum += mite->second;
	}

	if (1 - sum>1e-6)
		res[seed] += 1 - sum;

}

double extractPart(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double theta, uint seed, bool *visited) {
	uint *nbs = P.nbs;
	double sum = 0;
	uint *indexs = P.indexs;
	res.clear();

	queue<uint> que;
	memset(visited, 0, sizeof(bool)*(P.nodeSize + 1));
	que.push(seed);
	visited[seed] = true;

	//unweighted case
	while (sum < theta && !que.empty()) {
		uint cnodeid = que.front();
		que.pop();
		if (cnodeid > P.nodeSize) {
			continue;
		}
		map<uint, double>::const_iterator mite = x.find(cnodeid);
		if (mite == x.end())
			continue;
		if (indexs[cnodeid] == indexs[cnodeid + 1]) {
			continue;
		}
		else
			for (uint i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
				uint size = indexs[cnodeid + 1] - indexs[cnodeid];
				uint nb = nbs[i];
				if (nb<P.nodeSize + 1 && !visited[nb]) {
					que.push(nb);
					visited[nb] = true;
				}
			}
		res[cnodeid] = mite->second;
		sum += mite->second;
	}
	if(sum<theta)
		for (map<uint, double>::const_iterator mite = x.begin(); mite != x.end(); mite++) {
			res[mite->first] = mite->second;
			sum += mite->second;
			if (sum >= theta)
				break;
		}
	return sum;
}

void calConductance(const sparserow &graph, vector<double> &conductances, const vector<int_double_pair_sorted_by_score> tep, const sparserow &P) {
	int size = tep.size();
	if (size <= 0)
		return;

	set<uint> nodes;
	for (uint ti = 0; ti < tep.size(); ti++) {
		uint id = tep[ti].id;
		if (id > P.nodeSize) {
			if(conductances.size()==0)
				conductances.push_back(1);
			else
				conductances.push_back(conductances[conductances.size() - 1]);
			continue;
		}
		for (uint oindex = P.indexs[id]; oindex < P.indexs[id + 1]; oindex++) {
			uint oid = P.nbs[oindex];
			if(oid<=graph.nodeSize)
				nodes.insert(oid);
		}
		conductances.push_back(calConductance(graph, nodes));
	}

}

double calConductance(const sparserow &graph,const set<uint> &cmty){
	if (cmty.size() <= 0)
		return 1;

	//init
	//S = {indexs[0]}
	//hS = V-hS[0];
	uint totaldegree = graph.edgeSize*2;
	uint as = 0;
	uint fs = 0;
	for (set<uint>::const_iterator site = cmty.begin(); site != cmty.end();site++) {
		uint node = *site;
		as += graph.degree(node);
		for (int nbindex = graph.indexs[node]; nbindex < graph.indexs[node + 1]; nbindex++) {
			uint nb = graph.nbs[nbindex];
			if (cmty.find(nb) == cmty.end())
				fs++;
		}
	}
	if (as == 0)
		return 1;
	if (as * 2 < totaldegree)
		return 1.0*fs / as;
	return  1.0*fs / (totaldegree - as);

}

double calConductance(const sparserow &graph, const vector<uint> &cmty) {
	set<uint> cmt;
	for (auto c : cmty) {
		cmt.insert(c);
	}
	return calConductance(graph, cmt);
}



void calConductance(const sparserow &graph, vector<double> &conductances, const vector<int_double_pair_sorted_by_score> tep) {
	int size = tep.size();
	if (size <= 0)
		return;

	//init
	//S = {indexs[0]}
	//hS = V-hS[0];
	int ahs = graph.indexs[graph.nodeSize + 1];

	int as = graph.degree(tep[0].id);
	ahs = ahs - as;
	int s2hs = as;
	int bestindex = 1;
	set<uint> already;
	already.insert(tep[0].id);
	conductances.push_back(0);
	for (int i = 1; i < size && i<max_cmty_size; i++) {

		if (tep[i].id > graph.nodeSize) {
			conductances.push_back(conductances[conductances.size() - 1]);
		}

		as += graph.degree(tep[i].id);
		ahs -= graph.degree(tep[i].id);

		int nbIns = 0;

		for (uint k = graph.indexs[tep[i].id]; k < graph.indexs[tep[i].id + 1]; k++) {
			if (already.find(graph.nbs[k]) != already.end()) {
				nbIns++;
			}
		}

		s2hs -= nbIns;


		// add edge from index[j] to ahs

		s2hs += graph.degree(tep[i].id) - nbIns;
		already.insert(tep[i].id);
		float tempscore = 1.0*s2hs / min(as, ahs);
		conductances.push_back(tempscore);
	}
}

/*
void createGraph(const string path, sparserow & graph) {
	uint hashroot = 39569;
	set<uint> sedges;
	sedges.insert(0);
	ifstream fin(path);
	string line;
	vector<mtuple> edges;
	int maxid = -1;
	while(!fin.eof()) {
		getline(fin, line);
		if (line.length() < 2 || line.c_str()[0]=='#')
			continue;
		vector<string> ss;
		split(ss, line, '\t');
		int rid = atoi(ss[0].c_str());
		int cid = atoi(ss[1].c_str());
		if (cid > maxid) {
			maxid = cid;
		}

		if (rid > maxid) {
			maxid = rid;
		}

		double w = 1.0;
		if (ss.size() >= 3) {
			w = atof(ss[2].c_str());
		}
		mtuple ele;
		ele.row = rid;
		ele.col = cid;
		ele.w = w;
		uint hashc = rid*hashroot + cid;
		if (sedges.find(hashc) == sedges.end()) {
			edges.push_back(ele);
			sedges.insert(hashc);
		}
		ele.row = cid;
		ele.col = rid;

		hashc = cid*hashroot + rid;
		if (sedges.find(hashc) == sedges.end()) {
			edges.push_back(ele);
			sedges.insert(hashc);
		}

	}
	
	sort(edges.begin(), edges.end());

	fin.close();
	
	graph.nodeSize = maxid;
	graph.edgeSize = edges.size();
	graph.nbs = (unsigned int*)malloc(sizeof(unsigned int)*graph.edgeSize);
	graph.indexs = (unsigned int*)malloc(sizeof(unsigned int)*(graph.nodeSize+2));
	//graph.ws = (double*)malloc(sizeof(double)*graph.edgeSize);

	cout << graph.nodeSize << endl;


	graph.indexs[0] = 0;
	int old = 0;
	for (size_t i = 0; i < graph.edgeSize; i++) {
		graph.nbs[i] = edges[i].row;
		//graph.ws[i] = edges[i].w;
		if (edges[i].col != old) {
			for (unsigned int next = old + 1; next <= edges[i].col; next++) {
				graph.indexs[next] = i;
			}
		}
		old = edges[i].col;
	}
	graph.indexs[maxid + 1] = graph.edgeSize;
}
void createGraph(ifstream &fin, sparserow & graph, char endchar) {
	createGraph(fin, graph, endchar, false);
}

// Create graph from file stream. 
// skip the lines startswith '#'
// end with the lines startswith '-'
void createGraph(ifstream &fin, sparserow & graph, char endchar, bool self) {
	uint hashroot = 39569;
	set<uint> sedges;
	string line="";
	vector<mtuple> edges;
	int maxid = -1;
	string title;
	graph.weighted = false;
	while (!fin.eof()) {
		getline(fin, line);
		if (line.c_str()[0] == endchar) {
			if(DEBUG_CREATE_GRAPH)
				cout << "begin to read " << line << endl;
			title = line;
			break;
		}
		if (line.length() < 2 || line.c_str()[0] == '#')
			continue;
		vector<string> ss;
		split(ss, line, '\t');
		if (ss.size() <2) {
			cout << "skip :"<<line << endl;
			continue;
		}
		int rid = atoi(ss[0].c_str());
		int cid = atoi(ss[1].c_str());
		if (cid > maxid) {
			maxid = cid;
		}

		double w = 1.0;
		if (ss.size() >= 3) {
			w = atof(ss[2].c_str());
			graph.weighted = true;
		}
		mtuple ele;
		ele.row = rid;
		ele.col = cid;
		ele.w = w;
		uint hashc = rid*hashroot + cid;
		if (sedges.find(hashc) == sedges.end()) {
			edges.push_back(ele);
			sedges.insert(hashc);
		}
		if (!graph.directed) {
			ele.row = cid;
			ele.col = rid;

			if (rid > maxid) {
				maxid = rid;
			}

			hashc = cid*hashroot + rid;
			if (sedges.find(hashc) == sedges.end()) {
				edges.push_back(ele);
				sedges.insert(hashc);
			}
		}

	}

	if (self) {
		for (uint i = 0; i <= maxid; i++) {
			mtuple ele;
			ele.row = i;
			ele.col = i;
			ele.w = 1;
			edges.push_back(ele);
		}
	}
	sort(edges.begin(), edges.end());
	if (edges.size() == 0) {
		cout << title<<"\tedge size zero" << endl;
		return;
	}
	graph.nodeSize = maxid;
	graph.edgeSize = edges.size();
	graph.nbs = (unsigned int*)malloc(sizeof(unsigned int)*graph.edgeSize);
	graph.indexs = (unsigned int*)malloc(sizeof(unsigned int)*(graph.nodeSize+2));
	if (graph.weighted) {
		graph.ws = (double*)malloc(sizeof(double)*(graph.edgeSize));
		graph.volumns = (double*)malloc(sizeof(double)*graph.nodeSize+1);
	}
	else {
		graph.ws = NULL;
		graph.volumns = NULL;
	}

	if (DEBUG_CREATE_GRAPH)
		cout << "graph.nodeSize\t"<<graph.nodeSize << endl;

	graph.indexs[0] = 0;
	int old = 0;
	if (!graph.weighted) {
		for (unsigned int i = 0; i < graph.edgeSize; i++) {
			graph.nbs[i] = edges[i].row;
			if (edges[i].col != old) {
				for (unsigned int next = old + 1; next <= edges[i].col; next++) {
					graph.indexs[next] = i;
				}
			}
			old = edges[i].col;
		}
	}
	else {
		double volumn = 0;
		for (unsigned int i = 0; i < graph.edgeSize; i++) {
			graph.nbs[i] = edges[i].row;
			graph.ws[i] = edges[i].w;
			volumn += edges[i].w;
			if (edges[i].col != old) {
				graph.volumns[old] = volumn;
				volumn = 0;
				for (unsigned int next = old + 1; next <= edges[i].col; next++) {
					graph.indexs[next] = i;
				}
			}
			old = edges[i].col;
		}
		graph.volumns[maxid] = volumn;
	}

	graph.indexs[maxid + 1] = graph.edgeSize;
	
	if (DEBUG_CREATE_GRAPH) {
		cout << "done\t" << title << "\r";
	}


}

void createGraph(const sparserow & from, sparserow & to, bool reverse) {
	if (!reverse) {
		cout << "maybe wrong here" << endl;
		return;
	}
	vector<mtuple> edges;
	to.directed = from.directed;
	to.weighted = from.weighted;
	to.edgeSize = from.edgeSize;
	uint maxid = 0;
	for (uint i = 0; i <= from.nodeSize; i++) {
		uint begin = from.indexs[i];
		uint end = from.indexs[i + 1];
		for (uint j = begin; j < end; j++) {
			uint id2 = from.nbs[j];
			mtuple ele;
			ele.row = i;
			ele.col = id2;
			if (from.weighted) {
				ele.w = from.ws[j];
			}
			if (ele.col > maxid)
				maxid = ele.col;
			edges.push_back(ele);
		}
	}
	sort(edges.begin(), edges.end());
	if (edges.size() == 0) {
		cout << "\tedge size zero" << endl;
		return;
	}

	to.nodeSize = maxid;
	to.nbs = new uint[to.edgeSize];
	to.indexs = new uint[to.nodeSize + 2];
	cout << "node size " << to.nodeSize << endl;
	to.indexs[0] = 0;
	int old = 0;
	cout << "edges size\t" << edges.size() << endl;
	for (unsigned int i = 0; i < to.edgeSize; i++) {
		to.nbs[i] = edges[i].row;
		if (to.weighted) {
			to.ws[i] = edges[i].w;
		}
		if (edges[i].col != old) {
			for (unsigned int next = old + 1; next <= edges[i].col; next++) {
				to.indexs[next] = i;
			}
		}
		old = edges[i].col;
	}

	to.indexs[to.nodeSize+1] = to.edgeSize;
}

*/
/*

unsigned int calTriangle(int i, int j, sparserow & graph) {
unsigned int count = 0;
unsigned int *indexs = graph.indexs;
unsigned int *nbs = graph.nbs;

size_t begin1 = indexs[i], begin2 = indexs[j];
size_t end1 = indexs[i + 1],end2 = indexs[j+1];
size_t t1 = begin1;
size_t  t2 = begin2;
while (t1 != end1 && t2 != end2) {
if (nbs[t1] == nbs[t2]) {
count++;
t1++;
t2++;
}
else if (nbs[t1] < nbs[t2]) {
t1++;
}
else
t2++;
}

return count;
}



unsigned int calRectangle(int i, int j, sparserow & graph) {
unsigned int count = 0;
unsigned int *indexs = graph.indexs;
unsigned int *nbs = graph.nbs;
size_t begin1 = indexs[i], begin2 = indexs[j];
size_t end1 = indexs[i + 1], end2 = indexs[j + 1];
size_t t1 = begin1;
size_t  t2 = begin2;
vector<int> vnbs;
while (t1 != end1 && t2 != end2) {
if (nbs[t1] == nbs[t2]) {
t1++;
t2++;
vnbs.push_back(nbs[t1]);
}
else if (nbs[t1] < nbs[t2]) {
t1++;
}
else
t2++;
}
if (vnbs.size() <= 1)
return 0;
for (int nbindex1 = 0; nbindex1 < vnbs.size(); nbindex1++) {
for (int nbindex2 = nbindex1 + 1; nbindex2 < vnbs.size(); nbindex2++) {
unsigned int nb1 = vnbs[nbindex1];
unsigned int nbbegin = indexs[nb1];
unsigned int nbend = indexs[nb1 + 1];


for (int tempi = nbbegin; tempi != nbend; tempi++) {
if (nbs[tempi] == vnbs[nbindex2]) {
count += 1;
break;
}
}
}
}
return count;
}


bool isTriangle(int node1, int node2, int node3, sparserow & graph) {
unsigned int *ai = graph.indexs;
unsigned int *aj = graph.nbs;

unsigned int nbbegin = ai[node1];
unsigned int nbend = ai[node1 + 1];

int count = 0;
for (int tempi = nbbegin; tempi != nbend; tempi++) {
if (aj[tempi] == node2 || aj[tempi] == node3) {
count += 1;
if (count == 2)
break;
}
}

if (count != 2)
return false;

for (int tempi = ai[node2]; tempi != ai[node2+1]; tempi++) {
if (aj[tempi] == node3) {
return true;
}
}
return false;

}



unsigned int calFiveangle(int i, int j, sparserow & graph) {
unsigned int count = 0;
unsigned int *ai = graph.indexs;
unsigned int *aj = graph.nbs;
size_t begin1 = ai[i], begin2 = ai[j];
size_t end1 = ai[i + 1], end2 = ai[j + 1];
size_t t1 = begin1;
size_t  t2 = begin2;
vector<int> nbs;
while (t1 != end1 && t2 != end2) {
if (aj[t1] == aj[t2]) {
t1++;
t2++;
nbs.push_back(aj[t1]);
}
else if (aj[t1] < aj[t2]) {
t1++;
}
else
t2++;
}
if (nbs.size() <= 2)
return 0;
for (int nbindex1 = 0; nbindex1 < nbs.size(); nbindex1++) {
unsigned int nb1 = nbs[nbindex1];
for (int nbindex2 = nbindex1 + 1; nbindex2 < nbs.size(); nbindex2++) {
unsigned int nb2 = nbs[nbindex2];
for (int nbindex3 = nbindex2 + 1; nbindex3 < nbs.size(); nbindex3++) {
unsigned int nb3 = nbs[nbindex3];
if (isTriangle(nb1, nb2, nb3, graph)) {
count += 1;
break;
}
}
}
}
return count;
}



//count the number of triangles for each pair of nodes

void createTriGraph(sparserow & graph, sparserow & Trigraph, const string path) {

	unsigned int edgeSize = graph.edgeSize;
	unsigned int nodeSize = graph.nodeSize;

	unsigned int * indexs = graph.indexs;
	unsigned int * edges = graph.nbs;

	vector<mtuple> triedges;
	int count = 0;
	for (unsigned int i = 0; i < nodeSize; i++) {
		for (unsigned int jindex = indexs[i]; jindex < indexs[i + 1]; jindex++) {
			unsigned int j = edges[jindex];
			unsigned int trc = calTriangle(i, j, graph);
			if (trc == 0)
				continue;
			mtuple ele;
			ele.row = i;
			ele.col = j;
			ele.w = trc;
			triedges.push_back(ele);
		}
	}

	sort(triedges.begin(), triedges.end());

	ofstream fout(path);
	for (vector<mtuple>::iterator vite = triedges.begin(); vite != triedges.end(); vite++) {
		fout << vite->row << "\t" << vite->col << "\t" << vite->w << endl;
	}
	fout.close();

}



* count the number of triangles for each pair of nodes

void createFourGraph(sparserow & graph, sparserow & Trigraph, const string path) {

	unsigned int edgeSize = graph.edgeSize;
	unsigned int nodeSize = graph.nodeSize;

	unsigned int * indexs = graph.indexs;
	unsigned int * edges = graph.nbs;

	vector<mtuple> triedges;
	int count = 0;
	for (unsigned int i = 0; i < nodeSize; i++) {
		for (unsigned int jindex = indexs[i]; jindex < indexs[i + 1]; jindex++) {
			unsigned int j = edges[jindex];
			int trc = calRectangle(i, j, graph);
			if (trc == 0)
				continue;
			mtuple ele;
			ele.row = i;
			ele.col = j;
			ele.w = trc;
			triedges.push_back(ele);
		}
	}

	sort(triedges.begin(), triedges.end());

	ofstream fout(path);
	for (vector<mtuple>::iterator vite = triedges.begin(); vite != triedges.end(); vite++) {
		fout << vite->row << "\t" << vite->col << "\t" << vite->w << endl;
	}
	fout.close();

}


* count the number of triangles for each pair of nodes
void createFiveGraph(sparserow & graph, sparserow & Trigraph, const string path) {

	unsigned int m = graph.edgeSize;
	unsigned int n = graph.nodeSize;

	unsigned int * indexs = graph.indexs;
	unsigned int * edges = graph.nbs;

	vector<mtuple> triedges;
	int count = 0;
	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int jindex = indexs[i]; jindex < indexs[i + 1]; jindex++) {
			unsigned int j = edges[jindex];
			int trc = calFiveangle(i, j, graph);
			if (trc == 0)
				continue;
			mtuple ele;
			ele.row = i;
			ele.col = j;
			ele.w = trc;
			triedges.push_back(ele);
		}
	}

	sort(triedges.begin(), triedges.end());

	ofstream fout(path);
	for (vector<mtuple>::iterator vite = triedges.begin(); vite != triedges.end(); vite++) {
		fout << vite->row << "\t" << vite->col << "\t" << vite->w << endl;
	}
	fout.close();

}

* input: x and alpha
* output:alpha*P*x- through pointer
* calculate alpha*P*x
void multi(map<uint, double>& res, const sparserow& P, const map<uint, double>& x, double alpha) {
	uint *nbs = P.nbs;
	uint *indexs = P.indexs;
	res.clear();
	if (alpha < 1e-6)
		return;
	if (P.weighted) {
		double *ws = P.ws;
		double *volumns = P.volumns;
		for (map<uint, double>::const_iterator mite = x.begin(); mite != x.end(); mite++) {
			int cnodeid = mite->first;

			for (int i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
				uint nb = nbs[i];
				if (res.find(nb) == res.end()) {
					res[nb] = alpha*mite->second* ws[i] / volumns[cnodeid];
				}
				else
					res[nb] += alpha*mite->second* ws[i] / volumns[cnodeid];
			}
		}
		return;
	}
	//unweighted case
	for (map<uint, double>::const_iterator mite = x.begin(); mite != x.end(); mite++) {
		uint cnodeid = mite->first;
		if (cnodeid > P.nodeSize) {
			continue;
		}
		for (uint i = indexs[cnodeid]; i < indexs[cnodeid + 1]; i++) {
			uint size = indexs[cnodeid + 1] - indexs[cnodeid];
			uint nb = nbs[i];
			if (res.find(nb) == res.end()) {
				res[nb] = alpha*mite->second / size;
			}
			else
				res[nb] += alpha*mite->second / size;
		}
	}
}


void sparserow::check() {
	if (!directed) {
		for (uint i = 0; i < edgeSize; i++) {
			if (nbs[i] > nodeSize) {
				cout << "graph check Error\t" << i << "\t" << nbs[i] << endl;
			}
		}
	}
	for (uint i = 0; i <= nodeSize; i++) {
		if (indexs[i] > edgeSize) {
			cout << "graph indexs check Error\t" << i << "\t" << nbs[i] << endl;
		}
	}

	return;
}


* res = res + alpha*x;

void addArray(map<uint, double>& res, const map<uint, double>& x, double alpha) {
	map<uint, double>::const_iterator mite = x.begin();
	map<uint, double>::iterator index;

	for (; mite != x.end(); mite++) {
		index = res.find(mite->first);
		if (index == res.end()) {
			res.insert(map<int, float>::value_type(mite->first, alpha*mite->second));
		}
		else
			index->second += alpha*mite->second;
	}
}

*/