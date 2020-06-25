#include"stdafx.h"
#include"tools.h"
#include"MultiplexGraph.h"
#include"Similarities.h"
#define MAXITE 100
extern bool USETOPK;
extern double alpha;
extern double epsilon;
extern ofstream *fLog;
using namespace std;
extern uint topk;
extern int itenumber;
#define and &&
#define or ||
#define DEBUG true
extern double lamda;
extern Wmodel wmodel;
extern double selfsimi;
extern double epsilon_matrix, theta;
extern int min_cmty_size, max_cmty_size;

extern bool NeedFindCommunity, MC, phase1Only, SweepNormByDegree;
MultiplexGraph::MultiplexGraph() {
	nmbofLayer = 0;
	graphs = NULL;
	mute = false;
	xs = new map<uint, double>[nmbofLayer];
	distribution = geometric_distribution<int>(alpha);
}

void MultiplexGraph::setMute(bool m) {
	mute = m;
}

/*
* create the multiple graph from a single file. 
*/
MultiplexGraph::MultiplexGraph(string path) {
	distribution = geometric_distribution<int>(1-alpha);

	nodesize = 0;
	ifstream fin(path);
	if (!fin) {
		cout << path << "not found" << endl;
		exit(1);
	}
	string temp;
	getline(fin, temp);
	getline(fin, temp);
	uint index_head = 0;
	int t_nmbofLayer = getUint(temp, index_head);
	if (nmbofLayer <=0) {
		nmbofLayer = t_nmbofLayer;
	}
	getline(fin, temp);
	index_head = 0;
	bool directed = getUint(temp, index_head) == 1 ? true : false;


	graphs = (sparserow*)malloc(sizeof(sparserow)*nmbofLayer);
	uint i = 0;
	while (true) {
		getline(fin, temp);
		if (temp.c_str()[0] == '-')
			break;
	}
	if (DEBUG_CREATE_GRAPH)
		cout << "begin to read " << temp << endl;
	for (; i < nmbofLayer; i++) {
		if (DEBUG_CREATE_GRAPH)
			cout << "begin to read #############" << i <<"\r";
		sparserow graph;
		graph.directed = directed;
		graphs[i] = graph;
		createGraph(fin,graphs[i],'-');
		if (nodesize < graphs[i].nodeSize)
			nodesize = graphs[i].nodeSize;
		if (DEBUG_CREATE_GRAPH)
			cout << "read layer \t #############" << i << "done\r";
	}
	fin.close();
	visited = new bool[nodesize + 1];
}


uint MultiplexGraph::getNmbofLayer() {
	return nmbofLayer;
}

void MultiplexGraph::setNmbofLayer(uint m) {
	nmbofLayer = m;
}

MultiplexGraph::~MultiplexGraph() {

	//for (int i = 0; i < nmbofLayer; i++) {
	//	delete[] graphs[i].indexs;
	//	delete[] graphs[i].nbs;
	//	delete[] graphs[i].ws;
	//}

	//if(graphs!=NULL){
	//	delete[] graphs;
	//	graphs = NULL;
	//}
	//cout << "delete" << endl;
	//if (xs != NULL) {
	//	delete[] xs;
	//	xs = NULL;
	//}
	cout << "delete done" << endl;
}

MultiplexGraph::MultiplexGraph(string path, int layer) {
	distribution = geometric_distribution<int>(1 - alpha);

	nodesize = 0;
	ifstream fin(path);
	if (!fin) {
		cout << path << "not found" << endl;
		exit(1);
	}
	string temp;
	getline(fin, temp);
	getline(fin, temp);
	uint index_head = 0;
	nmbofLayer = getUint(temp, index_head);
	if (layer > 0)
		nmbofLayer = layer;
	getline(fin, temp);
	index_head = 0;
	bool directed = getUint(temp, index_head) == 1 ? true : false;


	graphs = (sparserow*)malloc(sizeof(sparserow)*nmbofLayer);
	uint i = 0;
	while (true) {
		getline(fin, temp);
		if (temp.c_str()[0] == '-')
			break;
	}
	for (; i < nmbofLayer; i++) {
		cout << "read file " << nmbofLayer << "\r";
		sparserow graph;
		graph.directed = directed;
		graphs[i] = graph;
		createGraph(fin, graphs[i], '-');
		if (nodesize < graphs[i].nodeSize)
			nodesize = graphs[i].nodeSize;
	}
	fin.close();
	visited = new bool[nodesize + 1];
	cout << "read file done" << endl;
}


/*
* MultiWalker Model on Multiple Network
 one walker for each layer.
 return vector are determined by influent nodes of each layer
*/
void MultiplexGraph::basicMultiWalker(uint seed, vector<uint> &cmty, map<uint,double> &prxmts) {
	
	delete[] xs;
	xs = new map<uint, double>[nmbofLayer];
	//initially, all influenced nodes are seed nodes
	uint *inodes = new uint[nmbofLayer];
	memset(inodes, seed, nmbofLayer);

	uint *ninodes = new uint[nmbofLayer];
	map<uint, double> v, nv;
	v[seed] = 1.0;
	uint t = 0;
	map<uint, double> *nxs = new map<uint, double>[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		xs[i].insert(map<int, double>::value_type(seed, 1.0));
	}
	int count = 0;
	for(;count<100;count++) {
		nv.clear();
		for (uint i = 0; i < nmbofLayer; i++) {
			cout << "i\t\t\t"<< count<<"\t" <<i<< "\r";
			multi(nxs[i], graphs[i], xs[i], alpha);
			addArray(nxs[i], v, (1 - alpha));
			norm(nxs[i]);
			uint mindex = 0;
			double maxv = -1;
			for (map<uint, double>::iterator mite = nxs[i].begin(); mite != nxs[i].end(); mite++) {

				if (mite->second > maxv) {
					maxv = mite->second;
					mindex = mite->first;
				}
			}

			if (nv.find(mindex) != nv.end()) {
				nv[mindex] += 1.0 / nmbofLayer;
			}
			else
				nv[mindex] = 1.0 / nmbofLayer;
		}
		double diff = mabs(nxs, xs, nmbofLayer);
		//if (DEBUG)
		//	cout << diff << endl;
		if ( diff< epsilon) {
			copyArray(prxmts, *xs);
			break;
		}
		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;
		v.swap(nv);
	}

	cout << "begin to find cmty\r";
	//findCmty(graphs[0],xs[0],cmty);
	cout << "find cmty done\r";
	delete[] ninodes;
	delete[] inodes;
	delete[] nxs;
}


void MultiplexGraph::rwr(uint seed, vector<uint> &cmty, map<uint, double> &prxmts) {

	double ** weights = new double*[nmbofLayer];
	for (int i = 0; i < nmbofLayer; i++) {
		weights[i] = new double[nmbofLayer];
	}
	for (int i = 0; i < nmbofLayer; i++) {
		for (int j = i + 1; j < nmbofLayer; j++) {
			weights[i][j] = ((rand() % 10000) / 10000.0);
			weights[j][i] = weights[i][j];
		}
		weights[i][i] = 1;
	}

	for (uint i = 0; i < nmbofLayer; i++) {
		double sumi = 0;
		for (uint j = 0; j < nmbofLayer; j++) {
			sumi += weights[j][i];
		}
		for (uint j = 0; j < nmbofLayer; j++) {
			weights[j][i] /= sumi;
			cout << weights[j][i] << " ";
		}
		cout << endl;
	}

	xs = new map<uint, double>[nmbofLayer];
	map<uint, double> *nxs = new map<uint, double>[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		xs[i].insert(map<int, double>::value_type(seed, 1.0));
	}
	map<uint, double> es;
	es[seed] = 1;
	int count = 0;
	for (;;count++) {
		for (uint i = 0; i < nmbofLayer; i++) {
			nxs[i].clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				map<uint, double>temp;
				multi(temp, graphs[j], xs[i], 1.0 / nmbofLayer);// weights[j][i]);
				double sum = sumArray(temp);
				addArray(nxs[i], temp, alpha);
			}
			addArray(nxs[i], es, (1 - alpha));
			norm(nxs[i]);
		}
		double diff = mabs(nxs, xs, nmbofLayer);
		if (diff< epsilon||count==100) {
			for (uint i = 0; i < nmbofLayer; i++)
				addArray(prxmts, xs[i], 1.0 / nmbofLayer);
			norm(prxmts);
			break;
		}
		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;
	}
	cout << "begin to find cmty\r";
	//findCmty(graphs[0],xs[0], cmty);
	delete[] nxs;
}




void MultiplexGraph::calConductance(const sparserow &graph, vector<double> &conductances, const vector<int_double_pair_sorted_by_score> tep) {
	int size = tep.size();
	if (size <= 0)
		return;

	//init
	//S = {indexs[0]}
	//hS = V-hS[0];
	int ahs = 0;

	ahs = graph.indexs[graph.nodeSize+1];

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


void MultiplexGraph::findCmty(map<uint, double> &prxmts, uint qlayer, vector<uint> &cmty) {
	if (SweepNormByDegree) {
		for (auto ite : prxmts) {
			prxmts[ite.first] = ite.second / (graphs[qlayer].degree(ite.first) + 1);
		}
	}
	findCmty(prxmts, graphs[qlayer], cmty);
}

void MultiplexGraph::findCmty(map<uint, double> &prxmts, sparserow & graph, vector<uint> &cmty) {
	vector<double> conductances;
	vector<int_double_pair_sorted_by_score> sortedxs;
	sort(sortedxs, prxmts);//based on similarity score, calculate the conductances of each layer.
	if (USETOPK) {
		for (uint i = 0; i < topk&&i < sortedxs.size(); i++)
			cmty.push_back(sortedxs[i].id);
		return;
	}

	calConductance(graph, conductances, sortedxs);
	double bestscore = 1000;
	uint bestindex = 0;

	for (uint i = 1; i < sortedxs.size() && i<max_cmty_size; i++) {
		double score = conductances[i];
		//cout << score << endl;
		if (bestscore > score) {
			bestscore = score;
			bestindex = i;
		}
	}
	if (bestindex <= min_cmty_size)
		bestindex = min_cmty_size;
	for (uint i = 0; i < bestindex && i < sortedxs.size(); i++)
		cmty.push_back(sortedxs[i].id);
}


/*
* focus on qlayer
* sort xs[qlayer] as sortedx
* then based on sortedx, calcualate conductance score for each layer (ws >= 0.1).
* conductance = ws(q,i)C(i)
*/
void MultiplexGraph::findCmty(uint qlayer, double **weights, vector<uint> &cmty) {
	vector<double> * conductances = new vector< double>[nmbofLayer];// conductances on each layer
	vector<int_double_pair_sorted_by_score> sortedxs;
	sort(sortedxs, xs[qlayer]);//based on similarity score, calculate the conductances of each layer.
							   //	cout << sortedxs.size() << endl;

	for (uint i = 0; i < nmbofLayer; i++) {
		if (xs[i].size() == 0 || weights[qlayer][i] <= 0.1)
			continue;
		calConductance(graphs[i], conductances[i], sortedxs);
		//cout <<"conductances[i].size() " <<conductances[i].size() << endl;
	}

	double bestscore = 1000;
	uint bestindex = 0;

	for (uint i = 1; i < sortedxs.size() && i<topk; i++) {
		double score = 0;
		for (uint j = 0; j < nmbofLayer; j++) {
			if (xs[j].size() == 0 || weights[qlayer][j] <= 0.1)
				continue;
			score += weights[qlayer][j] * conductances[j][i];
		}
		if (bestscore > score) {
			bestscore = score;
			bestindex = i;
		}
	}
	if (bestindex <= 5)
		bestindex = 5;
	for (uint i = 0; i < bestindex && i < sortedxs.size(); i++)
		cmty.push_back(sortedxs[i].id);
	delete[] conductances;
}

//sweep with minimum condactance.
void MultiplexGraph::findCmty(map<uint, double> &xs, double * weights, vector<uint> &cmty) {
	sparserow graph;
	graph.weighted = true;
	uint hashroot = 39569;
	set<uint> sedges;
	sedges.insert(0);
	string line;
	vector<mtuple> edges;
	uint maxid = 1;
	for (uint i = 0; i < nmbofLayer; i++) {
		sparserow &sgph = graphs[i];
		if (maxid< sgph.nodeSize) {
			maxid = sgph.nodeSize;
		}
		for (uint node = 0; node<sgph.nodeSize; node++) {
			int rid = node;
			for (uint nbindex = sgph.indexs[node]; nbindex < sgph.indexs[node + 1]; nbindex++) {
				int cid = sgph.nbs[nbindex];
				//double w = sgph.ws[nbindex];
				mtuple ele;
				ele.row = rid;
				ele.col = cid;
				//ele.w = w;
				uint hashc = rid*hashroot + cid;
				if (sedges.find(hashc) == sedges.end()) {
					edges.push_back(ele);
					sedges.insert(hashc);
				}
			}
		}
	}

	sort(edges.begin(), edges.end());

	graph.nodeSize = maxid;
	graph.edgeSize = edges.size();
	graph.nbs = (unsigned int*)malloc(sizeof(unsigned int)*graph.edgeSize);
	graph.indexs = (unsigned int*)malloc(sizeof(unsigned int)*(graph.nodeSize + 2));
	//graph.ws = (double*)malloc(sizeof(double)*graph.edgeSize);

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

	//findCmty(xs, qlayer, cmty);

	free(graph.indexs);
	free(graph.nbs);
	//free(graph.ws);
}


void MultiplexGraph::check() {
	for (uint i = 0; i < nmbofLayer; i++) {
		cout << "checking " << i << "\r";
		graphs[i].check();
	}
}

void MultiplexGraph::findCmty_Density(uint i, vector<uint> &cmty) {
	struct mpair {
		uint id;
		double score;
		bool operator< (mpair i) { return (score>i.score); }
	};

	vector<mpair> tep;

	for (map<uint, double>::iterator it = xs[0].begin(); it != xs[0].end(); ++it) {
		mpair mp;
		mp.id = it->first;
		mp.score = it->second;
		tep.push_back(mp);
	}
	sort(tep.begin(), tep.end());

	int size = tep.size();
	if (size <= 0)
		return;


	int bestsplit = 0;
	float best = 0;

	//init
	//S = {indexs[0]}
	//hS = V-hS[0];
	int ahs = 0;
	sparserow & graph = graphs[i];


	int as = 0;
	int bestindex = 1;
	float bestsocre = as;
	set<uint> before;
	before.insert(tep[0].id);

	for (int i = 1; i < size && i<200; i++) {
		
		for (int k = graph.indexs[tep[i].id]; k < graph.indexs[tep[i].id + 1]; k++) {
			if (before.find(graph.nbs[k])!=before.end()) {
				as++;
			}
		}

		float tempscore = as /(i+1);
		if (tempscore > bestsocre) {
			bestsocre = tempscore;
			bestindex = i;
		}
		before.insert(tep[i].id);
	}

	for (int k = 0; k < bestindex; k++) {
		if (tep[k].score < 0.001)
			break;
		cmty.push_back(tep[k].id);
	}

}


/*
* PPR
*/
void MultiplexGraph::PPR(int layer, map<uint, double> &seedVector, map<uint, double> &prxmts) {
	prxmts.clear();
	
	map<uint, double> xs = seedVector;
	map<uint, double> nxs;
	
	for (;;) {
		multi(nxs, graphs[layer], xs, alpha);
		addArray(nxs, seedVector, (1 - alpha));
		norm(nxs);

		double diff = mabs(nxs, xs);
		if (diff< epsilon) {
			copyArray(prxmts, xs);
			break;
		}
		xs.swap(nxs);
	}
}



sparserow * MultiplexGraph::getGraph(int layerId) {
	if (layerId >= nmbofLayer)
		return NULL;
	return &graphs[layerId];
}


uint MultiplexGraph::getNodeSize() {
	return this->nodesize;
}


//debug
extern int times;
extern double xssum;
void MultiplexGraph::ApproximateMultiWalker(uint seed, uint qlayer, double **weights, double bound, double theta, vector<uint> &sizes) {
	map<uint, double> prxmts;
	for (uint i = 0; i < nmbofLayer; i++) {
		for (uint j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = 1.00;
	}

	map<uint, double> *xs = new map<uint, double>[nmbofLayer];
	map<uint, double> *nxs = new map<uint, double>[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		xs[i].insert(map<int, double>::value_type(seed, 1.0));
	}

	map<uint, double> es;
	es[seed] = 1;

	int T = int(log(bound*(1 - lamda) / nmbofLayer) / log(lamda));
	if (itenumber != -1) {
		T = itenumber;
	}

	bool converge = false;
	int ite = 0;
	for (; ite<T; ite++) {
		for (uint i = 0; i < nmbofLayer; i++) {
			nxs[i].clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[j][i] < 1e-4)
					continue;
				map<uint, double>temp;
				if(theta==1)
					multi(temp, graphs[j], xs[i], weights[j][i], true);
				else
					amulti(temp, graphs[j], xs[i], weights[j][i], theta, seed, visited, true);
				double sum = sumArray(temp);
				addArray(nxs[i], temp, alpha);
			}
			addArray(nxs[i], es, (1 - alpha));
			norm(nxs[i]);
		}
		double diff = mabs(nxs, xs, nmbofLayer);
		if (diff< epsilon) {
			catAndNormArray(prxmts, nxs[qlayer]);
			converge = true;
			break;
		}
		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;


		double decay = lamda;
		if (ite >= 0)
			decay = pow(lamda, ite+1);
		
		for (uint i = 0; i < nmbofLayer; i++) {
			for (uint j = i + 1; j < nmbofLayer; j++) {
				switch (wmodel) {
				case MODEL_PC:weights[i][j] += decay*PospearsonCorrelation(seed, xs[i], xs[j]); break;
				case MODEL_COSINE:weights[i][j] += decay*cosineSimilarity(xs[i], xs[j], seed); break;
				default:
					weights[i][j] += decay*PospearsonCorrelation(seed, xs[i], xs[j]);
				}
				weights[j][i] = weights[i][j];
			}
			weights[i][i] += decay;
		}
		if (!mute)
			cout << "pow(lamda , ite)\t" << ite << "\t" << pow(lamda, ite) << endl;
		for (uint i = 0; i < nmbofLayer; i++) {
			double sumi = 0;
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[i][j] < 1e-4)
					weights[i][j]=0;
				sumi += weights[i][j];
			}
			for (uint j = 0; j < nmbofLayer; j++) {
				weights[i][j] /= sumi;
				if (!mute)
					cout << weights[i][j] << " ";
			}
			if (!mute)
				cout << endl;
		}
		if (!mute)
			cout << "------------------------------" << endl;
	}
	if (!mute)
		cout << "begin to find cmty\r";
	if(fLog)
		(*fLog) << "T\t" << T << "\tdiff\t";
	map<uint, double> & x = xs[qlayer];
	map<uint, double> & nx = nxs[qlayer];

	int p1 = x.size();
	cout << "phase1 " << "\t T"<<T<<"\t ite " <<ite<<"\txsize "<<x.size() << endl;

	if (itenumber == -1) converge = false;

	if (!(converge || phase1Only)) {
		for (ite = T; ite < MAXITE; ite++) {
			nx.clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[j][qlayer] < 1e-4)
					continue;
				map<uint, double>temp;
				multi(temp, graphs[j], x, weights[qlayer][j], true);
				addArray(nx, temp, alpha);
			}
			addArray(nx, es, (1 - alpha));
			norm(nx);
			double diff = mabs(nx, x);
			if (fLog) {
				(*fLog) << diff << " ";
				fLog->flush();
			}
			if (diff < epsilon)
				break;
			x.swap(nx);
		}
		catAndNormArray(prxmts, xs[qlayer]);
	}
	cout << "phase2 " << x.size() << endl;

	int p2 = x.size();
	sizes.push_back(p1);
	sizes.push_back(p2);
	sizes.push_back(ite);

	delete[] nxs;
	delete[] xs;
}




void MultiplexGraph::ApproximateMultiWalker_1st(uint seed, double **weights, double bound, double theta) {
	for (uint i = 0; i < nmbofLayer; i++) {
		for (uint j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = 1.00;
	}

	map<uint, double> *xs = new map<uint, double>[nmbofLayer];
	map<uint, double> *nxs = new map<uint, double>[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		xs[i].insert(map<int, double>::value_type(seed, 1.0));
	}

	map<uint, double> es;
	es[seed] = 1;

	int T = int(log(bound*(1 - lamda) / nmbofLayer) / log(lamda));
	if (itenumber != -1) {
		T = itenumber;
	}

	bool converge = false;
	int ite = 0;
	for (; ite<T; ite++) {
		for (uint i = 0; i < nmbofLayer; i++) {
			nxs[i].clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[j][i] < 1e-4)
					continue;
				map<uint, double>temp;
				if (theta == 1)
					multi(temp, graphs[j], xs[i], weights[j][i], true);
				else
					amulti(temp, graphs[j], xs[i], weights[j][i], theta, seed, visited, true);
				double sum = sumArray(temp);
				addArray(nxs[i], temp, alpha);
			}
			addArray(nxs[i], es, (1 - alpha));
			norm(nxs[i]);
		}
		double diff = mabs(nxs, xs, nmbofLayer);
		if (diff< epsilon)
			break;

		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;


		double decay = lamda;
		if (ite >= 0)
			decay = pow(lamda, ite + 1);

		for (uint i = 0; i < nmbofLayer; i++) {
			for (uint j = i + 1; j < nmbofLayer; j++) {
				switch (wmodel) {
				case MODEL_PC:weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]); break;
				case MODEL_COSINE:weights[i][j] += decay * cosineSimilarity(xs[i], xs[j], seed); break;
				default:
					weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]);
				}
				//cout << ite<<"\t"<<i << "\t" << j << "\t"<<cosineSimilarity(xs[i], xs[j], seed) << endl;
				weights[j][i] = weights[i][j];
			}
			weights[i][i] += decay;
		}
		if (!mute)
			cout << "pow(lamda , ite)\t" << ite << "\t" << pow(lamda, ite) << endl;
		for (uint i = 0; i < nmbofLayer; i++) {
			double sumi = 0;
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[i][j] < 1e-4)
					weights[i][j] = 0;
				sumi += weights[i][j];
			}
			for (uint j = 0; j < nmbofLayer; j++) {
				weights[i][j] /= sumi;
				if (!mute)
					cout << weights[i][j] << " ";
			}
			if (!mute)
				cout << endl;
		}
		if (!mute)
			cout << "------------------------------" << endl;
	}




	if (!mute)
		cout << "begin to find cmty\r";
	if (fLog)
		(*fLog) << "T\t" << T << "\tdiff\t";

	delete[] nxs;
	delete[] xs;
}



void MultiplexGraph::localCmtyDetection(uint querynode, int qlayer, vector<uint> &sizes, double **weights) {
	//init;
	for (int i = 0; i < nmbofLayer; i++) {
		for (int j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = selfsimi;
	}

	ApproximateMultiWalker(querynode, qlayer, weights, epsilon_matrix, theta, sizes);
}

void MultiplexGraph::singlelocalCmtyDetection(uint seed, int qlayer, vector<uint> &cmty) {
	cmty.clear();
	if (qlayer >= nmbofLayer) {
		cout << "multiWalker " << qlayer << " is large than number of layer " << nmbofLayer << endl;
		return;
	}
	map<uint, double> x;

	if (MC) {
		const int nrolls = 4 * nodesize;
		map<int, int> esreults;
		for (int i = 0; i < nrolls; i++) {
			int length = distribution(generator);
			int iid = singleSample(seed, length,qlayer);
			if (iid > -1) {
				esreults[iid] += 1;
			}
		}
		x.clear();
		for (map<int, int>::iterator mite = esreults.begin(); mite != esreults.end(); mite++) {
			x[mite->first] = 1.0*mite->second / nrolls;
		}
	}
	else {
		//initialize
		map<uint, double> vs;
		vs[seed] = 1.0;
		x.insert(map<int, double>::value_type(seed, 1.0));
		map<uint, double> nx;

		int count = 0;
		while (true) {
			count++;
			nx.clear();
			multi(nx, graphs[qlayer], x, alpha);//one step in layer j
			addArray(nx, vs, (1 - alpha));
			norm(nx);
			double diff = mabs(nx, x);
			if (diff < epsilon || count == MAXITE) {
				break;
			}
			x.swap(nx);
		}
		if (!mute)
			cout << "begin to find cmty\r";
	}
	//for (auto tp : x) {
	//	cout << tp.first << "\t" << tp.second << endl;
	//}
	if(NeedFindCommunity)
		findCmty( x, graphs[qlayer], cmty);
}

int MultiplexGraph::getDegree(int qlayer, int seed) {
	if (qlayer >= nmbofLayer) {
		cout << "multiWalker " << qlayer << " is large than number of layer " << nmbofLayer << endl;
		return -1;
	}
	if (seed > graphs[qlayer].nodeSize)
		return 0;
	return graphs[qlayer].indexs[seed + 1] - graphs[qlayer].indexs[seed];
}

int MultiplexGraph::Sample(uint seed, uint step, double *sum) {

	uint mnow = seed;
	int s = 0;
	for (; s < step; s++) {
		double r = abs(((double)rand() / (RAND_MAX)));
		//linear find
		uint i = 0;
		for (; i < nmbofLayer; i++) {
			if (sum[i] > r) {
				break;
			}
		}

		sparserow  &graph = graphs[i];
		if (mnow > graph.nodeSize) {
			return -1;
		}
		uint nsize = graph.indexs[mnow + 1] - graph.indexs[mnow];

		if (nsize == 0) {
			return -1;
		}

		if (graph.weighted) {
			r = abs(((double)rand() / (RAND_MAX)))*graph.volumns[mnow] - 1e-5;
			double sum = 0;
			int nindex = graph.indexs[mnow];
			for (; nindex < graph.indexs[mnow + 1]; nindex++) {
				sum += graph.ws[nindex];
				if (sum > r)
					break;
			}
			mnow = graph.nbs[nindex];
		}
		else {
			int ra = rand() % nsize;
			ra = abs(ra);
			uint nindex = (ra)+graph.indexs[mnow];
			mnow = graph.nbs[nindex];
		}
	}

	return mnow;
}

int MultiplexGraph::singleSample(uint seed, uint step, int qlayer) {

	uint mnow = seed;
	int s = 0;
	for (; s < step; s++) {
		sparserow  &graph = graphs[qlayer];
		if (mnow > graph.nodeSize) {
			return -1;
		}
		uint nsize = graph.indexs[mnow + 1] - graph.indexs[mnow];

		if (nsize == 0) {
			return -1;
		}

		if (graph.weighted) {
			double r = abs(((double)rand() / (RAND_MAX)))*graph.volumns[mnow] - 1e-5;
			double sum = 0;
			int nindex = graph.indexs[mnow];
			for (; nindex < graph.indexs[mnow + 1]; nindex++) {
				sum += graph.ws[nindex];
				if (sum > r)
					break;
			}
			mnow = graph.nbs[nindex];
		}
		else {
			int ra = rand() % nsize;
			ra = abs(ra);
			uint nindex = (ra)+graph.indexs[mnow];
			mnow = graph.nbs[nindex];
		}
	}

	return mnow;
}


void MultiplexGraph::ErrorInA1(uint seed, uint qlayer, int T) {
	double **static_weights = new double*[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		static_weights[i] = new double[nmbofLayer];
		for (uint j = 0; j < nmbofLayer; j++)
			static_weights[i][j] = 0;
		static_weights[i][i] = 1.00;
	}

	map<uint, double> *xs = new map<uint, double>[nmbofLayer];
	map<uint, double> *nxs = new map<uint, double>[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		xs[i].insert(map<int, double>::value_type(seed, 1.0));
	}

	map<uint, double> es;
	es[seed] = 1;

	int ite = 0;
	for (; ite<T; ite++) {
		for (uint i = 0; i < nmbofLayer; i++) {
			nxs[i].clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (static_weights[j][i] < 1e-2)
					continue;
				map<uint, double>temp;
				multi(temp, graphs[j], xs[i], static_weights[j][i], true);
				double sum = sumArray(temp);
				addArray(nxs[i], temp, alpha);
			}
			addArray(nxs[i], es, (1 - alpha));
			norm(nxs[i]);
		}
		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;

		double decay = lamda;
		if (ite >= 0)
			decay = pow(lamda, ite + 1);

		for (uint i = 0; i < nmbofLayer; i++) {
			for (uint j = i + 1; j < nmbofLayer; j++) {
				switch (wmodel) {
				case MODEL_PC:static_weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]); break;
				case MODEL_COSINE:static_weights[i][j] += decay * cosineSimilarity(xs[i], xs[j], seed); break;
				default:
					static_weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]);
				}
				//cout << "ite "<<ite << "\t" << i << "\t" << j << "\t" << cosineSimilarity(xs[i], xs[j], seed) << endl;
				static_weights[j][i] = static_weights[i][j];
			}
			static_weights[i][i] += decay;
		}

		for (uint i = 0; i < nmbofLayer; i++) {
			double sumi = 0;
			for (uint j = 0; j < nmbofLayer; j++) {
				sumi += static_weights[i][j];
			}
			for (uint j = 0; j < nmbofLayer; j++) {
				static_weights[i][j] /= sumi;
				if (!mute)
					cout << static_weights[i][j] << " ";
			}
			if (!mute)
				cout << endl;
		}
		if (!mute)
			cout << "------------------------------" << endl;
	}


	map<uint, double> &x = xs[qlayer];
	map<uint, double> &nx = nxs[qlayer];
	for (int ite = T; ite < MAXITE; ite++) {
		nx.clear();
		for (uint j = 0; j < nmbofLayer; j++) {
			if (static_weights[j][qlayer] < 1e-3)
				continue;
			map<uint, double>temp;
			multi(temp, graphs[j], x, static_weights[qlayer][j], true);
			addArray(nx, temp, alpha);
		}
		addArray(nx, es, (1 - alpha));
		norm(nx);
		double diff = mabs(nx, x);
		if (diff < epsilon)
			break;
		x.swap(nx);
	}
	map<uint, double> prxmts;
	catAndNormArray(prxmts, xs[qlayer]);

	for (uint i = 0; i < nmbofLayer; i++) {
		xs[i].clear();
		nxs[i].clear();
		xs[i].insert(map<int, double>::value_type(seed, 1.0));
	}

	double **weights = new double*[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		weights[i] = new double[nmbofLayer];
		for (uint j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = 1.00;
	}
	double delta1 = 0;
	double delta_ini = 0;
	DiffNorm1(static_weights[qlayer], weights[qlayer], delta1, delta_ini);
	cout << "0\t" << delta1 << "\t" << delta_ini << endl;

	for (ite=0; ite<T; ite++) {
		for (uint i = 0; i < nmbofLayer; i++) {
			nxs[i].clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[j][i] < 1e-2)
					continue;
				map<uint, double>temp;
				multi(temp, graphs[j], xs[i], weights[j][i], true);
				double sum = sumArray(temp);
				addArray(nxs[i], temp, alpha);
			}
			addArray(nxs[i], es, (1 - alpha));
			norm(nxs[i]);
		}
		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;

		double decay = lamda;
		if (ite >= 0)
			decay = pow(lamda, ite + 1);

		for (uint i = 0; i < nmbofLayer; i++) {
			for (uint j = i + 1; j < nmbofLayer; j++) {
				switch (wmodel) {
				case MODEL_PC:weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]); break;
				case MODEL_COSINE:weights[i][j] += decay * cosineSimilarity(xs[i], xs[j], seed); break;
				default:
					weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]);
				}
				//cout << ite << "\t" << i << "\t" << j << "\t" << cosineSimilarity(xs[i], xs[j], seed) << endl;
				weights[j][i] = weights[i][j];
			}
			weights[i][i] += decay;
		}

		for (uint i = 0; i < nmbofLayer; i++) {
			double sumi = 0;
			for (uint j = 0; j < nmbofLayer; j++) {
				sumi += weights[i][j];
			}
			for (uint j = 0; j < nmbofLayer; j++) {
				weights[i][j] /= sumi;
	
			}
		}

		map<uint, double> tx, tnx;
		copyArray(tx,xs[qlayer]);
		copyArray(tnx, xs[qlayer]);
		for (int ite2 = T; ite2 < MAXITE; ite2++) {
			tnx.clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[j][qlayer] < 1e-3)
					continue;
				map<uint, double>temp;
				multi(temp, graphs[j], tx, weights[qlayer][j], true);
				addArray(tnx, temp, alpha);
			}
			addArray(tnx, es, (1 - alpha));
			norm(tnx);
			double diff = mabs(tnx, tx);
			if (diff < epsilon)
				break;
			tx.swap(tnx);
		}
		map<uint, double> prxmts2;
		catAndNormArray(prxmts2, tx);


		double delta1 = 0;
		double delta_ini = 0;
		DiffNorm1(prxmts, prxmts2, delta1, delta_ini);
		cout << ite << " x \t" << delta1 <<"\t" << delta_ini << endl;
	}

	delete[] nxs;
	delete[] xs;
}


void MultiplexGraph::ErrorInA2(uint seed, uint qlayer, int T, float theta) {
	double **static_weights = new double*[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		static_weights[i] = new double[nmbofLayer];
		for (uint j = 0; j < nmbofLayer; j++)
			static_weights[i][j] = 0;
		static_weights[i][i] = 1.00;
	}
	map<uint, double> prxmts;
	map<uint, double> *xs = new map<uint, double>[nmbofLayer];
	map<uint, double> *nxs = new map<uint, double>[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		xs[i].insert(map<int, double>::value_type(seed, 1.0));
	}

	map<uint, double> es;
	es[seed] = 1;

	bool converge = false;
	int ite = 0;
	for (; ite<T; ite++) {
		for (uint i = 0; i < nmbofLayer; i++) {
			nxs[i].clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (static_weights[j][i] < 1e-3)
					continue;
				map<uint, double>temp;
				multi(temp, graphs[j], xs[i], static_weights[j][i], true);
				double sum = sumArray(temp);
				addArray(nxs[i], temp, alpha);
			}
			addArray(nxs[i], es, (1 - alpha));
			norm(nxs[i]);
		}
		double diff = mabs(nxs, xs, nmbofLayer);
		if (diff< epsilon) {
			catAndNormArray(prxmts, nxs[qlayer]);
			converge = true;
			break;
		}
		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;


		double decay = lamda;
		if (ite >= 0)
			decay = pow(lamda, ite + 1);

		for (uint i = 0; i < nmbofLayer; i++) {
			for (uint j = i + 1; j < nmbofLayer; j++) {
				switch (wmodel) {
				case MODEL_PC:static_weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]); break;
				case MODEL_COSINE:static_weights[i][j] += decay * cosineSimilarity(xs[i], xs[j], seed); break;
				default:
					static_weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]);
				}
				static_weights[j][i] = static_weights[i][j];
			}
			static_weights[i][i] += decay;
		}
		for (uint i = 0; i < nmbofLayer; i++) {
			double sumi = 0;
			for (uint j = 0; j < nmbofLayer; j++) {
				sumi += static_weights[i][j];
			}
			for (uint j = 0; j < nmbofLayer; j++) {
				static_weights[i][j] /= sumi;
			}
		}
	}
	map<uint, double> & x = xs[qlayer];
	map<uint, double> & nx = nxs[qlayer];

	for (int ite = T; ite < MAXITE; ite++) {
		nx.clear();
		for (uint j = 0; j < nmbofLayer; j++) {
			if (static_weights[j][qlayer] < 1e-3)
				continue;
			map<uint, double>temp;
			multi(temp, graphs[j], x, static_weights[qlayer][j], true);
			addArray(nx, temp, alpha);
		}
		addArray(nx, es, (1 - alpha));
		norm(nx);
		double diff = mabs(nx, x);
		if (diff < epsilon)
			break;
		x.swap(nx);
	}
	catAndNormArray(prxmts, xs[qlayer]);



	for (uint i = 0; i < nmbofLayer; i++) {
		xs[i].clear();
		nxs[i].clear();
		xs[i].insert(map<int, double>::value_type(seed, 1.0));
	}

	double **weights = new double*[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		weights[i] = new double[nmbofLayer];
		for (uint j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = 1.00;
	}

	converge = false;
	for (ite = 0; ite<T; ite++) {
		for (uint i = 0; i < nmbofLayer; i++) {
			nxs[i].clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[j][i] < 1e-3)
					continue;
				map<uint, double>temp;
				amulti(temp, graphs[j], xs[i], weights[j][i], theta, seed, visited, true);
				double sum = sumArray(temp);
				addArray(nxs[i], temp, alpha);
			}
			addArray(nxs[i], es, (1 - alpha));
			norm(nxs[i]);
		}
		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;

		double decay = lamda;
		if (ite >= 0)
			decay = pow(lamda, ite + 1);

		for (uint i = 0; i < nmbofLayer; i++) {
			for (uint j = i + 1; j < nmbofLayer; j++) {
				switch (wmodel) {
				case MODEL_PC:weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]); break;
				case MODEL_COSINE:weights[i][j] += decay * cosineSimilarity(xs[i], xs[j], seed); break;
				default:
					weights[i][j] += decay * PospearsonCorrelation(seed, xs[i], xs[j]);
				}
				//cout << ite << "\t" << i << "\t" << j << "\t" << cosineSimilarity(xs[i], xs[j], seed) << endl;
				weights[j][i] = weights[i][j];
			}
			weights[i][i] += decay;
		}

		for (uint i = 0; i < nmbofLayer; i++) {
			double sumi = 0;
			for (uint j = 0; j < nmbofLayer; j++) {
				sumi += weights[i][j];
			}
			for (uint j = 0; j < nmbofLayer; j++) {
				weights[i][j] /= sumi;

			}
		}

	}
	map<uint, double> prxmts2;
	double delta1 = 0;
	double delta_ini = 0;
	DiffNorm1(static_weights[qlayer], weights[qlayer], delta1, delta_ini);
	cout << "\n"<<theta << "\t"<< delta_ini<<"\t"<<delta1 << endl;
	if (!(converge)) {
		for (int ite = T; ite < MAXITE; ite++) {
			nx.clear();
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[j][qlayer] < 1e-3)
					continue;
				map<uint, double>temp;
				multi(temp, graphs[j], x, weights[qlayer][j], true);
				addArray(nx, temp, alpha);
			}
			addArray(nx, es, (1 - alpha));
			norm(nx);
			double diff = mabs(nx, x);
			if (diff < epsilon)
				break;
			x.swap(nx);
		}
		catAndNormArray(prxmts2, xs[qlayer]);

		double diff = mabs(nx, x);
		times++;
		xssum += diff;
	}

	double x_norm_1 = 0,x_norm_i = 0;
	DiffNorm1(prxmts, prxmts2, x_norm_1, x_norm_i);
	cout << theta << "\t" << x_norm_i <<"\t"<< x_norm_1<< endl;
	//map<uint, double> prxmts3;
	//for (int i = 0; i < nodesize; i++)
	//	prxmts3[i] = 1.0 / nodesize;
	//DiffNorm1(prxmts, prxmts3, x_norm_1, x_norm_i);
	//cout << "baseline\t" << x_norm_i << "\t" << x_norm_1 << endl;

	////baseline 2
	//x.clear(); x[seed] = 1.0;
	//for (int ite = T; ite < MAXITE; ite++) {
	//	nx.clear();
	//	map<uint, double>temp;
	//	multi(nx, graphs[qlayer], x, alpha, true);
	//	addArray(nx, es, (1 - alpha));
	//	norm(nx);
	//	double diff = mabs(nx, x);
	//	if (diff < epsilon)
	//		break;
	//	x.swap(nx);
	//}
	//prxmts3.clear();
	//catAndNormArray(prxmts3, x);
	//DiffNorm1(prxmts, prxmts3, x_norm_1, x_norm_i);
	//cout << "baseline\t" << x_norm_i << "\t" << x_norm_1 << endl;

	delete[] nxs;
	delete[] xs;
}


double MultiplexGraph::DiffNorm1(double *static_weights, double *weights, double &norm_1, double &norm_i) {
	int * pointers = new int[nmbofLayer];
	int max_node_size = 0;

	for (int i = 0; i < nmbofLayer; i++) {
		pointers[i] = 0;
		if (graphs[i].nodeSize > max_node_size)
			max_node_size = graphs[i].nodeSize;
	}
	double delta_1 = 0;
	double delta_infinity = 0;
	for (int nid = 0; nid < max_node_size; nid++) {
		//find the next one
		for (;;) {
			bool next = true;
			int nid2 = max_node_size + 1;
			for (int lid = 0; lid < nmbofLayer; lid++) {
				if (nid > graphs[lid].nodeSize || pointers[lid] >= graphs[lid].indexs[nid+1])
					continue;
				next = false;
				if (graphs[lid].nbs[pointers[lid]] < nid2)
					nid2 = graphs[lid].nbs[pointers[lid]];
			}
			if (next)
				break;
			double ele1 = 0;
			double ele2 = 0;
			for (int lid = 0; lid < nmbofLayer; lid++) {
				if (nid > graphs[lid].nodeSize || pointers[lid] >= graphs[lid].indexs[nid+1])
					continue;
				if (graphs[lid].nbs[pointers[lid]] == nid2) {

					double weight = 0;
					if (graphs[lid].weighted)
						weight = graphs[lid].ws[pointers[lid]] / graphs[lid].volumns[nid];
					else
						weight = 1.0 / (graphs[lid].indexs[nid + 1] - graphs[lid].indexs[nid]);
					ele1 += weight * static_weights[lid];
					ele2 += weight * weights[lid];
					pointers[lid] += 1;
				}
			}
			double diff = ele1 - ele2;
			if (diff < 0)
				diff = 0 - diff;
			delta_1 += diff;
			if (diff > delta_infinity) {
				delta_infinity = diff;
			//	cout << nid << "\t" << nid2 << "\t"<<delta1 << endl;
			}
		}
	}
	norm_1 = delta_1;
	norm_i = delta_infinity;
	return delta_1;
}


double MultiplexGraph::DiffNorm1(map<uint,double> &m1, map<uint, double> &m2, double &norm_1, double &norm_i) {

	double delta_1 = 0;
	double delta_infinity = 0;
	map<uint, double>::iterator mite1 = m1.begin();
	map<uint, double>::iterator mite2 = m2.begin();
	while (mite1 != m1.end() && mite2 != m2.end()) {
		double diff = 0;
		if (mite1->first == mite2->first) {
			diff = (mite1->second) - (mite2->second);
			if (diff < 0)
				diff = 0 - diff;
			mite1++;
			mite2++;
		}
		else if (mite1->first > mite2->first) {
			diff = (mite2->second);
			mite2++;
		}
		else {
			diff = mite1->second;
			mite1++;
		}
		if (diff > delta_infinity)
			delta_infinity = diff;
		delta_1 += diff;
	}
	while (mite1 != m1.end()) {
		double diff = (mite1->second);
		if (diff > delta_infinity)
			delta_infinity = diff;
		delta_1 += diff;
		mite1++;
	}
	while (mite2 != m2.end()) {
		double diff = (mite2->second);
		if (diff > delta_infinity)
			delta_infinity = diff;
		delta_1 += diff;
		mite2++;
	}


	norm_1 = delta_1;
	norm_i = delta_infinity;
	return delta_1;
}

MultiplexGraph::MultiplexGraph(string dataDir, int session, int layer_begin, int layer_end, bool binary) {
	distribution = geometric_distribution<int>(1 - alpha);
	nodesize = 0;
	nmbofLayer = layer_end - layer_begin + 1;
	graphs = (sparserow*)malloc(sizeof(sparserow)*nmbofLayer);
	string path = dataDir + "sessions/session" + to_string(session) + ".bin";
	ifstream fse(path, ios::in | ios::binary);
	if (!fse) {
		cout << "not found " << path << endl;
		exit(1);
	}
	fse.read((char*)graphs, sizeof(sparserow));
	graphs[0].indexs = new unsigned int[graphs[0].nodeSize + 2];
	graphs[0].nbs = new unsigned int[graphs[0].edgeSize];

	fse.read((char*)graphs[0].indexs, sizeof(unsigned int)*(graphs[0].nodeSize + 2));
	fse.read((char*)graphs[0].nbs, sizeof(unsigned int)*(graphs[0].edgeSize));
	if (graphs[0].weighted) {
		graphs[0].ws = new double[graphs[0].edgeSize];
		graphs[0].volumns = new double[graphs[0].nodeSize + 1];

		fse.read((char*)graphs[0].ws, sizeof(double)*(graphs[0].edgeSize));
		fse.read((char*)graphs[0].volumns, sizeof(double)*(graphs[0].nodeSize + 1));
	}
	fse.close();
	nodesize = graphs[0].nodeSize;

	for (int gid = layer_begin; gid < layer_end; gid++) {
		cout << "read file " << gid << "\r";
		string path = dataDir + "layers/layer" + to_string(gid) + ".bin";

		int lid = gid - layer_begin + 1;
		ifstream fin(path, ios::binary);
		if (!fin) {
			cout << "not found " << path << endl;
			continue;
		}
		fin.read((char*)&graphs[lid], sizeof(sparserow));
		graphs[lid].indexs = new unsigned int[graphs[lid].nodeSize + 2];
		graphs[lid].nbs = new unsigned int[graphs[lid].edgeSize];
		fin.read((char*)graphs[lid].indexs, sizeof(unsigned int)*(graphs[lid].nodeSize + 2));
		fin.read((char*)graphs[lid].nbs, sizeof(unsigned int)*(graphs[lid].edgeSize));
		if (graphs[lid].weighted) {
			graphs[lid].ws = new double[graphs[lid].edgeSize];
			graphs[lid].volumns = new double[graphs[lid].nodeSize + 1];
			fin.read((char*)graphs[lid].ws, sizeof(double)*(graphs[lid].edgeSize));
			fin.read((char*)graphs[lid].volumns, sizeof(double)*(graphs[lid].nodeSize + 1));
		}
		fin.close();
		if (nodesize < graphs[lid].nodeSize)
			nodesize = graphs[lid].nodeSize;
	}

	visited = new bool[nodesize + 1];
	cout << "read file done" << endl;
}
