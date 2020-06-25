#include"MultiGraph.h"
#include"Similarities.h"
using namespace std;
#define MC false
bool mute = false;
extern double alpha;
extern double epsilon;
extern double lamda;
extern Wmodel wmodel;
extern double beta;
extern double selfsimi;
extern double epsilon_matrix;
extern double theta;
extern uint MAXITE1st;
extern bool SweepNormByDegree, USETOPK;
extern uint min_cmty_size, max_cmty_size;
extern int ite1st;
#define MAXITE 30
#include <assert.h>
/*
Initiate a zero multigraph
*/
MultiGraph::MultiGraph() {
	nmbofLayer = 0;
	graphs = NULL;
	interGraphs = NULL;
	interFlags = NULL;
	visited = NULL;
}



/*
Input: filepath1: each layer, same format as multiplex network.
filepath2: interGraph Network. each network is divided by "=layer1,layer2="//
*/
MultiGraph::MultiGraph(string &path, string&interpath) {
	distribution = geometric_distribution<int>(1 - alpha);
	ifstream fin(path.c_str());
	assert (!fin.fail());    

	if(!fin){
		cout<<"network file "<<path<<" not found"<<endl;
		exit(1);
	}
	string temp;
	getline(fin, temp);
	uint index_head = 0;
	nmbofLayer = getUint(temp, index_head);
	getline(fin, temp);
	index_head = 0;
	bool directed = getUint(temp, index_head)==1? true: false;

	graphs = new sparserow[nmbofLayer];// (sparserow*)malloc(sizeof()*);
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
			cout << "begin to read #############" << i << "\r";
		sparserow graph;
		graph.directed = directed;
		graphs[i] = graph;
		createGraph(fin, graphs[i], '-');
		if (DEBUG_CREATE_GRAPH){
			graphs[i].check();
			cout << "read layer \t #############" << i << "done\r";
		}
	}
	fin.close();

	ifstream finter(interpath.c_str());
	assert (!finter.fail());    
	if(!finter){
		cout<<"cross network file "<<interpath<<" not found"<<endl;
		exit(1);
	}
	interFlags = new bool*[nmbofLayer];
	interGraphs = new sparserow*[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		interGraphs[i] = new sparserow[nmbofLayer];
		interFlags[i] = new bool[nmbofLayer];
		memset(interFlags[i], false, sizeof(bool)*nmbofLayer);
	}


	if (DEBUG_CREATE_GRAPH)
		cout << "begin to read " << temp << endl;
	 while (!finter.eof()) {
		getline(finter, temp);
		if (temp.length() <= 1)
			continue;
		uint layer1 = 0;
		uint layer2 = 0;
		uint index_head = 0;
		if (temp[0] == '=') {
			layer1 = getUint(temp, index_head);
			layer2 = getUint(temp, index_head);
		}
		sparserow graph;
		if (layer1 < nmbofLayer && layer2 < nmbofLayer && !interFlags[layer1][layer2] && !interFlags[layer1][layer2]) {
			interGraphs[layer1][layer2] = graph;
			interGraphs[layer1][layer2].directed = true;
			createGraph(finter, interGraphs[layer1][layer2], '-');//
			if(DEBUG_CREATE_GRAPH){
				cout << "interGraphs[" << layer1 << "][" << layer2 << "].check()" << endl;
				interGraphs[layer1][layer2].check();
			}
			createGraph(interGraphs[layer1][layer2], interGraphs[layer2][layer1], true);
			if(DEBUG_CREATE_GRAPH){
				cout << "interGraphs[" << layer2 << "][" << layer1 << "].check()" << endl;
				interGraphs[layer2][layer1].check();
			}

			interFlags[layer1][layer2] = true;
			interFlags[layer2][layer1] = true;
		}
		else {
			cout << "skip\t" << temp << endl;
			temp = "";
			while (!finter.eof()) {
				getline(finter, temp);
				if (temp.size() >= 1 && temp[0] == '-')
					break;
			}
		}
		if (DEBUG_CREATE_GRAPH)
			cout << "read layer \t #############" << i << "done\r";
	}
	finter.close();

	int begin, end;
	int nodesize = graphs[0].nodeSize;
	for (int i = 0; i < nmbofLayer; i++) {
		if (nodesize < graphs[i].nodeSize)
			nodesize = graphs[i].nodeSize;
	}
	visited = new bool[nodesize + 1];
}



MultiGraph::~MultiGraph() {
	if (graphs != NULL) {
		for (uint i = 0; i < nmbofLayer; i++) {
			if (graphs[i].indexs != NULL) {
				delete[] graphs[i].indexs;
				graphs[i].indexs = NULL;
			}
			if(graphs[i].nbs != NULL) {
				delete[] graphs[i].nbs;
				graphs[i].nbs = NULL;
			}
			if (graphs[i].ws != NULL) {
				delete[] graphs[i].ws;
				graphs[i].ws = NULL;
			}
			if (graphs[i].volumns != NULL) {
				delete[] graphs[i].volumns;
				graphs[i].volumns = NULL;
			}
		}
		delete[] graphs;
		graphs = NULL;
	}

	if (interFlags != NULL) {
		for (uint i = 0; i < nmbofLayer; i++) {
			for (uint j = 0; j < nmbofLayer; j++) {
				if (!interFlags[i][j])
					continue;
				if (interGraphs[i][j].indexs != NULL) {
					delete[] interGraphs[i][j].indexs;
					interGraphs[i][j].indexs = NULL;
				}
				if (interGraphs[i][j].nbs != NULL) {
					delete[] interGraphs[i][j].nbs;
					interGraphs[i][j].nbs = NULL;
				}
				if (interGraphs[i][j].ws != NULL) {
					delete[] interGraphs[i][j].ws;
					interGraphs[i][j].ws = NULL;
				}
				if (interGraphs[i][j].volumns != NULL) {
					delete[] interGraphs[i][j].volumns;
					interGraphs[i][j].volumns = NULL;
				}
			}
			delete[] interFlags[i];
			delete[] interGraphs[i];
		}
	}
}


/*
* Shrink the sizeof Multi-Networks
*/
void MultiGraph::setNmbofLayer(uint m) {
	if (m < nmbofLayer)
		nmbofLayer = m;
	else
		cout << "setNmbofLayer ERROR!\t" << m << " is layer than original NumberofLayer" << endl;
}


uint MultiGraph::getNmbofLayer() {
	return nmbofLayer;
}


uint MultiGraph::getNodeSize(uint layer) {
	if (layer >= nmbofLayer) {
		cout << "getNodeSize ERROR!\t" << layer << " is layer than NumberofLayer " <<nmbofLayer<< endl;
		return 0;
	}
	return graphs[layer].nodeSize;
}

sparserow * MultiGraph::getLayer(uint layer) {
	if (layer >= nmbofLayer) {
		cout << "getLayer ERROR!\t" << layer << " is layer than NumberofLayer " << nmbofLayer << endl;
		return NULL;
	}
	return &graphs[layer];
}


void MultiGraph::ApproximateMultiWalker(uint seed, uint qlayer, double **weights, double bound, double theta, map<uint,double> &prxmts) {
	prxmts.clear();
	for (uint i = 0; i < nmbofLayer; i++) {
		for (uint j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = 1.00;
	}

	map<uint, double> *xs = new map<uint, double>[nmbofLayer];
	map<uint, double> *nxs = new map<uint, double>[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		for (uint j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = selfsimi;
	}

	//initialize
	xs = new map<uint, double>[nmbofLayer];
	xs[qlayer].insert(map<int, double>::value_type(seed, 1.0));
	for (uint i = 0; i < nmbofLayer; i++) {
		if (interFlags[qlayer][i]) {
			multi(xs[i], interGraphs[qlayer][i], xs[qlayer], 1, false);
			norm(xs[i]);//delete the elements with score less than MINSCORE;
		}
	}

	map<uint, double> *vs = new map<uint, double>[nmbofLayer];
	for (uint i = 0; i < nmbofLayer; i++) {
		copyArray(vs[i], xs[i]);
	}


	int T = int(log(bound*(1 - lamda) / nmbofLayer) / log(lamda));
	if (T> MAXITE1st) {
		T = MAXITE1st;
	}

	if(ite1st>=0)
		T = ite1st;

	bool converge = false;
	int ite = 0;
	for (; ite<T; ite++) {
		for (uint i = 0; i < nmbofLayer; i++) {
			if (xs[i].size() == 0)
				continue;
			nxs[i].clear();
			map<uint, double> part1;
			double left = extractPart(part1, graphs[i], xs[i], theta, seed, visited);
			for (uint j = 0; j < nmbofLayer; j++) {
				if (weights[j][i] < 1e-2) {
					continue;
				}
				map<uint, double>temp1, temp2;
				if (i != j) {
					multi(temp1, interGraphs[i][j], part1, 1, false);//from layer i to layer j
					multi(temp2, graphs[j], temp1, beta, true);//one step RWR in layer j
					addArray(temp2, temp1, 1 - beta);//one step RWR in layer j
					multi(temp1, interGraphs[j][i], temp2, weights[i][j], false);//go back from j to i
					addArray(temp1, vs[i], weights[i][j] * (1 - left));
				}
				else {
					multi(temp1, graphs[j], part1, weights[i][i], true);
					addArray(temp1, vs[i], weights[i][j] * (1 - left));
				}
				addArray(nxs[i], temp1, alpha);
			}
			addArray(nxs[i], vs[i], (1 - alpha));
			norm(nxs[i]);
		}

		double diff = mabs(nxs, xs, nmbofLayer);
		if (diff< epsilon) {
			for (uint i = 0; i < nmbofLayer; i++)
				addArray(prxmts, xs[i], 1.0 / nmbofLayer);
			norm(prxmts);
			converge = true;
			break;
		}
		map<uint, double> *temp = xs;
		xs = nxs;
		nxs = temp;

		//double decay = 1;
		//if (ite >= 1)
		//	decay = pow(lamda, ite);

		double decay = lamda;
		if (ite >= 0)
			decay = pow(lamda, lamda+1);

		///////////////
		for (uint i = 0; i < nmbofLayer; i++) {
			if (xs[i].size() == 0)
				continue;
			for (uint j = 0; j < nmbofLayer; j++) {
				if (xs[j].size() == 0 || (!interFlags[i][j]))
					continue;

				map<uint, double> temp, xsj2i;
				copyArray(temp, xs[j]);
				addArray(temp, vs[j], (alpha - 1));
				multi(xsj2i, interGraphs[j][i], temp, 1, false);//go back from layer j to layer i

				switch (wmodel) {
				case MODEL_PC:weights[i][j] += decay*PospearsonCorrelation(seed, xs[i], xsj2i); break;
				case MODEL_COSINE:weights[i][j] += decay*cosineSimilarity(xs[i], xsj2i, vs[i]); break;
				case MODEL_EC:weights[i][j] += decay*edgeCorrectness(graphs[i], xs[i], xsj2i, seed, theta); break;
				case MODEL_ICS: weights[i][j] += decay*ICS(graphs[i], xs[i], xsj2i, seed, theta); break;
				default:
					weights[i][j] += decay*PospearsonCorrelation(seed, xs[i], xsj2i);
				}
			}
			weights[i][i] += decay*selfsimi;
		}

		if (!mute)
			cout << "pow(lamda , ite)\t" << ite << "\t" << pow(lamda, ite) << endl;
		for (uint i = 0; i < nmbofLayer; i++) {
			double sumi = 0;
			for (uint j = 0; j < nmbofLayer; j++) {
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


	map<uint, double> & x = xs[qlayer];
	map<uint, double> & nx = nxs[qlayer];



	if (weights[qlayer][qlayer] < 0.9)
		used++;

	//debug
	if (!converge) {

		if (!MC) {
			for (int ite = T; ite < MAXITE; ite++) {
				if(!mute)
					cout << "ite " << ite << "\r"; 
				nx.clear();
				for (uint j = 0; j < nmbofLayer; j++) {
					if (weights[qlayer][j] < 1e-3)
						continue;

					map<uint, double>temp1, temp2;
					if (qlayer != j) {
						multi(temp1, interGraphs[qlayer][j], x, 1, false);//from layer i to layer j
						multi(temp2, graphs[j], temp1, beta, true);//one step RWR in layer j
						addArray(temp2, temp1, 1 - beta);//one step RWR in layer j
						multi(temp1, interGraphs[j][qlayer], temp2, weights[qlayer][j], false);//go back from j to i
					}
					else
						multi(temp1, graphs[qlayer], xs[qlayer], weights[qlayer][qlayer], true);//one step in layer j
					addArray(nx, temp1, alpha);
				}
				addArray(nx, vs[qlayer], (1 - alpha));
				norm(nx);
				double diff = mabs(nx, x);
				if (diff < epsilon)
					break;
				x.swap(nx);
			}
			catAndNormArray(prxmts, x);
		}
		else
		{
			const int nrolls = 10000;// *graphs[qlayer].nodeSize;
			map<int, int> esreults;
			double *sum = weights[qlayer];
			for (uint j = 1; j < nmbofLayer; j++) {
				sum[j] += sum[j - 1];
			}
			sum[nmbofLayer - 1] = 1.1;
			for (int i = 0; i < nrolls;i++ ) {
				int length = distribution(generator);
				int iid = Sample(qlayer, seed, length, sum);
				if (iid > -1) {
					esreults[iid] += 1;
				}
				//if (i % 1000000 == 0) {
				//	x.clear();
				//	for (map<int, int>::iterator mite = esreults.begin(); mite != esreults.end(); mite++) {
				//		x[mite->first] = 1.0*mite->second / nrolls;
				//	}
				//	catAndNormArray(nx, x);
				//	cout << "i\t" << i << "\t" << mabs(prxmts, nx) << "\r";
				//}
			}
			x.clear();
			for (map<int, int>::iterator mite = esreults.begin(); mite != esreults.end(); mite++) {
				x[mite->first] = 1.0*mite->second / nrolls;
			}
		}

	}

	//catAndNormArray(nx, x);
	//cout <<"query\t"<< seed<<"\t"<<mabs(prxmts, nx) << endl;
	catAndNormArray(prxmts, x);
	delete[] nxs;
	delete[] xs;
}


void MultiGraph::StaticMultiWalker(uint seed, uint qlayer, double **weights, map<uint, double> &prxmts) {
	map<uint, double>  x,v;
	map<uint, double>  nx;
	x[seed] = 1;
	v[seed] = 1;



	for (int ite = 0; ite<MAXITE; ite++) {
		nx.clear();
		for (uint j = 0; j < nmbofLayer; j++) {
			if (weights[qlayer][j] < 1e-3)
				continue;

			map<uint, double>temp1, temp2;
			if (qlayer != j) {
				temp2= x;
				multi(temp1, interGraphs[qlayer][j], x, 1, false);//from layer i to layer j
				multi(temp2, graphs[j], temp1, beta, true);//one step RWR in layer j
				addArray(temp2, temp1, 1 - beta);//one step RWR in layer j
				multi(temp1, interGraphs[j][qlayer], temp2, weights[qlayer][j], false);//go back from j to i
			}
			else
				multi(temp1, graphs[qlayer], x, weights[qlayer][qlayer], true);//one step in layer j
			addArray(nx, temp1, alpha);
		}
		addArray(nx, v, (1 - alpha));
		norm(nx);
		double diff = mabs(nx, x);
		if (diff < epsilon)
			break;
		x.swap(nx);
	}
	catAndNormArray(prxmts, x);
}

sparserow * MultiGraph::getInterGraph(uint l1, uint l2) {
	if (l1 >= nmbofLayer || l2>nmbofLayer) {
		cout << "getLayer ERROR!\t" << l1<<"or "<<l2<< " is layer than NumberofLayer " << nmbofLayer << endl;
		return NULL;
	}
	if(interFlags[l1][l2])
		return &interGraphs[l1][l2];
	return NULL;
}


void MultiGraph::setMute(bool m) {
	mute = m;
}

void MultiGraph::findCmty(map<uint,double> &prxmts, uint qlayer, vector<uint> &cmty) {

	if (SweepNormByDegree) {
		for (auto tp : prxmts) {
			prxmts[tp.first] = tp.second / (graphs[qlayer].degree(tp.first) + 1);
		}
	}

	vector<double> conductances;
	vector<int_double_pair_sorted_by_score> sortedxs;
	sort(sortedxs, prxmts);//based on similarity score, calculate the conductances of each layer.
	if (USETOPK) {
		for (uint i = 0; i < max_cmty_size&&i < sortedxs.size(); i++)
			cmty.push_back(sortedxs[i].id);
		return;
	}
	//cout << "begin to calculate conductance" << endl;
	calConductance(graphs[qlayer], conductances, sortedxs);
	//cout << "conductance done" << endl;

	double bestscore = 1000;
	uint bestindex = 0;

	for (uint i = 1; i < sortedxs.size() && i<max_cmty_size; i++) {
		double score = conductances[i];
		if (bestscore > score) {
			bestscore = score;
			bestindex = i;
		}
	}
	if (bestindex <=min_cmty_size )
		bestindex = min_cmty_size;
	for (uint i = 0; i < bestindex && i < sortedxs.size(); i++)
		cmty.push_back(sortedxs[i].id);
}

void MultiGraph::localCmtyDetection(uint querynode, int qlayer, vector<uint> &cmty, double **weights) {
	cmty.clear();

	//init;
	for (int i = 0; i < nmbofLayer; i++) {
		for (int j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = selfsimi;
	}

	map<uint, double> prxmts;
	ApproximateMultiWalker(querynode, qlayer, weights, epsilon_matrix, theta, prxmts);
	if (!mute)
		cout << "begin to find cmty\r";
	findCmty(prxmts, qlayer, cmty);
}

void MultiGraph::StaticlocalCmtyDetection(uint querynode, int qlayer, vector<uint> &cmty, double **weights) {
	cmty.clear();

	//init;
	for (int i = 0; i < nmbofLayer; i++) {
		for (int j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0.5;
		weights[i][i] = 0.5;
	}

	map<uint, double> prxmts;
	StaticMultiWalker(querynode, qlayer, weights, prxmts);
	findCmty(prxmts, qlayer, cmty);
}


int MultiGraph::getDegree(int qlayer, int seed) {
	if (qlayer >= nmbofLayer) {
		cout << "multiWalker " << qlayer << " is large than number of layer " << nmbofLayer << endl;
		return 0;
	}
	if (seed > graphs[qlayer].nodeSize)
		return 0;
	return graphs[qlayer].indexs[seed + 1] - graphs[qlayer].indexs[seed];
}


int MultiGraph::Sample (int qlayer, uint seed, uint step, double *sum) {

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

		uint nid = mnow;
		if (i != qlayer) {
			sparserow &goMatrix = interGraphs[qlayer][i];
			if (mnow> goMatrix.nodeSize) {
				break;
			}

			uint nsize = goMatrix.indexs[mnow + 1] - goMatrix.indexs[mnow];
			if (nsize == 0) {
				return -1;
			}
			//unweighted
			int ra = rand() % nsize;
			ra = abs(ra);
			uint nindex = (ra)+goMatrix.indexs[mnow];
			nid = goMatrix.nbs[nindex];

			sparserow  &graph = graphs[i];
			nsize = graph.indexs[nid + 1] - graph.indexs[nid];
			if (nsize == 0) {
				return -1;
			}

			//unweighted
			uint onestepId = 0;
			if (!graph.weighted) {
				int ra = rand() % nsize;
				ra = abs(ra);
				uint nindex = (ra)+graph.indexs[nid];
				onestepId = graph.nbs[nindex];
			}
			else {
				if (graph.indexs[nid] == graph.indexs[nid + 1]) {
					onestepId = nid;
				}
				else {
					double r = abs(((double)rand() / (RAND_MAX)))*graph.volumns[nid] - 1e-3;
					double sum = 0;
					uint nindex = 0;
					for (nindex = graph.indexs[nid]; nindex < graph.indexs[nid + 1]; nindex++) {
						sum += graph.ws[nindex];
						if (sum > r)
							break;
					}
					onestepId = graph.nbs[nindex];
				}
			}

			sparserow &backMatrix = interGraphs[i][qlayer];
			if (onestepId> backMatrix.nodeSize) {
				return -1;
			}

			nsize = backMatrix.indexs[onestepId + 1] - backMatrix.indexs[onestepId];
			//cout << "----nsize#" << nsize <<"#\r";
			if (nsize == 0) {
				return -1;
			}
			//unweighted
			ra = rand() % nsize;
			ra = abs(ra);
			nindex = (ra)+backMatrix.indexs[onestepId];
			mnow = backMatrix.nbs[nindex];
		}
		else {
			sparserow  &graph = graphs[qlayer];
			uint nsize = graph.indexs[nid + 1] - graph.indexs[nid];
			if (nsize == 0) {
				return -1;
			}

			if (!graph.weighted) {
				int ra = rand() % nsize;
				ra = abs(ra);
				int nindex = (ra)+graph.indexs[nid];
				mnow = graph.nbs[nindex];
			}
			else {
				if (graph.indexs[nid] == graph.indexs[nid + 1])
					break;
				double r = abs(((double)rand() / (RAND_MAX)))*graph.volumns[nid] - 1e-3;
				double sum = 0;
				uint nindex = 0;
				for (nindex = graph.indexs[nid]; nindex < graph.indexs[nid + 1]; nindex++) {
					sum += graph.ws[nindex];
					if (sum > r)
						break;
				}
				mnow = graph.nbs[nindex];
			}
		}
	}
	return mnow;
}

int MultiGraph::singleSample(uint seed, uint step, int qlayer) {

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

