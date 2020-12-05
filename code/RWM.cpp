//Main entrance of the RWM
//Release Version
//C++ implementation for submission "Local Community Detection in Multiple Networks"
//The demo dataset used is 20news-6ng. For other datasets, please refer to references.

#include "Graph.h"
#include"debug.h"
#include"Similarities.h"
#include"MultiGraph.h"
#include <ctime>
using namespace std;

double alpha = 0.9;
double epsilon = 1e-3;
uint qlayer = 0;
double beta = 1;
double lamda = 0.4;
double selfsimi = 1;
double epsilon_matrix = 0.01;
double theta = 0.7;
bool NeedFindCommunity = true, MC = true, phase1Only = false, SweepNormByDegree = true ;
int min_cmty_size = 5;
int max_cmty_size = 200;
bool USETOPK = false;
Wmodel wmodel = MODEL_COSINE;
uint MAXITE1st=5;
uint step =1;
ofstream flog("./tuning.txt",ios::app);
int ite1st =-1;
string settings = "";
int RWM_mute = 0;

void readCmt(string cmty, vector<set<unsigned int>> &gts) {
	ifstream fin(cmty);
	string line;
	assert (!fin.fail()); 

	if(!fin){
		cout<<"community file "<<cmty<<" not found"<<endl;
		exit(1);
	}
	vector<string> ss;
	while (!fin.eof()) {
		getline(fin, line);
		ss.clear();
		split(ss, line, '\t');
		set<unsigned int> cmt;
		gts.push_back(cmt);
		set<unsigned int> &cmt1 = gts[gts.size()-1];
		for (vector<string>::iterator vite = ss.begin(); vite != ss.end(); vite++) {
			cmt1.insert(atoi(vite->c_str()));
		}
	}
	fin.close();
}

void evaluate(const vector<unsigned int> &res, const set<unsigned int> &gt, double &pre, double &rec, double &f1) {
	assert(res.size() * gt.size() != 0);
	unsigned int com = 0;
	for (vector<unsigned int>::const_iterator vite = res.begin(); vite != res.end(); vite++) {
		unsigned int node = *vite;
		if (gt.find(node) != gt.end()) {
			com += 1;
		}
	}
	if (com == 0) {
		pre = 0;
		rec = 0;
		f1 = 0;
		return;
	}
	pre = com*1.0 / res.size();
	rec = com*1.0 / gt.size();
	f1 = pre*rec * 2 / (pre + rec);
}


void generalCase(string path,string path2,string cmtyPath) {

	cout<<"begin to read network"<<endl;
	MultiGraph mgraph = MultiGraph(path,path2);
	mgraph.setMute(!DEBUGG);
	cout<<"Read network done"<<endl;
	
	uint qnodesize = mgraph.getNodeSize(qlayer);// D

	int *labels = new int[qnodesize + 1];
	for (uint i = 0; i < qnodesize + 1; i++)
		labels[i] = -1;
	vector<set<unsigned int>> gts;
	readCmt(cmtyPath, gts);
	uint labelsize = 0;
	for (uint i = 0; i < gts.size(); i++) {
		set<unsigned int>::iterator site = gts[i].begin();
		for (; site != gts[i].end(); site++) {
			if (labels[*site] == -1)
				labels[*site] = i;
			else
				labels[*site] = -3;// means have multiple labels, just ignore it.
			labelsize++;
		}
	}

	int nmbofLayer = mgraph.getNmbofLayer();

	double apre = 0, arecall = 0, af1 = 0, sumpre = 0, sumrecall = 0, sumf1 = 0;
	double apre2 = 0, arecall2 = 0, af12 = 0, sumpre2 = 0, sumrecall2 = 0, sumf12 = 0;
	uint count = 1;
	double ** weights = new double*[nmbofLayer];
	for (int i = 0; i < nmbofLayer; i++) {
		weights[i] = new double[nmbofLayer];
		for (int j = 0; j < nmbofLayer; j++)
			weights[i][j] = 0;
		weights[i][i] = selfsimi;
	}

	int step_t = step;
	clock_t    start, end;
	start = std::clock();
	for (unsigned int i = 1; i <= qnodesize; i+= step_t) {
		step_t =1;

		if (labels[i] <0 )
			continue;
	
		if (mgraph.getDegree(qlayer, i) == 0)
			continue;
		if (!RWM_mute)
			cout << i << "\r";
		step_t = step;
		unsigned int querynode = i;
		uint label = labels[querynode];
		double f1 = 0, pre = 0, rec = 0;
		vector<uint> res1, res2;
		set<uint> reset;


		vector<uint> cmty1;
		
		mgraph.localCmtyDetection(querynode, qlayer, cmty1, weights);
		for (vector<uint>::iterator vite = cmty1.begin(); vite != cmty1.end(); vite++) {

			if (*vite <= qnodesize && labels[*vite] != -1 && reset.find(*vite) == reset.end()) {
				res1.push_back(*vite);
				reset.insert(*vite);
			}
		}
		evaluate(res1, gts[label], pre, rec, f1);

		sumf1 += f1;
		sumpre += pre;
		sumrecall += rec;

		af1 = sumf1 / (count);
		arecall = sumrecall / (count);
		apre = sumpre / (count);

		if (!RWM_mute)
			cout << "query "<< i << " current average pre: " << apre << " rec: " << arecall << " f1: " << af1 << endl;
		// cout.flush();
		count+=1;
	}

	cout << "RWM:  Macro results" << " pre: " << apre << " rec: " << arecall << " f1: " << af1 << endl;
	end = std::clock();
	cout << (std::clock() - start) << "\t" << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) / count << std::endl;
	flog << settings <<endl;
	flog << "RWM:  Macro results" << " pre: " << apre << " rec: " << arecall << " f1: " << af1 << "\t"<<(std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) / count << std::endl;
	flog << endl;

	delete [] labels;
	for (uint i = 0; i < mgraph.getNmbofLayer(); i++)
		delete[] weights[i];
	delete[] weights;
}


void simpleCase(string path, string cmtyPath) {
	cout<<"not implemented yet"<<endl;

}

void print_help(){
			cout << "\nEXAMPLE:\n" <<
			"./RWM -d 6ng "<< endl;
		cout << "Parameters:\n" <<
			"-d dataset name, 6ng, 9ng,...\n" <<
			"-n the target network\n" <<
			"-a alpha: default: 0.9\n" <<
			"-l lamda: default: 0.7 \n" <<
			"-t theta: default: 0.9 \n" <<
			"-e epsilon: default: 0.01 \n" <<endl;
}

int	main_entrance(int argc, char** argv){
	srand(time(NULL));
	qlayer = 0;
	alpha = 0.9;
	lamda = 0.7;
	theta = 0.9;
	epsilon = 1e-2;

	vector<string> simpleDatasets; simpleDatasets.push_back("rm");simpleDatasets.push_back("brainNet");
	vector<string> generalDatasets;
	generalDatasets.push_back("6ng");
	generalDatasets.push_back("9ng");
	generalDatasets.push_back("dblp");
	generalDatasets.push_back("citeseer");


	string dataset = "6ng";

	string dirPath = "./data/"+dataset;
	string network_path = dirPath+"/networks.txt";
	string cross_network_path2 = dirPath+"/inter.txt";
	string cmtyPath = dirPath+"/cmty"+to_string(qlayer)+".txt";


	// 2. Parsing the input parameters.
	int nIdx = 1; // index for parsing the input parameters
	while (nIdx < argc){
		string sToken1 = argv[nIdx];
		string sToken2 = argv[nIdx + 1];
		nIdx = nIdx + 2;
		if (sToken1.compare("-d") == 0) { // Input file name for nodes of the networks
			dataset = sToken2;
			// flog <<"data "<<dataset <<" ";
			settings += "data "+dataset+" ";
			dirPath= "./data/"+sToken2;
			network_path = dirPath+"/networks.txt";
			cross_network_path2 = dirPath+"/inter.txt";
			cmtyPath = dirPath+"/cmty"+to_string(qlayer)+".txt";
		}

		else if (sToken1.compare("-n") == 0) { // the target network
			qlayer = atoi(sToken2.c_str());
			// flog <<"qlayer "<<qlayer <<" ";
			settings +="qlayer " + std::to_string(qlayer) +" ";
			cmtyPath = dirPath+"/cmty"+to_string(qlayer)+".txt";
		}

		else if (sToken1.compare("-a") == 0) { // parameter: aplha
			alpha = atof(sToken2.c_str());
			// flog <<"alpha "<<alpha <<" ";
			settings +="alpha " + std::to_string(alpha) +" ";
		}

		else if (sToken1.compare("-l") == 0) {// parameter: beta
			lamda = atof(sToken2.c_str());
			// flog <<"lamda "<<lamda <<" ";
			settings +="lamda " + std::to_string(lamda) +" ";
		}

		else if (sToken1.compare("-t") == 0) { // parameter: theta
			theta = atof(sToken2.c_str());
			// flog <<"theta "<<theta <<" ";
			settings +="theta " + std::to_string(theta) +" ";
		}

		else if (sToken1.compare("-e") == 0) {// parameter: gamma_W
			epsilon_matrix = atof(sToken2.c_str());
			// flog <<"epsilon_matrix "<<epsilon_matrix <<" ";
			settings +="epsilon_matrix " + std::to_string(epsilon_matrix) +" ";
		}
		
		else if (sToken1.compare("-s") == 0) {// parameter: gamma_W
			step = atoi(sToken2.c_str());
			// flog <<"step "<<step <<" ";
			settings +="step " + std::to_string(step) +" ";
		}

		else if (sToken1.compare("-m") == 0) {// parameter: gamma_W
			RWM_mute = atoi(sToken2.c_str());
		}

		else if (sToken1.compare("-min") == 0) {// parameter: gamma_W
			min_cmty_size = atoi(sToken2.c_str());
			// flog <<"min_cmty_size "<<min_cmty_size <<" ";
			settings +="min_cmty_size " + std::to_string(min_cmty_size) +" ";
		}
		
		else if (sToken1.compare("-max") == 0) {// parameter: gamma_W
			max_cmty_size = atoi(sToken2.c_str());
			// flog <<"max_cmty_size "<<max_cmty_size <<" ";
			settings +="max_cmty_size " + std::to_string(max_cmty_size) +" ";
		}

		else if (sToken1.compare("-ite1st") == 0) {// parameter: gamma_W
			ite1st = atoi(sToken2.c_str());
			// flog <<"ite1st "<<ite1st <<endl;
			settings +="ite1st " + std::to_string(ite1st) +" ";
		}

		else{
			print_help();
			return 1;
		}
	}
	if(std::find(simpleDatasets.begin(),simpleDatasets.end(),dataset)!=simpleDatasets.end())
		simpleCase(network_path,cmtyPath);
	else if(std::find(generalDatasets.begin(),generalDatasets.end(),dataset)!=generalDatasets.end())
		generalCase(network_path,cross_network_path2,cmtyPath);
	else
		cout<<"dataset error"<<endl;
	flog.close();
	return 0;

}



int main(int argc,char ** argv) {//**************************

	if (argc==1|| string(argv[1]).compare("-h") == 0 || string(argv[1]).compare("-help") == 0) {
		//sample case
		print_help();
	}
	else
		main_entrance(argc, argv);
}


