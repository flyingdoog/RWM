// MultiPPR.cpp : Defines the entry point for the console application.
// Coded By: Dongsheng Luo
// 12/24/2019
// dul262@psu.edu
// College of Information Science and Technology
// The Pennsylvania State University
// University Park PA 16802
// U.S.
#include "Graph.h"
#include"BC.h"
#include"preprocess.h"
#include"Similarities.h"
#include"MultiGraph.h"
#include <ctime>
#define usingCenterNode false
#define USTOPKPRE false
int times;
double xssum;
int step = 1;
using namespace std;
double alpha = 0.9;
int itenumber = -1;
double epsilon = 1e-3;
int min_cmty_size = 10;
int max_cmty_size = 500;

uint topk = 200;
uint qlayer = 0;
ofstream flog("./RWM_log.txt", ios::app);
int tlayer = 10;
double beta = 1;
double lamda = 0.4;
double selfsimi = 1;
bool countTime = false;
double epsilon_matrix = 0.01;
double theta = 0.7;
bool NeedFindCommunity = true, MC = true, phase1Only = false, singleVersion = false, SweepNormByDegree = true ;
bool USETOPK = false;
vector<uint> qlayers;
Wmodel wmodel = MODEL_COSINE;
#define DEBUGG true
ofstream *fLog;


void MultiplexNetworkTIME(string dataset) {

	string path = "C:\\Users\\dul262\\Desktop\\MultiLayerNet\\" + dataset + ".txt";

	MultiplexGraph mgraph = MultiplexGraph(path);
	mgraph.setMute(true);
	if (mgraph.getNmbofLayer()>tlayer)
		mgraph.setNmbofLayer(tlayer);
	int nmbofLayer = mgraph.getNmbofLayer();
	double ** weights = new double*[mgraph.getNmbofLayer()];
	for (int i = 0; i < mgraph.getNmbofLayer(); i++) {
		weights[i] = new double[mgraph.getNmbofLayer()];
		for (int j = 0; j < mgraph.getNmbofLayer(); j++)
			weights[i][j] = 0;
		weights[i][i] = 1;
	}

	mgraph.check();
	int nodesize = mgraph.getNodeSize();
	cout << "begin" << endl;
	clock_t    start, end;
	start = std::clock();
	uint s1 = 0;
	uint s2 = 0;
	uint ite = 0;
	int count = 0;
	for (uint query=1;query<nodesize;query += (nodesize/100)){
		vector<uint> sizes;
		mgraph.localCmtyDetection(query, qlayer, sizes, weights);
		count++;
		s1 += sizes[0];
		s2 += sizes[1];
		ite += sizes[2];
	}
	end = std::clock();
	cout << (s1*1.0 / count) << "\t" << (s2*1.0 / count) <<"\t"<< (ite*1.0 / count)  << "\t"<<((end-start)*1.0/ count)<<endl;
	for (uint i = 0; i < mgraph.getNmbofLayer(); i++)
		delete[] weights[i];
	delete[] weights;
}


int main(int argc, char ** argv) {

	srand(time(NULL));
	fLog = NULL;
	qlayer = 0;
	alpha = 0.9;
	lamda = 0.7;
	theta = 0.9;
	epsilon = 1e-4;
	topk = 500;
	countTime = false; 
	min_cmty_size = 30;
	max_cmty_size = 2000;
	phase1Only = true;
	//itenumber = 100;

	string dataset = "Syn\\differentNodesize\\data\\10k";
	MultiplexNetworkTIME(dataset);
}
