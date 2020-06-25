#pragma once
#include"tools.h"
#include"Graph.h"
#include<random>
#ifndef _HEAD_MULTIPGTAPH
#define _HEAD_MULTIPGTAPH

class MultiGraph {
private:
	uint nmbofLayer;
	sparserow * graphs;//store each layer
	sparserow ** interGraphs; // store the graphs between layers. if not exist, then interGraphs[i][j]=NULL;
	bool **interFlags;
	bool *visited;
	//void findCmty(uint qlayer, double **weights, vector<uint>& cmty);
	void findCmty(map<uint,double> &, uint qlayer,vector<uint>& cmty);
	default_random_engine generator;
	geometric_distribution<int> distribution;
	int singleSample(uint seed, uint step, int qlayer);
	int Sample(int qlayer, uint seed, uint step, double *sum);
public:
	MultiGraph();
	MultiGraph(string &, string&);
	~MultiGraph();
	int used = 0;
	int MultiGraph::getDegree(int qlayer, int seed);

	void setNmbofLayer(uint);
	uint getNmbofLayer();
	uint getNodeSize(uint );
	
	void singlelocalCmtyDetection(uint querynode, int qlayer, vector<uint> &cmty);
	sparserow * getLayer(uint);
	sparserow * getInterGraph(uint,uint);
	void localCmtyDetection(uint querynode, int qlayer, vector<uint> &cmty, double **);
	void StaticlocalCmtyDetection(uint querynode, int qlayer, vector<uint> &cmty, double **);
	void setMute(bool m);
	void ApproximateMultiWalker(uint seed, uint qlayer, double **weights,  double bound, double theta, map<uint, double> &prxmts);
	void StaticMultiWalker(uint seed, uint qlayer, double **weights, map<uint, double> &prxmts);
};
#endif // !_HEAD_MULTIPGTAPH

