#pragma once
#include"tools.h"
#include"Graph.h"
#include <random>
#ifndef _HEAD_MULTIPLEXGTAPH
#define _HEAD_MULTIPLEXGTAPH
class MultiplexGraph {
private:
	uint nmbofLayer;
	uint nodesize;
	sparserow * graphs;
	bool mute;
	map<uint, double> *xs;
	bool *visited;
	void ApproximateMultiWalker(uint seed, uint qlayer, double **weights,double bound, double theta, vector<uint> &prxmts);

	void findCmty(map<uint, double> &prxmts, uint qlayer, vector<uint> &cmty);

	default_random_engine generator;
	geometric_distribution<int> distribution;
	int Sample(uint seed, uint step, double *sum);
	int singleSample(uint seed, uint step, int qlayer);

	void calConductance(const sparserow &graph, vector<double> &conductances, const vector<int_double_pair_sorted_by_score> tep);
	double DiffNorm1(double *static_weights, double *weights,double&, double&);
	double DiffNorm1(map<uint, double> &m1, map<uint, double> &m2, double &norm_1, double &norm_i);

public:
	MultiplexGraph();
	void setNmbofLayer(uint);
	uint getNodeSize();
	sparserow * getGraph(int);
	MultiplexGraph(string);
	MultiplexGraph(string, int);
	MultiplexGraph(string dataDir, int session, int layer_begin, int layer_end, bool binary);
	void ApproximateMultiWalker_1st(uint seed, double **weights, double bound, double theta);

	uint getNmbofLayer();
	~MultiplexGraph();
	void setMute(bool m);
	void check();
	void findCmty(map<uint, double> &xs, double * weights, vector<uint> &cmty);
	void findCmty_Density(uint i, vector<uint> &cmty);
	void findCmty(map<uint, double> &xs, sparserow &graphs, vector<uint> &cmty);
	void findCmty(uint qlayer, double **weights, vector<uint> &cmty);
	int getDegree(int qlayer, int s);
	void ErrorInA1(uint seed, uint qlayer, int T);
	void ErrorInA2(uint seed, uint qlayer, int T, float theta);

	//void multiWalker(uint seed, vector<uint> *cmty,double **);
	void rwr(uint seed, vector<uint> &cmty, map<uint, double> &prxmts);
	void basicMultiWalker(uint seed, vector<uint> &cmty, map<uint, double> &prxmts);
	void AmultiWalker(uint seed, vector<uint> &cmty, map<uint, double> &prxmts);
	void PPR(int layer, map<uint, double> &seedVector, map<uint, double> &prxmts);
	void localCmtyDetection(uint querynode, int qlayer, vector<uint> &cmty, double **weights);
	void singlelocalCmtyDetection(uint seed, int qlayer, vector<uint> &cmty);

};
#endif // !_HEAD_MULTIPLEXGTAPH
