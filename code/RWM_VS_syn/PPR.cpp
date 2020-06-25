#include "stdafx.h"
/**
* @file pprclus_mex.cc
* Implement a PPR clustering scheme.
*
* mex pprclus_mex.cc CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims
*/


#include <vector>
#include <queue>
#include <utility> // for pair sorting
#include <assert.h>
#include <limits>
#include <algorithm>
#include<iostream>

#include <unordered_set>
#include <unordered_map>

#include"Graph.h"
#ifndef __APPLE__
#define __STDC_UTF_16__ 1
#endif
typedef unsigned int	mwIndex;
using namespace std;

/** A replacement for std::queue<int> using a circular buffer array */
class array_queue {
public:
	std::vector<int> array;
	size_t max_size;
	size_t head, tail;
	size_t cursize;
	array_queue(size_t _max_size)
		: max_size(_max_size), array(_max_size), head(0), tail(0), cursize(0)
	{}

	void empty() {
		head = 0;
		tail = 0;
		cursize = 0;
	}

	size_t size() {
		return cursize;
	}

	void push(int i) {
		assert(size() < max_size);
		array[tail] = i;
		tail++;
		if (tail == max_size) {
			tail = 0;
		}
		cursize++;
	}

	int front() {
		assert(size() > 0);
		return array[head];
	}

	void pop() {
		assert(size() > 0);
		head++;
		if (head == max_size) {
			head = 0;
		}
		cursize--;
	}
};

int sr_degree(const sparserow &s, unsigned int u) {

	return (s.indexs[u + 1] - s.indexs[u]);
}

double sr_volumn(const sparserow &s, unsigned int u) {

	double sum = 0;
//	for (int i = s.indexs[u]; i<s.indexs[u + 1]; i++)
//		sum += s.ws[i];
	//std:cout<<sum<<endl;
	return sum;
}

template <class Queue>
int compute_local_pagerank(const sparserow &s, sparsevec& r, sparsevec& p,
	double alpha, double epsilon, int max_push_count, Queue& q)
{
	for (sparsevec::map_type::iterator it = r.map.begin(), itend = r.map.end();
		it != itend; ++it) {
		if (it->second > epsilon*sr_degree(s, it->first)) {
			q.push(it->first);
		}
	}

	int push_count = 0;
	while (q.size()>0 && push_count < max_push_count) {
		push_count += 1;
		mwIndex u = q.front();
		q.pop();
		mwIndex du = sr_degree(s, u);
		double moving_probability = r.map[u] - 0.5*epsilon*(double)du;
		r.map[u] = 0.5*epsilon*(double)du;
		p.map[u] += (1. - alpha)*moving_probability;

		double neighbor_update = alpha*moving_probability / (double)du;

		for (mwIndex nzi = s.indexs[u]; nzi<s.indexs[u + 1]; nzi++) {
			mwIndex x = s.nbs[nzi];
			mwIndex dx = sr_degree(s, x);
			double rxold = r.get(x);
			double rxnew = rxold + neighbor_update;
			r.map[x] = rxnew;
			if (rxnew > epsilon*dx && rxold <= epsilon*dx) {
				q.push(x);
			}
		}
	}

	return (push_count);
}



struct greater2nd {
	template <typename P> bool operator() (const P& p1, const P& p2) {
		return p1.second > p2.second;
	}
};

void cluster_from_sweep(const sparserow & G, sparsevec& p,
	std::vector<unsigned int>& cluster, double *outcond, double* outvolume,
	double *outcut)
{
	typedef std::vector< std::pair<int, double> > vertex_prob_type;
	vertex_prob_type prpairs(p.map.begin(), p.map.end());
	std::sort(prpairs.begin(), prpairs.end(), greater2nd());

	// compute cutsize, volume, and conductance
	std::vector<double> conductance(prpairs.size());
	std::vector<unsigned int> volume(prpairs.size());
	std::vector<unsigned int> cutsize(prpairs.size());

	size_t i = 0;
	tr1ns::unordered_map<int, size_t> rank;
	for (vertex_prob_type::iterator it = prpairs.begin(), itend = prpairs.end();
		it != itend; ++it, ++i) {
		rank[it->first] = i;
	}
	//printf("support=%i\n",prpairs.size());
	unsigned int total_degree = G.indexs[G.nodeSize];
	unsigned int curcutsize = 0;
	unsigned int curvolume = 0;
	i = 0;
	for (vertex_prob_type::iterator it = prpairs.begin(), itend = prpairs.end();
		it != itend; ++it, ++i) {
		unsigned int v = it->first;
		unsigned int deg = G.indexs[v + 1] - G.indexs[v];
		unsigned int change = deg;
		for (unsigned int nzi = G.indexs[v]; nzi<G.indexs[v + 1]; ++nzi) {
			unsigned int nbr = G.nbs[nzi];
			if (rank.count(nbr) > 0) {
				if (rank[nbr] < rank[v]) {
					change -= 2;
				}
			}
		}
		curcutsize += change;
		//if (curvolume + deg > target_vol) {
		//break;
		//}
		curvolume += deg;
		volume[i] = curvolume;
		cutsize[i] = curcutsize;
		if (curvolume == 0 || total_degree - curvolume == 0) {
			conductance[i] = 1;
		}
		else {
			conductance[i] = (double)curcutsize /
				(double)std::min(curvolume, total_degree - curvolume);
		}
		//printf("%5i : cut=%6i vol=%6i prval=%8g cond=%f\n", i, curcutsize, curvolume, it->second, conductance[i]);
	}
	// we stopped the iteration when it finished, or when it hit target_vol
	size_t lastind = i;
	double mincond = std::numeric_limits<double>::max();
	size_t mincondind = 0; // set to zero so that we only add one vertex 
	for (i = 0; i<lastind; i++) {
		if (conductance[i] < mincond) {
			mincond = conductance[i];
			mincondind = i;
		}
	}
	//printf("mincond=%f mincondind=%i\n", mincond, mincondind);
	if (lastind == 0) {
		// add a case 
		mincond = 0.0;
	}
	i = 0;
	for (vertex_prob_type::iterator it = prpairs.begin(), itend = prpairs.end();
		it != itend && i<mincondind + 1; ++it, ++i) {
		cluster.push_back(it->first);
	}
	if (outcond) { *outcond = mincond; }
	if (outvolume) { *outvolume = volume[mincondind]; }
	if (outcut) { *outcut = cutsize[mincondind]; }
}


/** Cluster will contain a list of all the vertices in the cluster
* @param set the set of starting vertices to use
* @param alpha the value of alpha in the PageRank computation
* @param target_vol the approximate number of edges in the cluster
* @param p the pagerank vector
* @param r the residual vector
* @param a vector which supports .push_back to add vertices for the cluster
* @param stats a structure for statistics of the computation
*/
template <class Queue>
int hypercluster_pagerank_multiple(const sparserow &G,
	const std::vector<unsigned int>& set, double alpha, double target_vol,
	sparsevec& p, sparsevec &r, Queue& q,
	std::vector<unsigned int>& cluster, local_pagerank_stats *stats)
{
	// reset data
	p.map.clear();
	r.map.clear();
	q.empty();

	assert(target_vol > 0);
	assert(alpha < 1.0); assert(alpha > 0.0);

	//r.map[start] = 1.0;
	//double maxdeg = 0;
	int maxdeg = 0;
	for (size_t i = 0; i<set.size(); ++i) {
		assert(set[i] >= 0); assert(set[i] < G.nodeSize);
		r.map[set[i]] = 1. / (double)(set.size());
		maxdeg = std::max(maxdeg, sr_degree(G, set[i]));
	}



	//double pr_eps = 1.0/std::max((double)sr_volumn(G,start)*(double)target_vol, 100.0);
	//double pr_eps = std::min(1.0/std::max(10.*target_vol, 100.0), 
	//1./(double)(set.size()*maxdeg + 1));
	double pr_eps = 1.0 / std::max(10.0*target_vol, 100.0);
	if (stats) { stats->eps = pr_eps; }

	//printf("find_cluster: target_vol=%7lli alpha=%5.3ld pr_eps=%ld\n", target_vol, alpha, pr_eps);

	// calculate an integer number of maxsteps
	double maxsteps = 1. / (pr_eps*(1. - alpha));
	maxsteps = std::min(maxsteps, 0.5*(double)std::numeric_limits<int>::max());

	int nsteps = compute_local_pagerank(G, r, p, alpha, pr_eps, (int)maxsteps, q);
	if (nsteps == 0) {
		p = r; // just copy over the residual
	}
	int support = r.map.size();
	if (stats) { stats->steps = nsteps; }
	if (stats) { stats->support = support; }

	//mexPrintf("setsize=%zu, nsteps=%i, support=%i\n", set.size(), nsteps, support);

	// scale the probablities by their degree
	for (sparsevec::map_type::iterator it = p.map.begin(), itend = p.map.end();
		it != itend; ++it) {
		it->second *= 1.0 / (double)std::max(sr_degree(G, it->first), 1);
	}

	double *outcond = NULL;
	double *outvolume = NULL;
	double *outcut = NULL;
	if (stats) { outcond = &stats->conductance; }
	if (stats) { outvolume = &stats->volume; }
	if (stats) { outcut = &stats->cut; }
	cluster_from_sweep(G, p, cluster, outcond, outvolume, outcut);
	return (0);
}

local_pagerank_stats pprgrow(const sparserow & G, std::vector<unsigned int>& clusters, double alpha,
	double targetvol, double &fcond, double &fcut,
	double &fvol, sparsevec &p)
{
	sparsevec r;
	std::queue<unsigned int> q;
	local_pagerank_stats stats;
	std::vector<unsigned int> bestclus;
	hypercluster_pagerank_multiple(G, clusters, alpha, targetvol,
		p, r, q, bestclus, &stats);
	clusters = bestclus;
	fcond = stats.conductance;
	fcut = stats.cut;
	fvol = stats.volume;
	return stats;
}

void copy_array_to_index_vector(const unsigned int* v, unsigned int length, std::vector<unsigned int>& vec)
{

	for (unsigned int i = 0; i < length; i++) {
		vec.push_back(v[i]);
	}
}

/*
*Input: alpha, targetvol
*
*/
local_pagerank_stats ppr(const double alpha, const double targetvol,const sparserow &r,const vector<unsigned int> seeds,vector<unsigned int> &cluster)
{
	//double alpha = 0.99;
	//double targetvol = 1000;

	//r.ai = mxgetjc(mat);//number of nonzero elements in each colomns
	//r.aj = mxgetir(mat);//nonzero elements index 
	//r.a = mxgetpr(mat);//real number
	
	double cond, cut, vol;
	sparsevec p;
	local_pagerank_stats stats = pprgrow(r, cluster, alpha, targetvol,
		cond, cut, vol, p);
	return stats;
	
}

//one single query node
local_pagerank_stats ppr(const double alpha,const double targetvol, const sparserow &r, const unsigned int seed, vector<unsigned int> &cluster)
{
	
	vector<unsigned int> seeds;
	seeds.push_back(seed);
	cluster.clear();
	cluster.push_back(seed);
	return ppr(alpha, targetvol, r, seeds, cluster);
}