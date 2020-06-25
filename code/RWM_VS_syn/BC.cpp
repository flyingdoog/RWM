/*
* Centrality Measures Calculator (Closeness, Node Betweenness, Edge Betweenness)
* Copyright (C) 2013 George Piskas
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*
* Contact: geopiskas@gmail.com
* Important notice: Works with Linux.
*/
#include"stdafx.h"
#include"BC.h"
#include <iostream>
#include <cstdio>
#include <fstream>

#include <vector>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <list>
#include"tools.h"
#include <cstring>
#include <sstream>
#include <algorithm>
#include <limits>
#include"Graph.h"
using namespace std;

// A simple neighbor struct, consisting of the target neighbor and the edge weight.
struct neighbor {
	int target;
	int weight;
	neighbor(int mTarget, int mWeight) : target(mTarget), weight(mWeight) {
	}
};

// The new adjacency list type.
typedef vector<vector<neighbor> > adjacency_list;

// Dijkstra's algorithm is used to calculate all the single source shortest paths in a weighted graph and the source's closeness.
float dijkstra_SSSP(int src, int n, stack<int> &visitStack, vector<int> &sigma, list<int> *pred, adjacency_list &adjList) {
	// Closeness counter.
	float closeness = 0;

	// Vector that holds the distances from the source.
	vector<int> dist;
	dist.resize(n, numeric_limits<int>::max());
	dist[src] = 0;

	// Queue used for the Dijkstra's algorithm.
	set<pair<int, int> > visitQueue;
	visitQueue.insert(make_pair(dist[src], src));

	// While there are still elements in the queue.
	while (!visitQueue.empty()) {
		// Pop the first.
		set<pair<int, int> >::iterator vPair = visitQueue.begin();
		int srcDist = vPair->first;
		int v = vPair->second;
		visitQueue.erase(vPair);
		visitStack.push(v);

		// Closeness part aggregation.
		closeness += srcDist;

		// Check the neighbors w of v.
		for (vector<neighbor>::iterator it = adjList[v].begin(); it != adjList[v].end(); it++) {
			int w = it->target;
			int newDist = srcDist + it->weight;
			// Node w found for the first time or the new distance is shorter?
			if (newDist < dist[w]) {
				visitQueue.erase(make_pair(dist[w], w));
				visitQueue.insert(make_pair(newDist, w));
				dist[w] = newDist;
				pred[w].clear();
				sigma[w] = 0;
			}
			// Is the shortest path to w via u?
			if (newDist == dist[w]) {
				pred[w].push_back(v);
				sigma[w] += sigma[v];
			}
		}
	}
	// Closeness part inversion.
	if (closeness != 0) {
		return 1.0 / closeness;
	}
	else {
		return 0;
	}
}

// BFS algorithm is used to calculate all the single source shortest paths in a non weighted graph and the source's closeness.
float bfs_SSSP(int src, int n, stack<int> &visitStack, vector<int> &sigma, list<int> *pred, adjacency_list &adjList) {
	// Closeness counter.
	float closeness = 0;

	// Vector that holds the distances from the source.
	vector<int> dist;
	dist.resize(n, -1);
	dist[src] = 0;

	// Queue used for the Bfs algorithm.
	queue<int> visitQueue;
	visitQueue.push(src);

	// While there are still elements in the queue.
	while (!visitQueue.empty()) {
		// Pop the first.
		int v = visitQueue.front();
		visitQueue.pop();
		visitStack.push(v);

		// Closeness part aggregation.
		closeness += dist[v];

		// Check the neighbors w of v.
		for (vector<neighbor>::iterator it = adjList[v].begin(); it != adjList[v].end(); it++) {
			int w = it->target;
			// Node w found for the first time?
			if (dist[w] < 0) {
				visitQueue.push(w);
				dist[w] = dist[v] + 1;
			}
			// Is the shortest path to w via u?
			if (dist[w] == dist[v] + 1) {
				pred[w].push_back(v);
				sigma[w] += sigma[v];
			}
		}

	}
	// Closeness part inversion.
	if (closeness != 0) {
		return 1.0 / closeness;
	}
	else {
		return 0;
	}
}

// Given two node indices, this function returns a string representation with
// the smallest index first, a dash, and the second index after.
string getEdgeTag(int n1, int n2) {
	ostringstream os;
	if (n1 <= n2) {
		os << n1 << "-" << n2;
	}
	else {
		os << n2 << "-" << n1;
	}
	return os.str();
}

// Reads an input file and fills up the adjacency list as well as the edges.
void readGraph(map<uint, uint> &newindex, adjacency_list &adjList, map<string, float> edgeList, set<uint> &scmty, sparserow & graph) {

	int e = 0; // Total number of edges (for statistics).
	size_t len = 0;
	adjList.reserve(scmty.size());
	int i = 0;
	for (set<uint>::iterator site = scmty.begin(); site != scmty.end();site++) {
		newindex[*site] = i++;
	}

	uint *indexs = graph.indexs;
	uint *nbs = graph.nbs;

	// Read the nodes and the edges, one by one, and fill up adjList and edgeBetweenness.
	int start, end, weight = 1;
	for (set<uint>::iterator site = scmty.begin(); site != scmty.end(); site++) {
		uint node = *site;
		if (newindex.find(node) == newindex.end())
			continue;
		start = newindex.find(node)->second;
		
		for (uint index = indexs[node]; index < indexs[node + 1]; index++) {
			uint nb = nbs[index];
			if (newindex.find(nb) != newindex.end()) {
				end = newindex.find(nb)->second;
				edgeList.insert(pair<string, float>(getEdgeTag(start, end), 0));
				adjList[start].push_back(neighbor(newindex.find(nb)->second, weight));
			}
		}
	}
}

// Clears the variables or re-initializes to 0, so that they are ready for the next loop.
void resetVariables(int src, int n, list<int> *pred, vector<int> &sigma, vector<float> &delta) {
	for (int i = 0; i < n; i++) {
		pred[i].clear();
	}

	sigma.clear();
	sigma.resize(n, 0);
	sigma[src] = 1;

	delta.clear();
	delta.resize(n, 0);
}

// Prints Closeness Centrality.
//void printCloseness(int n, vector<float> closeness, bool normalize) {
//	float nrml = 1;
//	if (normalize) {
//		nrml = 1.0 / (n - 1);
//	}
//	ofstream out;
//	out.open("out_closeness.txt");
//	cout << "> Closeness Centrality" << endl;
//	for (int i = 0; i < n; i++) {
//		cout << "Node " << i << ": " << closeness[i] / nrml << endl;
//		out << "Node " << i << ": " << closeness[i] / nrml << endl;
//	}
//	out.close();
//}

// Prints Edge Betweenness Centrality.
//void printEdgeBetweenness(int n, map<string, float> edgeBetweenness, bool normalize) {
//	float nrml = 1;
//	if (normalize) {
//		nrml = n * (n - 1);
//	}
//	ofstream out;
//	out.open("out_edge_tweenness.txt");
//	cout << endl << "> Edge Betweenness Centrality" << endl;
//	for (map<string, float>::iterator it = edgeBetweenness.begin(); it != edgeBetweenness.end(); it++) {
//		cout << "Edge " << it->first << ": " << it->second / nrml << endl;
//		out << "Edge " << it->first << ": " << it->second / nrml << endl;
//	}
//	out.close();
//}


int findCenterNode(set<uint> &cmty, sparserow & graph) {
	int n = cmty.size(); // Number of nodes
	bool isWeigthed = graph.weighted; // Weighted graph check.
	adjacency_list adjList; // Adjacency list.
	adjList.reserve(cmty.size());
	for (int i = 0; i < cmty.size(); i++) {
		vector<neighbor> temp;
		adjList.push_back(temp);
	}
	uint *indexs = new uint[cmty.size()];
	int i = 0;
	for (set<uint>::iterator site = cmty.begin(); site != cmty.end(); site++) {
		indexs[i++] = *site;
	}
	// Centrality measures.
	map<string, float> edgeBetweenness;
	vector<float> nodeBetweenness;
	nodeBetweenness.resize(n, 0);
	vector<float> closeness;
	closeness.resize(n, 0);
	map<uint, uint> newindex;
	// Input is read, and values are set to all the arguments.
	readGraph(newindex,adjList, edgeBetweenness,cmty, graph);

	list<int> *pred = new list<int>[n]; // List of predecessors of node v.
	vector<int> sigma;
	vector<float> delta;
	stack<int> visitStack; // Stack that holds the inverse order of visited nodes.

						   // For each node of the graph.
	for (int src = 0; src < n; src++) {
		// Prepare the variables for the next loop.
		resetVariables(src, n, pred, sigma, delta);

		if (isWeigthed) {
			// Closeness part. Using Dijkstra because graph is weighted.
			closeness[src] = dijkstra_SSSP(src, n, visitStack, sigma, pred, adjList);
		}
		else {
			// Closeness part.
			closeness[src] = bfs_SSSP(src, n, visitStack, sigma, pred, adjList);
		}

		// Get the inverse order of visited nodes.
		while (!visitStack.empty()) {
			int w = visitStack.top();
			visitStack.pop();

			// For each predecessors of node w, do the math!
			for (list<int>::iterator it = pred[w].begin(); it != pred[w].end(); it++) {
				int v = *it;
				float c = ((float)sigma[v] / (float)sigma[w]) * (1.0 + delta[w]);

				delta[v] += c;

				// Edge betweenness aggregation part.
				string tag = getEdgeTag(v, w);
				float tempC = edgeBetweenness[tag];
				edgeBetweenness.erase(tag);
				edgeBetweenness.insert(pair<string, float>(tag, tempC + c));
			}
			// Node betweenness aggregation part.
			if (w != src) {
				nodeBetweenness[w] += delta[w];
			}
		}

	}
	float nrml = 1;
	float maxf = -1;
	uint rslt = -1;
	nrml = (n - 1)*(n - 2);
	for (int i = 0; i < n; i++) {
		if (nodeBetweenness[i] > maxf) {
			maxf = nodeBetweenness[i];
			rslt = indexs[i];
		}
	}
	delete  indexs;
	
	return rslt;
}