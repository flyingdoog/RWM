#ifndef _HEAD_SIMILARITIES
#define _HEAD_SIMILARITIES


#include"tools.h"
#include"Graph.h"
float PospearsonCorrelation(map<uint, double>& m1, map<uint, double>& m2);
float pearsonCorrelation(map<uint, double>& m1, map<uint, double>& m2);
double pearsonCorrelation(double *m1, double * m2, size_t size);
double cosineSimilarity(map<uint, double>& m1, map<uint, double>& m2, uint seed, double thres);
double edgeCorrectness(sparserow &graph, map<uint, double>& m1, map<uint, double>& m2, uint seed, double thres);
double ICS(sparserow &graph, map<uint, double>& m1, map<uint, double>& m2, uint seed, double thres);
float PospearsonCorrelation(uint seed, map<uint, double>& m1, map<uint, double>& m2);
double cosineSimilarity(map<uint, double>& m1, map<uint, double>& m2, uint seed);
double cosineSimilarity(map<uint, double>& m1, map<uint, double>& m2, map<uint, double>& v1);
double cosineSimilarity(map<uint, double>& m1, map<uint, double>& m2);

#endif // !_HEAD_SIMILARITIES

