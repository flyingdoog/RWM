#include"Similarities.h"
extern double alpha;


/*
*
*/
float pearsonCorrelation(map<uint, double>& m1, map<uint, double>& m2) {
	if (m1.size() * m2.size() == 0)
		return 0;
	float sumx = 0, sumy = 0, sumx2 = 0, sumy2 = 0, sumxy = 0;
	int n = 0;
	map<uint, double>::iterator mite1 = m1.begin();
	map<uint, double>::iterator mite2 = m2.begin();

	while (mite1 != m1.end() && mite2 != m2.end()) {
		if (mite1->first == mite2->first) {
			sumxy += mite1->second*mite2->second;
			sumx += mite1->second;
			sumx2 += (mite1->second)*(mite1->second);
			sumy += mite2->second;
			sumy2 += (mite2->second)*(mite2->second);

			mite1++;
			mite2++;
		}
		else if (mite1->first > mite2->first) {
			sumy += mite2->second;
			sumy2 += (mite2->second)*(mite2->second);
			mite2++;
		}
		else {
			sumx += mite1->second;
			sumx2 += (mite1->second)*(mite1->second);
			mite1++;
		}
		n++;
	}
	while (mite1 != m1.end()) {
		sumx += mite1->second;
		sumx2 += (mite1->second)*(mite1->second);
		n++;
		mite1++;
	}
	while (mite2 != m2.end()) {
		sumy += mite2->second;
		sumy2 += (mite2->second)*(mite2->second);
		n++;
		mite2++;
	}

	float dev = sqrt(n*sumx2 - sumx*sumx)*sqrt(n*sumy2 - sumy*sumy);
	if (dev == 0) {
		//cout << " devided by ZERO!!!!! @ pearsonCorrelation @ tools.cpp" << endl;
		//return ERROR_CODE;
		return 1;
	}
	return (n*sumxy - sumx*sumy) / dev;

}

float PospearsonCorrelation(map<uint, double>& m1, map<uint, double>& m2) {
	float res1 = pearsonCorrelation(m1, m2);
	//cout << m1.size() << "\t" << m2.size() << "\t"<<res1<<endl;
	if (res1 == ERROR_CODE)
		return ERROR_CODE;
	return (res1>0) ? res1 : 0;
}


float PospearsonCorrelation(uint seed, map<uint, double>& m1, map<uint, double>& m2) {
	double m1seed = m1.find(seed)->second;
	double m2seed = m2.find(seed)->second;
	m1.erase(seed);
	m2.erase(seed);
	float res1 = pearsonCorrelation(m1, m2);
	//cout << m1.size() << "\t" << m2.size() << "\t"<<res1<<endl;
	m1[seed] = m1seed;
	m2[seed] = m2seed;
	if (res1 == ERROR_CODE)
		return ERROR_CODE;
	return (res1>0) ? res1 : 0;
}



double pearsonCorrelation(double *m1, double * m2, size_t size) {
	float sumx = 0, sumy = 0, sumx2 = 0, sumy2 = 0, sumxy = 0;
	int n = 0;
	size_t mite = 0;

	while (mite != size) {
		sumxy += m1[mite] * m2[mite];
		sumx += m1[mite];
		sumx2 += m1[mite] * m1[mite];
		sumy += m2[mite];
		sumy2 += m2[mite] * m2[mite];
		mite++;
	}

	float dev = sqrt(n*sumx2 - sumx*sumx)*sqrt(n*sumy2 - sumy*sumy);
	if (dev == 0) {
		return 0;
	}
	return (n*sumxy - sumx*sumy) / dev;

}


/*
// norm m1 by delete the query node and keep only top thres.
*/
void norm(map<uint, double> &n1, map<uint, double> &m1, uint seed, double thres) {
	struct mpair {
		uint id;
		double score;
		bool operator< (mpair i) { return (score>i.score); }
	};
	n1.clear();
	vector<mpair> tep;
	for (map<uint, double>::iterator it = m1.begin(); it != m1.end(); ++it) {
		mpair mp;
		mp.id = it->first;
		mp.score = it->second;
		if (mp.id == seed)
			continue;
		tep.push_back(mp);
	}
	sort(tep.begin(), tep.end());

	double sum = 0;
	for (uint i = 0; i < tep.size(); i++) {
		sum += tep[i].score;
		n1[tep[i].id] = tep[i].score;
		if (sum > thres)
			break;
	}
}


double cosineSimilarity(map<uint, double>& m1, map<uint, double>& m2) {
	double sumxy = 0;
	double sumx2 = 0, sumy2 = 0;
	map<uint, double>::iterator mite1 = m1.begin();
	map<uint, double>::iterator mite2 = m2.begin();
	uint n = 0;
	while (mite1 != m1.end() && mite2 != m2.end()) {
		if (mite1->first == mite2->first) {
			sumxy += mite1->second*mite2->second;
			sumx2 += (mite1->second)*(mite1->second);
			sumy2 += (mite2->second)*(mite2->second);

			mite1++;
			mite2++;
		}
		else if (mite1->first > mite2->first) {
			sumy2 += (mite2->second)*(mite2->second);
			mite2++;
		}
		else {
			sumx2 += (mite1->second)*(mite1->second);
			mite1++;
		}
		n++;
	}
	while (mite1 != m1.end()) {
		sumx2 += (mite1->second)*(mite1->second);
		n++;
		mite1++;
	}
	while (mite2 != m2.end()) {
		sumy2 += (mite2->second)*(mite2->second);
		n++;
		mite2++;
	}

	if (sumx2*sumy2 == 0)
		return 0;
	return sumxy / sqrt(sumx2) / sqrt(sumy2);
}

double cosineSimilarity(map<uint, double>& m1, map<uint, double>& m2, uint seed, double thres) {
	map<uint, double> n1, n2;
	norm(n1, m1, seed, thres);
	norm(n2, m2, seed, thres);
	return cosineSimilarity(n1, n2);
}

double cosineSimilarity(map<uint, double>& m1, map<uint, double>& m2, map<uint, double>& v1) {
	addArray(m1, v1, alpha - 1);
	double res = cosineSimilarity(m1, m2);
	float sumx2 = 0;
	map<uint, double>::iterator mite1 = m1.begin();

	for (; mite1 != m1.end(); mite1++) {
		sumx2 += (mite1->second)*(mite1->second);
	}

	res -= sqrt(sumx2 / m1.size());
	addArray(m1, v1, 1 - alpha);
	return (res>0) ? res : 0;
}

double cosineSimilarity(map<uint, double>& m1, map<uint, double>& m2, uint seed) {
	double m1seed = m1.find(seed)->second;
	double m2seed = m2.find(seed)->second;
	m1.erase(seed);
	m2.erase(seed);
	float res1 = cosineSimilarity(m1, m2);
	//cout << m1.size() << "\t" << m2.size() << "\t"<<res1<<endl;
	m1[seed] = m1seed;
	m2[seed] = m2seed;

	float sumx2 = 0;
	map<uint, double>::iterator mite1 = m1.begin();

	for (; mite1 != m1.end(); mite1++) {
		sumx2 += (mite1->second)*(mite1->second);
	}

	res1 -= sqrt(sumx2 / m1.size());

	return (res1>0) ? res1 : 0;

}


double edgeCorrectness(sparserow &graph, map<uint, double>& m1, map<uint, double>& m2, uint seed, double thres) {
	map<uint, double> n1, n2;
	norm(n1, m1, seed, thres);
	norm(n2, m2, seed, thres);
	uint root = graph.nodeSize;
	set<uint> edges1;
	for (map<uint, double>::iterator vite = n1.begin(); vite != n1.end(); vite++) {
		uint id = vite->first;
		for (uint nbindex = graph.indexs[id]; nbindex < graph.indexs[id + 1]; nbindex++) {
			uint nb = graph.nbs[nbindex];
			if (n1.find(nb) == n1.end() || nb<id)
				continue;
			edges1.insert(id* root + nb);
		}
	}
	if (edges1.size() == 0)
		return 0;
	uint count = 0;
	for (map<uint, double>::iterator vite = n2.begin(); vite != n2.end(); vite++) {
		uint id = vite->first;
		for (uint nbindex = graph.indexs[id]; nbindex < graph.indexs[id + 1]; nbindex++) {
			uint nb = graph.nbs[nbindex];
			if (n2.find(nb) == n2.end() || nb<id)
				continue;
			uint key = id* root + nb;
			if (edges1.find(key) != edges1.end())
				count++;
		}
	}
	return 1.0*count / edges1.size();

	return 0;
}


double ICS(sparserow &graph, map<uint, double>& m1, map<uint, double>& m2, uint seed, double thres) {
	return edgeCorrectness(graph, m2, m1, seed, thres);
}
