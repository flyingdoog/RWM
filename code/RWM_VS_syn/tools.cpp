#include"stdafx.h"
#include"tools.h"
void getcontext(string& origin) {
	const char * chs = origin.c_str();
	int begin = 0, end = 0;
	for (unsigned int i = 0; i<origin.length(); i++)
		if (chs[i] == '>') {
			begin = i + 1;
			break;
		}
	for (unsigned int i = begin; i<origin.length(); i++)
		if (chs[i] == '<') {
			end = i;
			break;
		}
	origin = string(origin.substr(begin, (end - begin)));
}

void split(vector<string>& result, string &origin, char c) {
	if (origin.length() <= 1)
		return;
	int begin = 0;
	int next = 0;
	int end = origin.length();
	while (next<end) {
		if (origin[next] == c) {
			if (next - begin>0)
				result.push_back(origin.substr(begin, (next - begin)));
			begin = next + 1;
		}
		next++;
	}
	if (next - begin>0)
		result.push_back(origin.substr(begin, (next - begin)));
}

void split(set<string>& result, string &origin, char c) {
	if (origin.length() <= 1)
		return;
	int begin = 0;
	int next = 0;
	int end = origin.length();
	while (next<end) {
		if (origin[next] == c) {
			if (next - begin>0)
				result.insert(origin.substr(begin, (next - begin)));
			begin = next + 1;
		}
		next++;
	}
	if (next - begin>0)
		result.insert(origin.substr(begin, (next - begin)));
}


void split(vector<string> &result, string& origin, set<char> &sps) {
	int begin = 0;
	int next = 0;
	int end = origin.length();
	while (next<end) {
		if (sps.find(origin[next]) != sps.end()) {
			result.push_back(origin.substr(begin, (next - begin)));
			begin = next + 1;
		}
		next++;
	}
	result.push_back(origin.substr(begin, (next - begin)));
}

string normtitle(string &origin) {
	const char * chs = origin.c_str();
	char * res = new char[origin.length() + 1];
	int j = 0, i = 0;
	for (; i<origin.length(); i++) {
		if ((chs[i] >= 'a'&&chs[i] <= 'z') || chs[i] == '-')
			res[j++] = chs[i];
		else if ((chs[i] >= 'A'&&chs[i] <= 'Z'))
			res[j++] = chs[i] - 'A' + 'a';
	}
	for (j--; j>0; j--) {
		if (res[j] != 's')
			break;
	}
	res[j + 1] = 0;
	string re = string(res);
	delete[] res;
	return re;
}

string norm(string &origin) {
	const char * chs = origin.c_str();
	char * res = new char[origin.length() + 1];
	int j = 0, i = 0;
	for (i = 0; i<origin.length(); i++) {
		if ((chs[i] >= 'a'&&chs[i] <= 'z') || (chs[i] >= 'A'&&chs[i] <= 'Z'))
			break;
	}
	for (; i<origin.length(); i++) {
		if (chs[i] >= 'a'&&chs[i] <= 'z')
			res[j++] = chs[i];
		else if ((chs[i] >= 'A'&&chs[i] <= 'Z'))
			res[j++] = chs[i] - 'A' + 'a';
		else if (chs[i] == ' '&&res[j - 1] != ' ')
			res[j++] = ' ';
	}
	for (j--; j>0; j--) {
		if (res[j] != ' ')
			break;
	}
	res[j + 1] = 0;
	string re = string(res);
	delete[] res;
	return re;
}

bool equals(string &s1, string &s2) {
	if (s1.empty() || s2.empty())
		return false;
	if (s1.length() != s2.length())
		return false;
	unsigned int index = 0;
	while (index++<s1.length())
		if (s1[index] != s2[index])
			return false;
	return true;


}


void getFiles(string path, vector<string> &names) {
	intptr_t hFile = 0;
	struct _finddata_t fileinfo;
	string p;
	int count = 0;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1) {
		do {
			// delete . and ..
			if (count++>1)
				names.push_back(fileinfo.name);
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}

void delNumber(string &orgin) {
	const char* chs = orgin.c_str();
	char * newchs = new char[orgin.length() + 1];
	int leng = orgin.length();
	int j = 0;
	for (int i = 0; i<leng; i++) {
		while (chs[i] <= '9'&&chs[i] >= '0') i++;
		newchs[j++] = chs[i];
	}
	newchs[j] = 0;
	orgin = string(newchs);
	delete[] newchs;
}

double mabs(const map<uint,double> &m1, const map<uint, double> &m2) {
	double sum = 0;
	map<uint, double>::const_iterator mite1 = m1.begin();
	map<uint, double>::const_iterator mite2 = m2.begin();
	
	while (mite1 != m1.end() && mite2!=m2.end()) {
		if (mite1->first == mite2->first) {
			sum += (mite1->second > mite2->second) ? mite1->second - mite2->second : mite2->second - mite1->second;
			mite1++;
			mite2++;
		}
		else if (mite1->first> mite2->first){
			sum += (mite2->second > 0) ? mite2->second : -1.0*mite2->second;
			mite2++;
		}
		else {
			sum += (mite1->second > 0) ? mite1->second : -1.0*mite1->second;
			mite1++;
		}
	}

	while (mite1 != m1.end()) {
		sum += (mite1->second > 0) ? mite1->second : -1.0*mite1->second;
		mite1++;
	}

	while (mite2 != m2.end()) {
		sum += (mite2->second > 0) ? mite2->second : -1.0*mite2->second;
		mite2++;
	}
	return sum;

}

double mabs(const map<uint, double> *xs1, const map<uint, double> *xs2, const  uint K) {
	double sum = 0;
	for (uint i = 0; i < K; i++) {
		sum += mabs(xs1[i], xs2[i]);
	}
	return sum;
}

void copyArray(map<uint, double> &to,const map<uint, double>&from) {
	to.clear();
	for (map<uint, double>::const_iterator mite = from.begin(); mite != from.end(); mite++)
		to[mite->first] = mite->second;
}

void catAndNormArray(map<uint, double>& to,const map<uint, double>& from) {
	map<uint, double>().swap(to);
	float sum = 0;
	for (map<uint, double>::const_iterator mite = from.begin(); mite != from.end(); mite++)
		if (mite->second > MINSCORE)
			sum += mite->second;

	for (map<uint, double>::const_iterator mite = from.begin(); mite != from.end(); mite++)
		if (mite->second>MINSCORE)
			to[mite->first] = mite->second / sum;
}


//normalize m
void norm(map<uint, double>& m) {
	map<uint, double> to;
	catAndNormArray(to, m);
	m.swap(to);
}



double sumArray(const map<uint, double> &xs) {
	double sum = 0;
	for (map<uint, double>::const_iterator mite = xs.begin(); mite != xs.end(); mite++)
		sum += mite->second;
	return sum;
}

uint getUint(const char* x, size_t size, uint &i) {
	uint number = 0;
	bool found = false;
	for (; i< size; i++) {
		if (x[i]<'0' || x[i] >'9') {
			if (!found)
				continue;
			else
				break;
		}
		found = true;
		number = number * 10 + (x[i] - '0');
	}
	return number;

}

uint getUint(string & x, uint& index) {
	return getUint(x.c_str(), x.size(), index);
}

void sort(vector<int_double_pair_sorted_by_score> &sortedxs, const map<uint, double> &x) {
	sortedxs.clear();
	for (map<uint, double>::const_iterator it = x.begin(); it != x.end(); ++it) {
		int_double_pair_sorted_by_score mp;
		mp.id = it->first;
		mp.score = it->second;
		sortedxs.push_back(mp);
	}

	sort(sortedxs.begin(), sortedxs.end());
}
