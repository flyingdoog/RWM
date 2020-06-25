#pragma once
#ifndef _tools_head
#define _tools_head
#include<string>
#include<vector>
#include<time.h>
#include <fstream>
#include<vector>
#include <cstring>
#include<iostream>
#include<algorithm>
#include <list>
#include<map>
#include<set>
#include<math.h>
#include<stack>
#include<sstream>
#include<queue>
#include"debug.h"
typedef unsigned int uint;
using namespace std;
#define MINSCORE 1e-8
#define ERROR_CODE -1234567
#define MINI(x,y) (((x)<(y))? (x):(y))
#define MAX(x,y) (((x)>(y))? (x):(y))

enum Wmodel { MODEL_PC, MODEL_COSINE, MODEL_EC, MODEL_ICS};

struct int_double_pair_sorted_by_score {
	uint id;
	double score;
	bool operator< (int_double_pair_sorted_by_score i) { return (score>i.score); }
};

void sort(vector<int_double_pair_sorted_by_score> &sortedxs, const map<uint,double> &x);
void getcontext(string& origin);
void split(vector<string>&, string &origin, char c);
void split(vector<string>&, string &origin, set<char> &sps);
void split(set<string>&, string &origin, char c);
double sumArray(const map<uint, double> &xs1);
void delNumber(string &);
string norm(string &origin);
bool equals(string &s1, string& s2);
string normtitle(string &tword);
double mabs(const map<uint, double> &, const map<uint, double> &);
double mabs(const map<uint, double> *xs1, const map<uint, double> *, const  uint K);
void copyArray(map<uint, double> &to, const map<uint, double>&from);
void catAndNormArray(map<uint, double>& to,const map<uint, double>& from);
void norm(map<uint, double>&);
uint getUint(const char*,size_t size, uint&);
uint getUint(string &, uint&);
double mabs(const map<uint, double> &m1, const map<uint, double> &m2);
#endif