#include "global.h"

string inputfasta;
string oprefix;
string logfile;
fstream LOG;
string rtype = "BIPES";
int hbcut = 3;
int endgap = 1;
int ksize = 6;
double kcut = 0.5;
int thread = 4;
int scut = 10;
float distcut = 0.03;
int otunumberindex = 0;
string method = "cl";
int gapopen = -10;
int gapextent = -1;
int mmt = -1;
int match = 1;

int newotu = 1;

//store the hbcut index and iterator
int hbcount = 0; 

vector <pair<string,int> > tvector;
/*
hash_map <string, vector<char>, str_hash, compare_str> seqTab;
hash_map <string, vector<int>, str_hash, compare_str> seqCode;
hash_map <string, int, str_hash, compare_str> seqLen;
*/
vector < vector<char> > seqTab;
vector < vector<int> > seqCode;
vector < int > seqLen;


map <int, vector <pair<int,double> > > hbneedle;
map <int, vector <int> > hbseed;

hash_map <string, vector<int>, str_hash, compare_str> otu;
map <int, string> id_otu;
hash_map <string, int, str_hash, compare_str> seq_id;
hash_map <int, string> id_seq;
list <string> tlist;

double ktime = 0;
double ntime = 0;
