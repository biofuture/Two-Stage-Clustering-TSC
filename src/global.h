#ifndef GLOBAL_H
#define GLOBAL_H
 
 #include "common.h"

extern string inputfasta;
extern string oprefix;
extern string logfile;
extern fstream LOG;
extern  int hbcut;
extern  int endgap;
extern  int ksize;
extern  double kcut;
extern  int thread;
extern  int scut;
extern  float distcut;
extern   string rtype;
extern int otunumberindex;
extern string method;
extern int gapopen;
extern int gapextent;
extern int match;
extern int mmt;

extern int hbcount; 
extern int newotu;

extern vector <pair<string,int> > tvector;
/*
extern hash_map <string, vector<char>, str_hash, compare_str> seqTab;
extern hash_map <string, vector<int>, str_hash, compare_str> seqCode;
extern hash_map <string, int, str_hash, compare_str> seqLen; 
*/
extern vector < vector<char> > seqTab;
extern vector < vector<int> > seqCode;
extern vector < int > seqLen;

extern map <int, vector <pair<int,double> > > hbneedle;
extern map <int, vector <int> > hbseed;

extern hash_map <string, vector<int>, str_hash, compare_str> otu;
extern map <int, string> id_otu;
extern hash_map <string, int,str_hash, compare_str> seq_id;
extern hash_map <int, string> id_seq;
extern list <string> tlist;

extern double ktime;
extern double ntime;

#endif 
