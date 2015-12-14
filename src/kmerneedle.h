#ifndef KMERNEEDLE_H
#define KMERNEEDLE_H

#include "global.h"



int getKmerNumber(string sequence, int index, int kmerSize, int kmern);
void do_kmerneedle_set(map <int, vector < pair<int, double> > > &);
void do_kmerneedle_one2set(string ,  vector < pair<int,double> > &);
double NW_globaAlignment(string queryseq, string objseq);
inline double  cmpdouble(const pair<int,double>& x, const pair<int,double>& y);

#endif
