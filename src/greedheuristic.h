#ifndef GREEDHEURISTIC_H
#define GREEDHEURISTIC_H

#include "global.h"
#include "kmerneedle.h"

void greedyheuclust_excute();
void preclust();
void processlowabseq();
void outputclust();
void split_str(string pBuf, char szSub, vector<string> &vsArgv);

#endif
