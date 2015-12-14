#ifndef COMMON_H
#define COMMON_H

//The stream and io 
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <cstdlib>

//The vector and algorithm
#include <string>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <cmath>

//Time management
#include <time.h>
#include <ext/hash_map>



using namespace __gnu_cxx;
using namespace std;


#define MIN(x, y)       (((x) < (y)) ? (x) : (y))

struct str_hash{
	size_t operator()(const string& str) const
	{
		unsigned long __h = 0;
		for (size_t i = 0 ; i < str.size() ; i ++)
			__h = 5*__h + str[i];
		return size_t(__h);
	}
};

struct compare_str{
	bool operator()(const string p1, const string p2) const{
		return p1 == p2;
	}
};


#endif
