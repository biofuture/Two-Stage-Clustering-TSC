#include "global.h"
#include "kmerneedle.h"
#include "herarchicalcluster.h"
#include "greedheuristic.h"


void usage()
{
	cout 	<< "\tTwo Stage Clustering Version 1.1 CMS by GIT\n"
		<< "\tLast Modified 22/11/2011\n"
		<< "\tby Jiang Xiao-Tao\n"
		<< "\thust.cn@163.com\n"
		<< "\tUsage: ./TSC -i <.fasta> -o <out.prefix> -c [abundancecutoff]\n" 
		<< "\t-x [terminal gap] -n [numberofthread] -k [ksize] -f [kmerfilter]\n"
		<< "\t-s [NotMatchCut] -d [distance] -m [al/cl/sl]\n\n"
		<< "\t-i  <str>   input fasta file\n"
		<< "\t-o  <str>   output file prefix\n"
		<< "\t-c [number] the cutoff for high and low abundant default 3\n"
		<< "\t-g [gapopen] gapopen penality default -10\n"
		<< "\t-e [gapextension] default -1\n"
		<< "\t-a [match] substiution matrix match default 1\n"
		<< "\t-b [mismatch] substitution matrix mismatch default -1\n"
		<< "\t-k [number] the word size for kmer default 6\n"
		<< "\t-f [double] the kmer filter default 0.5\n"
		<< "\t-d [double] the distance to do OTU picking default 0.03\n"
		<< "\t-n [number] number of thread used default 4\n"
		<< "\t-x [0/1] for 454 data if 1 then caculate end gap\t if 0 remove end gap default 1\n"
		<< "\t-s [number] the kind number of the kmerdist not fit the distance require default 10\n"
		<< "\t-m [method] the clust method for high abundant sequences\n"
		<< "\t-r [data type] the data type to process deafult illumina data if using 454 data set -r 454\n"
		<< endl;
	exit(1);
}

/***-------------------------------------------------------------***/

int cmpint(const pair<string,int>& x, const pair<string,int>& y)
{
	return x.second > y.second;
}//

void unique(string fastr, string  opre, int cutoff)
{
	fstream IFA,OUQ,ONA;
	string	opreuq = opre + ".unique.fa";
	string	oprena = opre + ".unique.name";
	IFA.open(fastr.c_str(), ios::in);
	OUQ.open(opreuq.c_str(), ios::out);
	ONA.open(oprena.c_str(), ios::out);	

	if(IFA == NULL)
	{
		cout << "input fasta does not exist!"  << endl;
		exit(1);
	}

	if(OUQ == NULL)
	{
		cout << "can not write in this directory 1 !"  << endl;
		exit(1);
	}

	if(ONA == NULL)
	{
		cout << "can not write in this directory 2 !"  << endl;
		exit(1);
	}

	//caculate the diffr for this dataset
	int minlen = 100000;
	int maxlen = 0;
	int avelen = 0;
	int totallen = 0;
	int totalnum = 0;

	string id, seq;
	hash_map <string, vector<string>, str_hash, compare_str> uniq;
	while(getline(IFA, id ,'\n'))
	{
		id = id.erase(0,1);
		getline(IFA, seq, '\n');
		uniq[seq].push_back(id);
	}//while		
	id.clear();
	seq.clear();	

	int seqlens = 0;
	vector <string> namevec;
	hash_map <string, vector<string>, str_hash, compare_str>::iterator uqit = uniq.begin();
	for(;uqit != uniq.end();++uqit)
	{
		seqlens = (uqit->first).length();
		totallen += seqlens;//
		++totalnum;//
		if(minlen > seqlens)
		{
			minlen = seqlens;
		}

		if(maxlen < seqlens)
		{
			maxlen = seqlens;
		}

		namevec = uqit->second;	
		tvector.push_back(make_pair(uqit->first,namevec.size()));
		vector <string>::iterator cr = namevec.begin();
		ONA << *cr << "\t" <<  namevec.size() <<  "\t" << *cr;
		++cr;
		if(cr != namevec.end())
		{	
			for(; cr != namevec.end(); ++cr)
			{
				ONA << "," << *cr;
			}
		}
		ONA << endl;
		namevec.clear();
	}//for	
	//
	//
	avelen = floor(totallen / totalnum);

	LOG << "	Sequences State\n" << "	Maxlen\t" <<  maxlen <<  "	minlen\t" << minlen << "	average lenght\t" << avelen << endl;
	sort(tvector.begin(), tvector.end(), cmpint);	
	vector <pair<string,int> >::iterator tv = tvector.begin();
	int indexs =0 ;
	for(; tv != tvector.end(); ++tv)
	{
		OUQ << ">" << uniq[tv->first][0] << endl << tv->first << endl;

		//init the Kmer structure
		int sizes = pow(4.0, ksize);
		vector<char> KmerTab1(sizes);
		vector<int> KmerCode1;
		unsigned kmernum = 0;
		int len = tv->first.size();
		for(int k=0; k<len-ksize+1; ++k)
		{
			unsigned kmerNumber = getKmerNumber(tv->first, k, ksize, kmernum);
			++(KmerTab1[kmerNumber]);
			if (((int)KmerTab1[kmerNumber]) == 1)
			{
				KmerCode1.push_back(kmerNumber);
			}
			kmernum = kmerNumber;
		}
		/*
		   seqTab[tv->first] = KmerTab1;
		   seqCode[tv->first] = KmerCode1;
		   seqLen[tv->first] = tv->first.length();
		 */

		seq_id[tv->first] = indexs;
		id_seq[indexs] = tv->first;

		++indexs;
		seqTab.push_back(KmerTab1);
		seqCode.push_back(KmerCode1);
		seqLen.push_back(tv->first.length());

	}
	//cout << seqTab.size() << "\t" << seqCode.size() << endl;
	uniq.clear();
}//unique


int main(int argc, char *argv[])
{	
	if(argc < 3)
	{
		usage();
	}

	int opt;
	while( (opt=getopt(argc,argv,"i:o:c:x:n:k:f:s:d:m:g:e:a:b:r:")) != -1)
	{	
		switch(opt)
		{
			case 'i':inputfasta = optarg;
				 break;
			case 'o':oprefix = optarg;
				 break;
			case 'c':hbcut = atoi(optarg);
				 break;
			case 'x':endgap = atoi(optarg);
				 break;
			case 'n':thread = atoi(optarg);
				 break;
			case 'k':ksize = atoi(optarg);
				 break;
			case 'f':kcut = atof(optarg);
				 break;
			case 's':scut = atoi(optarg);
				 break;
			case 'd':distcut = atof(optarg);
				 break;
			case 'm':method = optarg;
				 break;
			case 'g':gapopen =  atoi(optarg);
				 break;
			case 'e':gapextent = atoi(optarg);
				 break;
			case 'b':mmt = atoi(optarg);
				 break;
			case 'a':match = atoi(optarg);
				 break;
			case 'r':rtype = atoi(optarg);
				 break;
			default:usage();
		}	
	}

	logfile = oprefix + ".log";
	LOG.open(logfile.c_str(), ios::out);
	if(LOG == NULL)
	{
		cout << "Cant not write"  << endl; exit(1);
	}
	time_t beg; time(&beg);
	time_t pstart; time(&pstart);
	struct tm * timeinfo; 
	timeinfo = localtime ( &pstart ); 
	string asc = asctime (timeinfo);

	LOG << "TSC begin: "  <<  asc  << endl;
	//unique and sort by abundance and then cut sequences into two part, init Kmer
	LOG << "1.	unique sort by abundance and init kmer data structure " << endl;
	unique(inputfasta,oprefix,hbcut);
	LOG << "	unique cost time " << difftime(time(NULL), pstart) << " s" << endl;

	//process high abundant sequences 
	//1.caculate the kmer needle dist for the high abundant sequences
	//2.denoise by 2% identity
	//3.do hcluster for the denoised seed
	time(&pstart);
	LOG << "2.	start to do herarchical clustering for the high abundant sequences" << endl;
	hc_excute();
	LOG << "	hc_excute total cost time " << difftime(time(NULL), pstart) << " s" << endl; 

	//process the low abundant sequences
	//1.do preculst (optional)
	//2.do heuristics clustering
	//3.output the otu

	time(&pstart);
	LOG << "3.	start to do greedy heuristic clustering for the low abundant sequences" << endl;
	greedyheuclust_excute();
	LOG << "	greedy clustering  total cost time " << difftime(time(NULL), pstart) << " s" << endl;

	time(&pstart);
	timeinfo = localtime ( &pstart ); 
	asc = asctime (timeinfo);
	LOG << "Total cost time " << difftime(pstart, beg) << "s" << endl;
	LOG << "TSC end: " << asc << endl;

	return 1;
}//main
