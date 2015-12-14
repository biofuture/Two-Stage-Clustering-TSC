#include "kmerneedle.h"
/*---------seqab lib---------------*/
#include <seqan/align.h>
#include <seqan/graph_algorithms.h>
#include <seqan/index.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

using namespace seqan;

int getKmerNumber(string sequence, int index, int kmerSize, int kmern)
{
     int power4s[14] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864 };

        int kmerNum = 0;
        if (index == 0)
        {
                for(int i=0;i<kmerSize;i++){
                        if(toupper(sequence[i+index]) == 'A')           {       kmerNum += (0 * power4s[kmerSize-i-1]); }
                        else if(toupper(sequence[i+index]) == 'C')      {       kmerNum += (1 * power4s[kmerSize-i-1]); }
                        else if(toupper(sequence[i+index]) == 'G')      {       kmerNum += (2 * power4s[kmerSize-i-1]); }
                        else if(toupper(sequence[i+index]) == 'U')      {       kmerNum += (3 * power4s[kmerSize-i-1]); }
                        else if(toupper(sequence[i+index]) == 'T')      {       kmerNum += (3 * power4s[kmerSize-i-1]); }
                }
        }
		else if (index > 0)
        {
                if(toupper(sequence[index-1]) == 'A')           {       kmern -= (0 * power4s[kmerSize-1]);     }
                else if(toupper(sequence[index-1]) == 'C')      {       kmern -= (1 * power4s[kmerSize-1]);     }
                else if(toupper(sequence[index-1]) == 'G')      {       kmern -= (2 * power4s[kmerSize-1]);     }
                else if(toupper(sequence[index-1]) == 'U')      {       kmern -= (3 * power4s[kmerSize-1]);     }
                else if(toupper(sequence[index-1]) == 'T')      {       kmern -= (3 * power4s[kmerSize-1]);     }
                kmern = kmern*4;
                if(toupper(sequence[index+kmerSize-1]) == 'A')          {       kmern += (0 * power4s[0]);      }
                else if(toupper(sequence[index+kmerSize-1]) == 'C')     {       kmern += (1 * power4s[0]);      }
                else if(toupper(sequence[index+kmerSize-1]) == 'G')     {       kmern += (2 * power4s[0]);      }
                else if(toupper(sequence[index+kmerSize-1]) == 'U')     {       kmern += (3 * power4s[0]);      }
                else if(toupper(sequence[index+kmerSize-1]) == 'T')     {       kmern += (3 * power4s[0]);      }
                kmerNum = kmern;
        }

        return kmerNum;
}

void do_kmerneedle_set(map<int, vector< pair<int,double> > > & needlemap)
{
	vector < pair<string, int> >::iterator itvec = tvector.begin();
	vector < pair<string, int> >::iterator cur = itvec;
	//vector < pair<string, int> >::iterator ptr = itvec;
	vector < pair<int, double> > tempneedle;
	vector < pair<int, double> > tempkmer;

/*
	fstream OKMER, ONEEDLE;
	string okmer = oprefix + ".hb.kdist";
	string oneedle = oprefix + ".hb.ndist";

	OKMER.open(okmer.c_str(), ios::out);
	ONEEDLE.open(oneedle.c_str(), ios::out);
	if(OKMER == NULL || ONEEDLE == NULL)
	{
		cout << "can not write " << endl;
		exit(1);
	}
*/
 	//cout << hbcount << endl;
	while(itvec->second >= hbcut && itvec != tvector.end())
	{
		//ptr=itvec;
		itvec++;
		hbcount++;
	}//while
	
	//cout << hbcut <<  "\t" << hbcount <<  "\t" <<  tvector.size() << endl;
	//store the iterator 
	//cout << hbcount <<  "\thbcount" << endl;
	int qlen;
	vector <char> KmerTab1;
	vector <int> KmerCode1;
	int index = 1;
	if(cur != itvec)
	{
		for(; cur != itvec; cur++)
		{
			//Kmer needle
			string queryseq = cur->first;

			qlen = seqLen[index-1];
			KmerTab1 = seqTab[index-1];
			KmerCode1 = seqCode[index-1];

			//cout << hbcount << "\t" << index  << endl;
			tempneedle.resize(hbcount-index);
			tempkmer.resize(hbcount-index);

			#pragma omp parallel for num_threads(thread)
			for(int k = 0; k < hbcount - index; k++)
			{
				string objectseq = tvector[index+k].first;
				int slen;
				vector <char> KmerTab2;
				vector <int> KmerCode2;
				slen = seqLen[index+k];
				KmerTab2 = seqTab[index+k];
				KmerCode2 = seqCode[index+k];

				int minlen = MIN(qlen,slen);
				minlen = minlen - ksize + 1;
				int uCommoncount = 0;
				for(int m = 0; m < KmerCode1.size(); ++m)
				{
					int kmernum = KmerCode1[m];
					int count1 = KmerTab1[kmernum];
					int count2 = KmerTab2[kmernum];
					uCommoncount += MIN(count1, count2);
				}
				double kdist = 1 - (double)uCommoncount/minlen;	
				
				//tempkmer[k] = make_pair(index+k,kdist);

				if(kdist <= kcut)
				{
					//caculate the needledist
					double score = NW_globaAlignment(queryseq,objectseq);
					tempneedle[k] = make_pair(index+k, score);
				}
				else
				{
					tempneedle[k] = make_pair(index+k,1);
				}
			}

			//read needleman dist to map
			for(vector <pair<int, double> >::iterator cr = tempneedle.begin(); cr != tempneedle.end(); cr++)
			{
				if(cr->second != 1)
				{
					needlemap[index-1].push_back(*cr);
	
				//	ONEEDLE << index-1 << "\t" << cr->first << "\t" << fixed << setprecision(6) << cr->second << endl;
				}
			}
/*
			for(vector <pair<int, double> >::iterator cr = tempkmer.begin(); cr != tempkmer.end(); cr++)
			{
					OKMER << index-1 << "\t" << cr->first << "\t" << cr->second << endl;
			}
*/			
			index++;

		}
	}
}//do_kmeneedleset


inline double  cmpdouble(const pair<int,double>& x, const pair<int,double>& y)   
{   
	return x.second < y.second;   
}//cmp


void do_kmerneedle_one2set(string qseq, vector < pair<int,double> > & needlevect)
{
		//caculate all the kmer
		int qid = seq_id[qseq];
		//vector < pair<int,double> > tmpOTU;
		map <double, vector<int> > tmpOTU;
		vector < pair<int,double> > cellOTU;
		vector < int > temp;
		int Step = 1600;
		int i = 0;
		
		int qlen = seqLen[qid];
		vector <char> KmerTab1;
		vector <int> KmerCode1;

		KmerTab1 = seqTab[qid];
		KmerCode1 = seqCode[qid];

		time_t kst; time(&kst);
		for(hash_map <string, vector<int>, str_hash, compare_str>::iterator cur = otu.begin(); cur != otu.end(); cur++)
		{
			vector <int> oneotu = cur->second;
			for(vector<int>::iterator curr = oneotu.begin(); curr != oneotu.end(); curr++)
			{
				i++;
				if(i < Step )
				{
					temp.push_back(*curr);
				}
				
				if(i == Step)
				{
					cellOTU.resize(Step);
					temp.push_back(*curr);	
					#pragma omp parallel for num_threads(thread)
					for(int in=0; in< Step; in++)
					{
						int sid = temp[in];
						string sseq = id_seq[sid];
						int slen;
						vector <char> KmerTab2;
						vector <int> KmerCode2;
						slen = seqLen[sid];
						KmerTab2 = seqTab[sid];
						KmerCode2 = seqCode[sid];

						int minlen = MIN(qlen,slen);
						minlen = minlen - ksize + 1;
						int uCommoncount = 0;
						for(int m = 0; m < KmerCode1.size(); ++m)
						{
							int kmernum = KmerCode1[m];
							int count1 = KmerTab1[kmernum];
							int count2 = KmerTab2[kmernum];
							uCommoncount += MIN(count1, count2);
						}
						double kdist = 1 - (double)uCommoncount/minlen;	
						cellOTU[in] = make_pair(sid, kdist);
					}
					
					for(vector <pair<int,double> >::iterator cu = cellOTU.begin();
							cu != cellOTU.end(); cu++)
					{
						if(cu->second <= kcut)
						{
							//tmpOTU.push_back(*cu);
							tmpOTU[cu->second].push_back(cu->first);
						}
					}
					temp.clear();
					cellOTU.clear();
					i=0;
				}//if
			}//for one otu
		}//for all otu
					
		if(i != 0 &&  i< Step)
		{
			cellOTU.resize(i);
			#pragma omp parallel for num_threads(thread)
			for(int in=0; in< i; in++)
					{
						int sid = temp[in];
						string sseq = id_seq[sid];
						int slen;
						vector <char> KmerTab2;
						vector <int> KmerCode2;
						slen = seqLen[sid];
						KmerTab2 = seqTab[sid];
						KmerCode2 = seqCode[sid];

						int minlen = MIN(qlen,slen);
						minlen = minlen - ksize + 1;
						int uCommoncount = 0;
						for(int m = 0; m < KmerCode1.size(); ++m)
						{
							int kmernum = KmerCode1[m];
							int count1 = KmerTab1[kmernum];
							int count2 = KmerTab2[kmernum];
							uCommoncount += MIN(count1, count2);
						}
						double kdist = 1 - (double)uCommoncount/minlen;	
						cellOTU[in] = make_pair(sid, kdist);
					}
					
					for(vector <pair<int,double> >::iterator cu = cellOTU.begin();
							cu != cellOTU.end(); cu++)
					{
						if(cu->second <= kcut)
						{
							//tmpOTU.push_back(*cu);
							tmpOTU[cu->second].push_back(cu->first);
						}
					}
					temp.clear();
					cellOTU.clear();
		}//if

time_t ken; time(&ken);
ktime += difftime(ken,kst);

		//Kmer finish 
		//start needleman
		int flag = 0;
		if(!tmpOTU.empty())
		{
			//sort(tmpOTU.begin(),tmpOTU.end(),cmpdouble);
			int cnum=0;
			//for(vector <pair<int, double> >::iterator cur = tmpOTU.begin(); cur != tmpOTU.end(); cur++)
			for(map <double, vector<int> >::iterator cur = tmpOTU.begin(); cur != tmpOTU.end(); cur++)
			{
			
				if(cnum >= scut)
				{
					break;
				}
				
				flag = 0;
				vector <int> tmk = cur->second;
				for(vector <int>::iterator curr = tmk.begin(); curr != tmk.end(); ++curr)
				{
					double score = NW_globaAlignment(qseq, id_seq[*curr]);
					
					if(score <= distcut)
					{
						//cout << cur->first << endl;
						needlevect.push_back(make_pair(*curr, score));
					}
					else
					{
						flag = 1;
					}
				}
				
				if(flag = 1)
				{
					++cnum;
				}

			}
			//cout << "This is END\\n"<< endl;
		}
time_t nen; time(&nen);
ntime += difftime(nen, ken);

}//do_kmerneedle_one2set


double NW_globaAlignment(string queryseq, string objseq)
{
	typedef String <Dna> TSequence;
	TSequence seq1 = queryseq;
	TSequence seq2 = objseq;
	double idscore;
	//Score < int > score(match,mism,gape,gapo);
	if(rtype == "454")
	{
		gapextent = 0;
	}

	Score < int > score(match,mmt,gapextent,gapopen);

	typedef StringSet<TSequence, Dependent<> > TStringSet;
	typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;
	TStringSet string_set;
	appendValue(string_set, seq1);
	appendValue(string_set, seq2);
	TAlignmentGraph alignment_graph(string_set);
	String<char> align_matrix;
	double sc = globalAlignment(alignment_graph, score,AlignConfig<false,false,true,true>(), NeedlemanWunsch());
	convertAlignment(alignment_graph, align_matrix);

	seqan::Iterator<seqan::String<char>, seqan::Rooted >::Type it2 = begin(align_matrix);
	string str1,str2;
	for (goBegin(it2); position(it2)<length(align_matrix)/2; goNext(it2))
	{
		str1 += value(it2);
		str2 += value(it2+(length(align_matrix)/2));
	}

//	cout << str1 << endl << str2 << endl;
	if(rtype == "BIPES")
	{
		double j = 0;
		for (int i = 0; i < str1.size(); i++)
		{
			if (str1[i] != str2[i] || str1[i] == '-' || str2[i] == '-')
			{
				++j;
			}
		}
		if (length(seq1) > length(seq2))
		{
			idscore =  j/length(seq1);
		}
		else
		{
			idscore =  j/length(seq2);
		}
	//	cout << idscore << endl;
		return idscore;
	}       
	else
	{
		//trim begin and end gap 
		double gap_alert1 = 0;
		double gap_alert2 = 0;
		double residuecount = 0;
		double distance = 0;

		int start_seq1 = 0;
		int start_seq2 = 0;
		int end_seq1 = length(str1);
		int end_seq2 = length(str2);
		if(end_seq1 != end_seq2)
		{
			cout << "wrong\n" <<  str1
				<< "\n" << str2 << endl;
			exit(1);
		}
		if(endgap == 0)
		{
			int p = 0;
			while(str1[p] == '-') p++;
			start_seq1 = p;

			p = 0;
			while(str2[p] == '-') p++;
			start_seq2 = p;

			p = end_seq1 - 1;
			while(str1[p] == '-') p--;
			end_seq1 = p;

			p = end_seq2 - 1;
			while(str2[p] == '-') p--;
			end_seq2 = p;
		}

		int start_seq = start_seq1;
		if(start_seq < start_seq2)
		{
			start_seq = start_seq2;
		}

		int end_seq = end_seq1;
		if(end_seq > end_seq2)
		{
			end_seq = end_seq2;
		}


		for (int k = start_seq; k<=end_seq; k++)
		{
			if(str1[k] == '-' && str2[k] == '-')
			{
				continue;
			}

			if(str1[k] == 'N' && str2[k] == 'N')
			{
				gap_alert1 = 0;
				gap_alert2 = 0;
				continue;
			}


			if(str1[k]== '.' && str2[k] == '.')
			{
				gap_alert1 = 0;
				gap_alert2 = 0;
				continue;
			}

			if (str1[k] == '-')
			{
				if(gap_alert1 == 0)
				{
					residuecount+=1.0;
					distance+=1.0;
				}
				gap_alert1 = 1;
				gap_alert2 = 0;
				continue;
			}

			if(str2[k]== '-')
			{
				if(gap_alert2 == 0)
				{
					residuecount+=1.0;
					distance+=1.0;
				}
				gap_alert1 = 0;
				gap_alert2 = 1;
				continue;
			}

			if(str1[k] != str2[k] )
			{
				distance += 1.0;
				residuecount+=1.0;
				gap_alert1 = 0;
				gap_alert2 = 0;
				continue;

			}

			if(str1[k] == str2[k] )
			{
				residuecount+= 1.0;
				gap_alert1 = 0;
				gap_alert2 = 0;
				continue;
			}

		}

		if(residuecount > 0)
		{
			idscore = distance / residuecount;
		}
	//		cout << idscore << endl;
		return idscore;
	}
}
