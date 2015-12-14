#include "greedheuristic.h"


/*split string by delimit*/
void split_str(string pBuf, char szSub, vector<string> &vsArgv)
{
	if (pBuf.empty())
	{
		//return -1;
		exit(-1);
	}

	vsArgv.clear();

	string strRaw = pBuf;
	string::size_type pPre = 0, pCur = 0;

	while((pCur = strRaw.find(szSub, pPre)) != string::npos)
	{
		string strSlice = strRaw.substr(pPre, pCur - pPre);
		if (!strSlice.empty())
		{
			vsArgv.push_back(strSlice);
		}
		pCur ++;
		pPre = pCur;
	}

	string strTail = strRaw.substr(pPre, strRaw.size() - pPre);
	if (!strTail.empty())
	{
		vsArgv.push_back(strTail);
	}

}

void greedyheuclust_excute()
{
	//do preclust to reduce the processed low abundant sequences
	LOG <<	"	[1] preclust start " << endl;
	time_t t1; time(&t1);
	preclust();
	LOG <<  "	cost time " << difftime(time(NULL), t1) << " s" << endl;

	//process the left low abundant sequences
	LOG <<	"	[2]. left low abundant sequences " << endl;
	time_t t2; time(&t2);
	processlowabseq();
	LOG <<  "	low ab cost time " << difftime(time(NULL), t2) << " s" << endl;
	//output the final clust result of the demanded distcut
	outputclust();

}//greedhclust


void preclust()
{
	//init the data structure
	int indexnum = 0;
	vector <pair<string, int> >::iterator tv = tvector.begin();
	for(; tv!=tvector.end(); tv++)
	{
		if(indexnum >= hbcount)
		{
			tlist.push_back(tv->first);
		}
		indexnum++;
	}
	int totallow = indexnum - hbcount;
	tvector.clear();
	//cout << tlist.size() << endl;
	//start the preclust process
	int tsize = tlist.size();
	//cout << totallow << "\t" << tsize <<endl;
	//	int k = 0;
	int indelnum = 0;
	int pre1num = 0;
	int max = 0 ;
	int min = 0;
	int sl = 0;
	int ql = 0;
	int Mismatch = 0;
	int diffr = 0;

	for(hash_map <string, vector<int> ,str_hash, compare_str>::iterator cur = otu.begin(); cur != otu.end(); cur++)
	{
		string otuname = cur->first;
		vector < int> otucell = cur->second;
		vector < int > tmlowseqid;

		for(vector <int>::iterator curr = otucell.begin(); curr != otucell.end(); curr++)
		{
			string sseq = id_seq[*curr];
			sl = sseq.size();
			for(list <string>::iterator lt = tlist.begin(); lt != tlist.end();)
			{
				string qseq = *lt;
				ql = qseq.size();
				Mismatch = 0;
				if(sl <= ql)
				{
					max = ql;
					min = sl;
				}
				else
				{
					max = sl;
					min = ql;
				}
				diffr = int(min * distcut); // diffr decided by the min len and the distcut
				Mismatch += (max - min);

				if((max - min) <= diffr)
				{
					for(int j=0; j< min; j++)
					{
						if(qseq[j] != sseq[j])
						{
							++Mismatch;
						}
						if(Mismatch > diffr)
						{
							Mismatch = max;
							break;
						}
					}
				}
				else
				{
					Mismatch = max;
				}

				if(Mismatch != 0 && Mismatch <= diffr)
				{
					++pre1num;
				}
				//indel check
				int indel = 0;
				if(sl == (ql+1) && qseq != sseq.substr(0,ql))
				{
					int j = 0;
					for(int i = 0; i < sl ; i++)
					{
						if(sseq[i] != qseq[j])
						{
							if(indel == 1)
							{
								indel = 2;
								break;
							}
							indel = 1;
							--j;
						}
						++j;
					}
				}

				if(ql == (sl+1) && sseq != qseq.substr(0,sl))
				{
					int j = 0;
					for(int i = 0; i < ql ; ++i)
					{
						if(qseq[i] != sseq[j])
						{
							if(indel == 1)
							{
								indel = 2;
								break;

							}
							indel = 1;
							--j;
						}
						++j;
					}
				}


				if(indel == 1)
				{
					Mismatch = 0;
					indelnum++;
				}

				if(Mismatch <= diffr)
				{
					int seqid = seq_id[qseq];
					tmlowseqid.push_back(seqid);
					id_otu[seqid] = otuname;
				//	tsize--;
					list <string>::iterator tmp = lt;
					lt++;
					tlist.erase(tmp);
				}
				else
				{
					lt++;
				}
			}
		}//for otucell
		//cout << indelnum <<  "\n" << endl;
		for(vector <int>::iterator cu = tmlowseqid.begin(); cu != tmlowseqid.end(); cu++)
		{
			otu[otuname].push_back(*cu);

		}
		//cout << endl;
		tmlowseqid.clear();
	}//for every otu
	//cout << tlist.size() << endl;
	LOG << "	Indel Number " << indelnum << "\tPreclust1\t" << pre1num <<  endl;
}//preclust

void lowabundantpreclust(map <int, vector<int> > & merge)
{
	list <string>::iterator cur = tlist.begin();
	list <string>::iterator ptr;
	int indelnum = 0;
	int preclust2 = 0;
	int max = 0;
	int min = 0;
	int ql = 0;
	int sl = 0;
	int Mismatch = 0;
	int diffr = 0;
	//cout << tlist.size() << endl;
	for(; cur != tlist.end(); cur++)
	{
		ptr=cur;
		string qseq = *cur;
		ql = qseq.length();
		int qid = seq_id[qseq];
		merge[qid].push_back(qid);
		ptr++;
		if(ptr != tlist.end())
		{
			for(; ptr != tlist.end();)
			{
				string sseq = *ptr;
				Mismatch = 0;
				sl = sseq.length();

				if(sl <= ql)
				{
					max = ql;
					min = sl;
				}
				else
				{
					max = sl;
					min = ql;
				}

				diffr = int(min * distcut); // diffr decided by the min len and the distcut

				Mismatch += (max - min);

				if((max - min) <= diffr)
				{
					for(int j=0; j< min; j++)
					{
						if(qseq[j] != sseq[j])
						{
							Mismatch++;
						}
						if(Mismatch > diffr)
						{
							Mismatch = max;
							break;
						}
					}
				}
				else
				{
					Mismatch = max;
				}

				if(Mismatch != 0 && Mismatch <= diffr)
				{
					++preclust2;
				}

				//indel check
				int indel = 0;
				if(sl == (ql+1) && qseq != sseq.substr(0,ql))
				{
					int j = 0;
					for(int i = 0; i < sl ; i++)
					{
						if(sseq[i] != qseq[j])
						{
							if(indel == 1)
							{
								indel = 2;
								break;
							}
							indel = 1;
							--j;
						}
						++j;
					}
				}

				if(ql == (sl+1) && sseq != qseq.substr(0,sl))
				{
					int j = 0;
					for(int i = 0; i < ql ; ++i)
					{
						if(qseq[i] != sseq[j])
						{
							if(indel == 1)
							{
								indel = 2;
								break;
							}
							indel = 1;
							--j;
						}
						++j;
					}
				}

				if(indel == 1)
				{
					Mismatch = 0;
					indelnum++;
				}

				if(Mismatch <= diffr)//there could be redundant for indel and mismatch for the diffr could > 0; if use Mismatch ==0 then the indel will include all the indel
				{
					int sid = seq_id[sseq];
					merge[qid].push_back(sid);
					//cout << qseq << "\n" << sseq << endl;
					list <string>::iterator temp = ptr;
					ptr++;
					tlist.erase(temp);
				}
				else
				{
					ptr++;
				}
			}//for
		}//if
	}
	LOG << "	low abundant indel	"  << indelnum  <<  "\tLow preclust2 number\t" <<  preclust2 <<  endl;
}


void processlowabseq()
{
	//Start to do greedy heuristic clustering for the left low abundant sequences in the tlist
	vector < pair<int, double> > ndvec;
	LOG << "	HB OTU number " << otu.size() << endl;
	LOG << "	Needle processed low abundant seq "<<  tlist.size() << endl;
	int cn = 0;
	time_t p1,p2;
	string hbo = "NO";
	double dis = -1;
	vector <string> lho;
	time(&p1);

	/*
	   const size_t arr_size=17;
	   int int_arr[arr_size]={29,31,32,34,35,36,37,39,40,45,55,59,75,86,144,170,205};;
	   vector<int> test(int_arr,int_arr+arr_size);
	 */
	//preclust low abundant sequences 
	map <int , vector<int> > lowpre;
	lowabundantpreclust(lowpre);

	LOG << "	low abundant preclust left sequences   " <<  tlist.size() << "\t" << lowpre.size() << endl;
	for(list <string>::iterator cur = tlist.begin(); cur != tlist.end(); cur++)
	{
		cn++;
		string qseq = *cur;
		int qid = seq_id[qseq];
		newotu = 1;

		do_kmerneedle_one2set(qseq, ndvec); //sub function from kmerneedle

		if(!ndvec.empty())
		{
			//merge the low otu first and if there only high abundant otu merge into the highest abundant one

			for(vector <pair<int, double> >::iterator cur = ndvec.begin(); cur != ndvec.end(); cur++)
			{
				string otuids = id_otu[cur->first];
				vector <string> temps;
				hbo = "NO";
				split_str(otuids,'_',temps);
				if( temps[0] == "OTUHB")
				{
					if(hbo !=  "NO")
					{
						if(cur->second < dis)
						{
							hbo = otuids;
						}
					}
					else
					{
						hbo = otuids;
						dis = cur->second;
					}
					//find the highest abundance otu
				}
				else if(temps[0] == "OTULB")
				{
					//cout << "herh"  << endl;
					if(find(lho.begin(),lho.end(),otuids) != lho.end())
					{
					}
					else
					{
						lho.push_back(otuids);
					}
				}
				else
				{
					cout << "wrong here " << temps[0]  << "\t" << otuids <<endl;
					exit(1);
				}
			}
			

			//we link the high abundant and low abundant, but we can not link high abundant and high abundant with low abundant
			if(!lho.empty() && hbo == "NO")
			{
				vector <string>::iterator lt = lho.begin();
				string  seed = *lt;
				lt++;
				if(lt != lho.end())
				{
					for(;lt != lho.end(); lt++)
					{
						vector <int> temp = otu[*lt];
						//otu[seed].insert(otu[seed].end(),temp.begin(),temp.end());
						for(vector <int>::iterator llt = temp.begin(); llt != temp.end(); llt++)
						{
							otu[seed].push_back(*llt);
							id_otu[*llt] = seed;
						}
						otu.erase(*lt);
					}
				}

				if(lowpre.find(qid) != lowpre.end())
				{
					vector <int> cellprelow = lowpre[qid];
					for(vector <int>::iterator cu = cellprelow.begin(); cu != cellprelow.end(); cu++)
					{
						otu[seed].push_back(*cu);
						id_otu[*cu] = seed;
					}
				}
				else
				{
					cout << "This low abundant is not seed" << endl;
					exit(1);
				}
				//otu[seed].push_back(qid);
				//id_otu[qid] = seed;

				newotu = 0;
			}
		
			if(!lho.empty() && hbo != "NO")
			{
				string seed = hbo;
				vector <string>::iterator lt = lho.begin();

				for(;lt != lho.end(); lt++)
				{
						vector <int> temp = otu[*lt];
						
						for(vector <int>::iterator llt = temp.begin(); llt != temp.end(); llt++)
						{
							otu[seed].push_back(*llt);
							id_otu[*llt] = seed;
						}
						otu.erase(*lt);
				}

				if(lowpre.find(qid) != lowpre.end())
				{
					vector <int> cellprelow = lowpre[qid];
					for(vector <int>::iterator cu = cellprelow.begin(); cu != cellprelow.end(); cu++)
					{
						otu[seed].push_back(*cu);
						id_otu[*cu] = seed;
					}
				}
				else
				{
					cout << "This low abundant is not seed" << endl;
					exit(1);
				}
				//otu[seed].push_back(qid);
				//id_otu[qid] = seed;

				newotu = 0;
			}

			if(lho.empty() && hbo !=  "NO")
			{
				if(lowpre.find(qid) != lowpre.end())
				{
					vector <int> cellprelow = lowpre[qid];
					for(vector <int>::iterator cu = cellprelow.begin(); cu != cellprelow.end(); cu++)
					{
						otu[hbo].push_back(*cu);
						id_otu[*cu] = hbo;
					}
				}
				else
				{
					cout << "This low abundant is not seed" << endl;
					exit(1);
				}
				//otu[hbo].push_back(qid);
				//id_otu[qid] = hbo;
				newotu = 0;
			}
			lho.clear();
		}

		if(newotu == 1)
		{
			otunumberindex++;
			string addotu;
			stringstream str(addotu);
			str << otunumberindex;
			addotu = str.str();
			addotu = "OTULB_" + addotu;

			if(lowpre.find(qid) != lowpre.end())
			{
				vector <int> cellprelow = lowpre[qid];
				for(vector <int>::iterator cu = cellprelow.begin(); cu != cellprelow.end(); cu++)
				{
					otu[addotu].push_back(*cu);
					id_otu[*cu] = addotu;
				}
			}
			else
			{
				cout << "This low abundant is not seed" << endl;
				exit(1);
			}
		}

		if(cn % 100 == 0)
		{
			time(&p2);
			LOG << "	"<<  cn << "\t" << difftime(p2,p1) << "\t" << ktime << "\t" << ntime  << "\t" << otu.size() << endl;
			p1 = p2;
		}

		ndvec.clear();
	}
}//processlowabseq


void outputclust()
{
	//Start to ouput the cluster result 
	fstream OTUR;
	string dist;
	stringstream str(dist);
	str << distcut;
	dist = str.str();
	string outr = oprefix + "_" + dist + "_otu." +  method + ".list";
	OTUR.open(outr.c_str(), ios::out);
	int valid = 0;

	hash_map <string, vector<int>, str_hash, compare_str >::iterator cur = otu.begin();
	for(;cur != otu.end(); cur++)
	{
		OTUR << cur->first << "\t" << cur->second.size() << "\t";
		vector <int>::iterator curr = cur->second.begin();
		OTUR << *curr ;
		valid++;
		curr++;
		if(curr != cur->second.end() )
		{
			for(;curr != cur->second.end(); curr++)
			{
				OTUR << ","<< *curr;
				valid++;
			}
		}
		OTUR << endl;
	}
	LOG  << "	seq num in otu " <<  valid  << "\t" << "otu number " <<  otu.size() << endl;
}//outputclust
