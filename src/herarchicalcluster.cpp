#include "herarchicalcluster.h"

//////////
typedef int ui;

map <int, int> tempseqmap;
vector < int > tempvec;

/*************************************************************************/
//data store structure 
//basic cell

struct Pcell{
	ui row;
	ui col;
	float dist;
	Pcell ** vectorMap;
	Pcell() : row(0),col(0),dist(0),vectorMap(NULL) {};
	Pcell(ui r, ui c, float d):row(r), col(c), dist(d), vectorMap(NULL) {};
};

typedef list<Pcell>::iterator MatData; //iterator
typedef vector<MatData> MatVec;
typedef vector <pair<int, double> > MatP;

list <Pcell> matrix;
ui numNodes = 0;
float smallDist = 1.0;
ui smallRow = 0;
ui smallCol = 0;
ui numberofdist = 0;
Pcell *smallCell; //the cell with the smallDist
vector <Pcell*> mins; //

vector <MatVec> seqVec; //store all the dist for a specific sequence
MatVec rowCells;
MatVec colCells;
ui nrowCells;
ui ncolCells;

map <int, vector<int> > seqname; //store the name vector
vector <int> binSeq; //store then number of seqs for every cell

void addCell(Pcell va)
{
	matrix.push_back(va);
	numNodes++;
	if(va.dist < smallDist)
	{
		smallDist = va.dist;
	}
}

MatData rmCell(MatData data)
{
	if(data->vectorMap != NULL)
	{
		*(data->vectorMap) = NULL;
		data->vectorMap = NULL;
	}
	data = matrix.erase(data);
	numNodes--;
	return(data);
}

void removeCell(const MatData& cell, int vrow, int vcol)
{
	ui drow = cell->row;
	ui dcol = cell->col;

	if (((vrow >=0) && (drow != smallRow)) ||
			((vcol >=0) && (dcol != smallCol))) {
		ui dtemp = drow;
		drow = dcol;
		dcol = dtemp;
	}

	ui crow;
	ui  ccol;
	int nCells;

	if (vrow < 0) {
		nCells = seqVec[drow].size();
		for (vrow=0; vrow<nCells;vrow++) {
			crow = seqVec[drow][vrow]->row;
			ccol = seqVec[drow][vrow]->col;
			if (((crow == drow) && (ccol == dcol)) ||
					((ccol == drow) && (crow == dcol))) {
				break;
			}
		}
	}

	seqVec[drow].erase(seqVec[drow].begin()+vrow);
	if (vcol < 0) {
		nCells = seqVec[dcol].size();
		for (vcol=0; vcol<nCells;vcol++) {
			crow = seqVec[dcol][vcol]->row;
			ccol = seqVec[dcol][vcol]->col;
			if (((crow == drow) && (ccol == dcol)) ||
					((ccol == drow) && (crow == dcol))) {
				break;
			}
		}
	}

	seqVec[dcol].erase(seqVec[dcol].begin()+vcol);

	rmCell(cell);

}//removeCell


Pcell * getsmallestCell()
{
	while((!mins.empty()) && (mins.back() == NULL) )
	{
		mins.pop_back();
	}

	if(mins.empty())
	{
		mins.clear();
		smallDist = matrix.begin()->dist;

		for(MatData currentCell = matrix.begin(); currentCell != matrix.end(); currentCell++)
		{
			float tmdist = currentCell->dist;
			if(tmdist < smallDist)
			{
				mins.clear();
				smallDist = tmdist;
				mins.push_back(&*currentCell);
			}
			else if(tmdist == smallDist)
			{
				mins.push_back(&*currentCell);
			}
		}

		random_shuffle(mins.begin(), mins.end());  //randomize the order of the iterators in t    he mins vector

		for(int i=0;i<mins.size();i++){
			mins[i]->vectorMap = &mins[i];
			//	cout << mins[i]->dist << endl;
		}
		//	cout << endl;
	}

	smallCell = mins.back();
	mins.pop_back();

	return smallCell;
}

void cluster_al()
{
	//int myva = 0;
	for(map <int, vector< pair<int, double> > >::iterator cur = hbneedle.begin(); cur != hbneedle.end(); ++cur)
	{
		int ida = cur->first;
		if(hbseed.find(ida) != hbseed.end())
		{
			int ida2 = tempseqmap[ida];
			MatP celldist = cur->second;
			for(MatP::iterator curr = celldist.begin(); curr != celldist.end(); ++curr)
			{
				if(hbseed.find(curr->first) != hbseed.end())
				{
					int idb2 = tempseqmap[curr->first];
					//myva++;
					if(ida2 > idb2)
					{
						Pcell value(ida2, idb2,curr->second);
						addCell(value);
					}
					else
					{
						Pcell value(idb2, ida2, curr->second);
						addCell(value);
					}
				}
			}//for
		}//if
	}//for
	//cout << "herh 1 "  <<  myva <<  endl;
	//init seqVec
	LOG << "	hclust total dist number need process " <<  matrix.size() <<  "\t NumberofNodes\t" << numNodes << endl;
	seqVec = vector <MatVec>(hbseed.size());
	for(MatData currentCell = matrix.begin(); currentCell != matrix.end(); ++currentCell)
	{
		//cout << "herh 111"  << endl;
		seqVec[currentCell->row].push_back(currentCell);
		seqVec[currentCell->col].push_back(currentCell);
	}//for

	float oldrnddist = 0.0000;
	float rnddist = 0.0000;

	double cutoff = distcut;
	for(int i = 0; i < hbseed.size(); ++i)
	{
		binSeq.push_back(1);
	}

	while( (smallDist <= cutoff) && (numNodes > 0) )
	{
		smallCell = getsmallestCell();
		smallRow = smallCell->row;
		smallCol = smallCell->col;
		smallDist = smallCell->dist;

		//if(smallRow >= numberofdist || smallCol >= numberofdist){ exit(1);}
		if(smallDist > cutoff){ break;}

		rowCells = seqVec[smallRow];
		colCells = seqVec[smallCol];
		nrowCells = rowCells.size();
		ncolCells = colCells.size();
		ui binRow = binSeq[smallRow];
		ui binCol = binSeq[smallCol];

		//cout << smallRow << "\t" << smallCol << "\t" << smallDist << "\t" << binRow << "\t" << binCol << endl;

		vector <int> foundCol(ncolCells,0);
		int search;
		//bool change = false;

		for(int i = nrowCells-1; i>=0; --i)
		{
			if(! (rowCells[i]->row == smallRow && rowCells[i]->col == smallCol))
			{
				
				if(rowCells[i]->row == smallRow)
				{
					search = rowCells[i]->col;
				}
				else
				{
					search = rowCells[i]->row;
				}
				bool merge = false;
				for(int j=0; j<ncolCells; ++j)
				{
					if(!(colCells[j]->row == smallRow && colCells[j]->col == smallCol))
					{
						if(colCells[j]->row == search || colCells[j]->col == search)
						{
							foundCol[j] = 1;
							merge = true;

							//change = updateDistence(colCells[j],rowCells[i],method);

							ui bintotal = binRow + binCol;
							colCells[j]->dist = (binCol * colCells[j]->dist + binRow * rowCells[i]->dist) / bintotal;

							//if(change)
							//{
							if(colCells[j]->vectorMap != NULL)
							{
								*(colCells[j]->vectorMap) = NULL;
								colCells[j]->vectorMap = NULL;
							}
							//}
							break;
						}

					}
				}//for
				if((!merge) && (method == "al"))
				{
					if(cutoff > rowCells[i]->dist)
					{
						cutoff = rowCells[i]->dist;
					}
				}
			
				removeCell(rowCells[i],i,-1);
			}//if
		}//for
	
		//update the binSeq
		binSeq[smallRow] = 0;
		binSeq[smallCol] = binRow + binCol;

		//cout << binCol << "\t" << binSeq[smallCol] << endl;

		seqname[smallCol].insert(seqname[smallCol].end(), seqname[smallRow].begin(), seqname[smallRow].end());
		seqname[smallRow].clear();

		for(int i=ncolCells-1; i>=0; --i)
		{
			if(foundCol[i] == 0)
			{	
				if(method == "al")
				{
					if(!((colCells[i]->row == smallRow) && (colCells[i]->col == smallCol)))
					{
						if(cutoff > colCells[i]->dist)
						{
							cutoff = colCells[i]->dist;
						}

					}
				}
				removeCell(colCells[i], -1 , i);
			}
		}	

		rowCells.clear();
		colCells.clear();
		nrowCells = 0;
		ncolCells = 0;
	}//while

	//output the HBCLUSER and init the HBOTU
	if(cutoff == distcut)
	{
		fstream HBCLUSTER;
		string hbclust = oprefix + ".hb.al.cluster";
		HBCLUSTER.open(hbclust.c_str(), ios::out);
		if(HBCLUSTER == NULL)
		{
			exit(1);
		}

		//int sta = 0;
		HBCLUSTER << distcut << "\t|";
		map <int,vector <int> >::iterator cu = seqname.begin();
		vector <int> tmv;
		for( ;cu != seqname.end(); ++cu)
		{	
			if( !(cu->second).empty())
			{
				tmv = cu->second;
				++otunumberindex;
				string addotu;
				stringstream str(addotu);
				str << otunumberindex;
				addotu = str.str();
				addotu = "OTUHB_" + addotu;
				//cout << addotu << endl;

				//otu[addotu].assign(tmv.begin(), tmv.end());
				for(vector <int>::iterator cur = tmv.begin(); cur != tmv.end(); ++cur)
				{
							//++sta;
							otu[addotu].push_back(*cur);
 							HBCLUSTER << *cur << " ";
							id_otu[*cur] = addotu;
							//cout << *cur << endl;
				}
				HBCLUSTER << "|";
				tmv.clear();
			} 
		}
		//cout << "number of seqname\t"<<  sta << endl;
	}
	else
	{
		cout << "The matrix is not compelte for the average linkage clustering" << endl; exit(1);
	}
	
}
////The above is for mothur h cluster
/*****************************************************************************************************************/
void hc_excute()
{
	//caculate the needleman dist and store into the hbneedle map
	time_t t1;time(&t1);
	LOG << "	[1] do kmer and needle caculation" << endl;
	do_kmerneedle_set(hbneedle);	
	LOG << "	cost time "<< difftime(time(NULL), t1) << " s" << endl; 
	//if(!hbneedle.empty())
	//{
		//denoise by 0.02
	LOG << "	[2] denoise by a mismatch "<< endl;
		time_t t2; time(&t2);
	denoise();
		LOG << "	cost time "<< difftime(time(NULL), t2) << " s" << endl; 
		//do hcluster for the denoised seed, first select the dist of the seeds
		LOG << "	[3] do hclust for the high abundant sequences "<< endl;
		time_t t3; time(&t3);
	hcluster();
		LOG << "	cost time "<< difftime(time(NULL), t3) << " s" << endl; 
	//}
}

void denoise()
{
	//do denoise 
	vector <int> mark(hbcount, 0);
	for(int j = 0; j < hbcount; ++j)
	{
		//if(j==431){ cout << mark[j] <<endl;}
		if(mark[j] == 0)
		{
			if(hbneedle.find(j) != hbneedle.end())
			{
				hbseed[j].push_back(j);
				vector < pair<int, double> > vp = hbneedle[j]; 
				for(vector <pair<int, double> >::iterator cr = vp.begin(); cr != vp.end(); ++cr)
				{
					if(cr->second <= 0.02)
					{
						if(mark[cr->first] == 0)
						{
							mark[cr->first] = 1;
							hbseed[j].push_back(cr->first);
						}
					}
				}
			}
			else
			{
				//generate new seed
				hbseed[j].push_back(j);
				//cout << "got here? " << endl;
				mark[j] = 1;
			}
		}//if
	}//for
}

/**********************************************************************************************/
//we init the high abundant otu in this step
void hcluster()
{

	fstream DENOISE;
	string denois = oprefix + "_0.02_denoise.list";
	DENOISE.open(denois.c_str(), ios::out);
	if(DENOISE == NULL)
	{
		LOG << "Denoise file can not write " << endl; exit(1);
	}
	for(map<int, vector<int> >::iterator cur = hbseed.begin(); cur != hbseed.end(); ++cur)
	{
		DENOISE << cur->first << "\t";
		vector <int>::iterator curr = cur->second.begin();
		DENOISE << *curr;
		++curr;
		if(curr != cur->second.end())
		{
			for(; curr != cur->second.end(); ++curr)
			{
				DENOISE << "," << *curr;
			}
		}
		DENOISE << endl;
	}//for

	LOG  <<  "	Total sequences number " << tvector.size() <<  "\n	High abundant sequence " << hbcount  << "\n	High abundant seeds " << hbseed.size() << endl;
	//find hbseed distence from the hbneedle
	//transform it into another data structure

	int pos=0;
	map<int, vector<int> >::iterator cur = hbseed.begin();
	for(;cur != hbseed.end(); ++cur)
	{
		tempseqmap[cur->first] = pos;
		tempvec.push_back(cur->first);
		//for mothur cluster
		//sotre the seqmap
		seqname[pos].assign((cur->second).begin(), (cur->second).end());
		pos++;
	}

	if(method == "cl")
	{
		map < double, multimap < int, int > > needledist;
		cur = hbseed.begin();

		for(;cur != hbseed.end(); cur++)
		{
			if(hbneedle.find(cur->first) != hbneedle.end())
			{
				for(vector <pair<int, double> >::iterator  lt = hbneedle[cur->first].begin(); lt != hbneedle[cur->first].end(); lt++)
				{
					//cout << "first 1" << endl;
					if(hbseed.find(lt->first) != hbseed.end())
					{
						int index1 = tempseqmap[cur->first];
						int index2 = tempseqmap[lt->first];
						if(index1 > index2)
						{
							int tm = index1;
							index1 = index2;
							index2 = tm;
						}
						if(lt->second <= distcut)
						{
							needledist[lt->second].insert(make_pair(index1,index2));
							//cout <<  index1 <<  "\t"  << index2 << "\t"  << lt->second << endl;
						}
					}
				}
			}
		}

		int seqnum = tempseqmap.size();
		vector < int > num_seqs(seqnum, 1);
		vector < int > parent(seqnum, -1);

		map < pair<int, int>, int > seqpair;
		map < pair<int, int>, double > pairdist;
		map < double, multimap < int, int > >::iterator pos1 = needledist.begin();
		map < int, int >::iterator pos2;
		int seqnum1 = seqnum ;

		double dist1;
		int l = 0;

		dist1 = pos1->first;
		LOG << "	hclust total dist number need process " <<  needledist.size() << endl;
		while(dist1<=distcut && pos1 != needledist.end())
		{
			l++;
			dist1 = pos1->first;
			//if(l % 100 == 0)
			//{
			//	LOG << "	" << l << endl;
			//}
			////#pragma omp parallel for if(i>20)num_threads(ompnums)
			//printf("----- %d\n\n",(pos1->second).size());
			for (pos2 = (pos1->second).begin(); pos2 != (pos1->second).end(); ++pos2)
			{
				int i = pos2->first;
				int j = pos2->second;
				//printf("%d %d ----- %f\n\n",i,j,dist1);
				while (parent[i] != -1)
				{
					i = parent[i];
				}
				while (parent[j] != -1)
				{
					j = parent[j];
				}
				if (i>j)
				{
					int p1 = i;
					i = j;
					j = p1;
				}
				if (seqpair.find(make_pair(i,j)) != seqpair.end())
				{
					seqpair[make_pair(i,j)]++;
				}
				else
				{
					seqpair.insert(make_pair(make_pair(i,j),1));
				}
				if (pairdist.find(make_pair(i,j)) != pairdist.end())
				{
					if (pairdist[make_pair(i,j)] < dist1)
					{
						pairdist[make_pair(i,j)] = dist1;
					}
				}
				else
				{
					pairdist.insert( make_pair(make_pair(i,j),dist1) );
				}
				//printf("%d %d = %d\n\n",i,j,seqpair[make_pair(i,j)]);
				if (seqpair[make_pair(i,j)] == num_seqs[i]*num_seqs[j])
				{

					num_seqs.push_back(num_seqs[i]+num_seqs[j]);				
					parent.push_back(-1);
					//printf("%d %d -> %d\n\n",i,j,seqnum1);
					parent[i] = seqnum1;
					parent[j] = seqnum1;
					seqpair.erase(make_pair(i,j));

					map < pair<int, int>, int >::iterator pos3;
					for (pos3 = seqpair.begin(); pos3 != seqpair.end();)
					{
						int p1 = pos3->first.first;
						int p2 = pos3->first.second;
						int p3 = pos3->second;
						pos3++;
						if (p1 == i || p1 == j)
						{
							if (seqnum1 > p2)
							{
								seqpair[make_pair(p2,seqnum1)] += p3;
							}
							else
							{
								seqpair[make_pair(seqnum1,p2)] += p3;
							}
							seqpair.erase(make_pair(p1,p2));
						}
						else if (p2 == i || p2 ==j)
						{
							if (seqnum1 > p1)
							{
								seqpair[make_pair(p1,seqnum1)] += p3;
							}
							else
							{
								seqpair[make_pair(seqnum1,p1)] += p3;
							}
							seqpair.erase(make_pair(p1,p2));
						}
					}
					seqnum1++;
					l++;
				}
				/*
				   map < pair<int, int>, int >::iterator pos3;
				   for (pos3 = seqpair.begin(); pos3 != seqpair.end();pos3++)
				   {
				   int p1 = pos3->first.first;
				   int p2 = pos3->first.second;
				   int p3 = pos3->second;
				   printf("%d %d *** %d\n",p1,p2,p3);
				   }
				   */
			}
			pos1++;
		}


		map <int, vector < int > > hcluster;
		for (int i=parent.size()-1; i>=0; --i)
		{
			while ((parent[i] != -1) && (parent[parent[i]] != -1))
			{
				parent[i] = parent[parent[i]];
			}
		}
		for (int i=0; i<seqnum; ++i)
		{
			if (hcluster.find(parent[i]) != hcluster.end())
			{
				hcluster[parent[i]].push_back(i);
			}
			else
			{
				vector < int > seq;
				seq.push_back(i);
				hcluster.insert(make_pair(parent[i], seq));
			}	
		}

		fstream HBCLUSTER;
		string hbclust = oprefix + ".hb.cl.cluster";
		HBCLUSTER.open(hbclust.c_str(), ios::out);
		if(HBCLUSTER == NULL)
		{
			exit(1);
		}

		map <int, vector < int > >::iterator pos3;
		int otu_nums = 0;
		HBCLUSTER << distcut << "\t|" ;

		for (pos3 = hcluster.begin(); pos3 != hcluster.end(); ++pos3)
		{
			if(pos3->first == -1)
			{
				for (int i=0; i<pos3->second.size(); ++i)
				{
					otunumberindex++;
					string addotu;
					stringstream str(addotu);
					str << otunumberindex;
					addotu = str.str();
					addotu = "OTUHB_" + addotu;

					if(hbseed.find(tempvec[(pos3->second)[i]]) != hbseed.end())
					{
						vector <int>::iterator tm = hbseed[tempvec[(pos3->second)[i]]].begin();
						for (; tm != hbseed[tempvec[(pos3->second)[i]]].end(); tm++ )
						{
							HBCLUSTER << *tm   << " ";
							otu[addotu].push_back(*tm);
							id_otu[*tm] = addotu;
						}
						HBCLUSTER <<  "|";
					}
					else
					{
						cout << "wrong" << endl;
						exit(1);
					}
					otu_nums++;
				}
			}
			else
			{
				otunumberindex++;
				string addotu;
				stringstream str(addotu);
				str << otunumberindex;
				addotu = str.str();
				addotu = "OTUHB_" + addotu;

				for (int i=0; i<pos3->second.size(); ++i)
				{
					if(hbseed.find(tempvec[(pos3->second)[i]]) != hbseed.end())
					{
						vector <int>::iterator tm = hbseed[tempvec[(pos3->second)[i]]].begin();
						for (; tm != hbseed[tempvec[(pos3->second)[i]]].end(); tm++ )
						{
							HBCLUSTER << *tm   << " ";
							otu[addotu].push_back(*tm);
							id_otu[*tm] = addotu;
						}
					}
				}
				otu_nums++;
				HBCLUSTER << "|";
			}
		}
		HBCLUSTER << endl;
		//cout << "There" <<  otu_nums <<  endl;
	}
	else if(method == "sl")
	{

		map <double, vector< pair<int,int> > > needledist;

		for(map <int,vector <pair<int,double> > >::iterator cur = hbneedle.begin(); cur != hbneedle.end(); cur++)
		{
			if(hbseed.find(cur->first) != hbseed.end())
			{
				int index1 = tempseqmap[cur->first];
				for(vector <pair<int, double> >::iterator lt = (cur->second).begin(); lt != (cur->second).end(); lt++)
				{
					if(hbseed.find(lt->first) != hbseed.end())
					{
						int index2 = tempseqmap[lt->first];
						if(index1 > index2)
						{
							int tm = index1;
							index1 = index2;
							index2 = tm;
						}
						if(lt->second <= distcut)
						{
							needledist[lt->second].push_back(make_pair(index1, index2));
						}
					}
				}
			}
		}

		int seqnum = tempseqmap.size();
		//cout << "seqnum\t"<< seqnum << endl;
		//vector < int > parent(seqnum,-1);
		vector < int > parent(seqnum,-1);
		int seqnum1 = seqnum;

		LOG << "	hclust total dist number need process " <<  needledist.size() << endl;
		//cout << "here 1 "  << needledist.size()  <<endl;
	if(! needledist.empty())
	{
		map < double, vector< pair<int, int> > >::iterator pos1= needledist.begin();
		vector< pair<int, int> >::iterator pos2;

		double dist1 = 0;
		//double mydist1 = 0.00;
		int l = 0;
		//cout << "\nSingle linkage clustering...\n" ;	
		//	for (pos1 = needledist.begin(); pos1 != needledist.end();)
		//	{
		dist1 = pos1->first;
		//2011-04-12printf("\r n=%-10d ",l);
		//fflush(stdout);	
		//int flag1 = 0;
		
		while ((dist1 <= distcut) && (pos1 != needledist.end()))
		{
			//cout <<  "\n"<< dist1 << endl;
			//#pragma omp parallel for if(i>20)num_threads(ompnums)
			for (pos2 = (pos1->second).begin(); pos2 != (pos1->second).end(); ++pos2)
			{
				int i = pos2->first;
				int j = pos2->second;
				l++;
			//2011-04-12	printf("\r n=%-10d ",l);
			//2011-04-12	fflush(stdout);	
				while (parent[i] != -1)
				{
					i = parent[i];
				}
				while (parent[j] != -1)
				{
					j = parent[j];
				}
				if (i != j)
				{
					if (i>j)
					{
						int p1 = i;
						i = j;
						j = p1;
					}
					//flag1 = 1;
					parent[i]=seqnum1;
					parent[j]=seqnum1;
					parent.push_back(-1);
					seqnum1++;
				}
			}
			pos1++;
			dist1 = pos1->first;

		}
		//	mydist1 += 0.01;
		//}
		/*
		double dist1 = pos1->first;
		int flag1 = 0;

		
		while ((dist1 <= distcut) && (pos1 != needledist.end()))
		{
			for (pos2 = (pos1->second).begin(); pos2 != (pos1->second).end(); ++pos2)
			{
				int i = pos2->first;
				int j = pos2->second;

				while (parent[i] != -1)
				{
					i = parent[i];
				}
				while (parent[j] != -1)
				{
					j = parent[j];
				}
				if (i != j)
				{
					if (i>j)
					{
						int p1 = i;
						i = j;
						j = p1;
					}
					flag1 = 1;
					parent[i]=seqnum1;
					parent[j]=seqnum1;
					parent.push_back(-1);
					seqnum1++;
				}
			}
			pos1++;
			dist1 = pos1->first;

		}//while
		*/
		}

		if (1)
		{

			map <int, vector < int > > SL_cluster;
			int otu_num = 0;
			for (int i=parent.size()-1; i>=0; --i)
			{
				while ((parent[i] != -1) && (parent[parent[i]] != -1))
				{
					parent[i] = parent[parent[i]];
				}
				if (parent[i] == -1)
				{
					otu_num++;
				}
			}

			for (int i=0; i< seqnum; ++i)
			{
				if (SL_cluster.find(parent[i]) != SL_cluster.end())
				{
						SL_cluster[parent[i]].push_back(i);
				}
				else
				{
					vector < int > seq;
					seq.push_back(i);
					SL_cluster.insert(make_pair(parent[i], seq));
				}	
			}

			fstream HBCLUSTER;
			string hbclust = oprefix + ".hb.sl.cluster";
			HBCLUSTER.open(hbclust.c_str(), ios::out);
			if(HBCLUSTER == NULL)
			{
				exit(1);
			}
			//fprintf(fid1, "%-.2f %-4d ", mydist1, otu_num);
			int otu_nums = 0;

			
			HBCLUSTER << distcut << "\t|" ;
			map <int, vector < int > >::iterator pos3;
			for (pos3 = SL_cluster.begin(); pos3 != SL_cluster.end(); ++pos3)
			{
				if (pos3->first == -1)
				{
					for (int i=0; i<pos3->second.size(); ++i)
					{
						
						otunumberindex++;
						string addotu;
						stringstream str(addotu);
						str << otunumberindex;
						addotu = str.str();
						addotu = "OTUHB_" + addotu;
/*
						HBCLUSTER << (pos3->second)[i] << " ";
						otu[addotu].push_back((pos3->second)[i]);
						id_otu[(pos3->second)[i]] = addotu;
						HBCLUSTER <<  "|";
*/

						if(hbseed.find(tempvec[(pos3->second)[i]]) != hbseed.end())
						{
							vector <int>::iterator tm = hbseed[tempvec[(pos3->second)[i]]].begin();
							for (; tm != hbseed[tempvec[(pos3->second)[i]]].end(); tm++ )
							{
								HBCLUSTER << *tm   << " ";
								otu[addotu].push_back(*tm);
								id_otu[*tm] = addotu;
							}
							HBCLUSTER <<  "|";
						}
						else
						{
							LOG << "	wrong" << endl;
							exit(1);
						}
						otu_nums++;

					}
				}
				else
				{

					otunumberindex++;
					string addotu;
					stringstream str(addotu);
					str << otunumberindex;
					addotu = str.str();
					addotu = "OTUHB_" + addotu;

					for (int i=0; i<pos3->second.size(); ++i)
					{
						
						if(hbseed.find(tempvec[(pos3->second)[i]]) != hbseed.end())
						{
							vector <int>::iterator tm = hbseed[tempvec[(pos3->second)[i]]].begin();
							for (; tm != hbseed[tempvec[(pos3->second)[i]]].end(); tm++ )
							{
								HBCLUSTER << *tm   << " ";
								otu[addotu].push_back(*tm);
								id_otu[*tm] = addotu;
							}
						}
						/*
						
						HBCLUSTER << (pos3->second)[i] << " ";
						otu[addotu].push_back((pos3->second)[i]);
						id_otu[(pos3->second)[i]] = addotu;
						*/

					}
					otu_nums++;
					HBCLUSTER << "|";

				}
				//printf(" dist=%-10.2f OTUnums=%-10d\n", mydist1, otu_num);
				//fprintf(fid1, "\n");
				if (otu_num == 1)
				{
					break;
				}
			}//if
			HBCLUSTER << endl;
		}
	}
	else if(method == "al")
	{
		cluster_al();
	}
	else
	{
		LOG  << "Your method not in our programe please input fn|av|nn" << endl;
		exit(1);
	}
	//delete some old data
	hbneedle.clear();
	hbseed.clear();
}//hc_excute
