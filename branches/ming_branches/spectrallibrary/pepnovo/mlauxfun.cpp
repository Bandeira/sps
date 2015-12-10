#include "mlauxfun.h"


#ifndef __INCLUDES_H__
static unsigned int ML_RAND_SEED;	

void randomSeed (unsigned int init)   {
	if (init != 0)
	{
		ML_RAND_SEED = init;
	}
	else
	{
		time_t ltime;
		unsigned int t=(unsigned int)time( &ltime );

		ML_RAND_SEED = t;
	}
}

unsigned int getRandomSeed() { return ML_RAND_SEED; }


/* Returns random uniform number */
double myRandom()  
{
  static unsigned int a = 1588635695, m = 4294967291U, q = 2, r = 1117695901;

   ML_RAND_SEED = a*(ML_RAND_SEED % q) - r*(ML_RAND_SEED / q);
   return ((double)ML_RAND_SEED / (double)m);
}


void error()
{
	cout << "Error!!!" << endl;
	exit(1);
}

void error(const char *msg)
{
	cout << "Error: " << msg << endl;
	exit(1);
}
#endif









struct ChoosePair {
	ChoosePair(size_t i,double d) : idx(i),val(d) {};
	bool operator < (const ChoosePair& other) const
	{
		return val<other.val;
	}
	size_t idx;
	double val;
};



// chooses k numbers from 0,...,n-1 (unique)
void chooseKFromN(size_t k, size_t n, vector<size_t>& idxs)
{

	
	if (k>n)
	{
		cout << "Error: choose " << k << " from " << n << " !" << endl;
		exit(1);
	}

	idxs.clear();
	idxs.resize(k);
	vector<ChoosePair> pairs;

	size_t i;
	for (i=0; i<n; i++)
		pairs.push_back(ChoosePair(i,myRandom()));
	
	sort(pairs.begin(),pairs.end());

	for (i=0; i<k; i++)
		idxs[i]=pairs[i].idx;

	sort(idxs.begin(),idxs.end());
}








