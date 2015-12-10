#ifndef __MLAUXFUN_H__
#define __MLAUXFUN_H__

#ifndef __INCLUDES_H__
#include "includes.h"
#endif



#ifndef __INCLLUDES_H__
void randomSeed(unsigned int init = 0);

unsigned int getRandomSeed();

/* Returns random uniform number */
double myRandom();


void error();
void error(const char *msg);
#endif

void chooseKFromN(size_t k, size_t n, vector<size_t>& idxs);


template<class T>
void mergeSortedVectors(const vector<T>& a, const vector<T>&b, vector<T>& m)
{
	typename vector<T>::const_iterator it_a = a.begin();
	typename vector<T>::const_iterator it_b = b.begin();

	m.clear();
	m.reserve(a.size()+b.size());

	while (it_a != a.end() && it_b != b.end())
	{
		if (*it_a < *it_b)
		{
			m.push_back(*it_a++);
			continue;
		}

		if(*it_a > *it_b)
		{
			m.push_back(*it_b++);
			continue;
		}

		m.push_back(*it_a++);
		it_b++;
	}

	while (it_a != a.end())
		m.push_back(*it_a++);

	while (it_b != b.end())
		m.push_back(*it_b++);
}

template<class T>
void intersecSortedVectors(const vector<T>& a, const vector<T>&b, vector<T>& intersec)
{
	intersec.clear();
	intersec.reserve(a.size()>b.size() ? b.size() : a.size());
	
	typename vector<T>::const_iterator it_a = a.begin();
	typename vector<T>::const_iterator it_b = b.begin();

	while (it_a != a.end() && it_b != b.end())
	{
		if (*it_a < *it_b)
		{	
			it_a++;
			continue;
		}

		if(*it_a > *it_b)
		{
			it_b++;
			continue;
		}

		intersec.push_back(*it_a++);
		it_b++;
	}
}



template<class T>
void createHistogram(vector<T>& vals, int num_bins, T min_val,
					  T max_val , ostream& os =cout)
{
	T bin_size = (max_val-min_val)/ (T)(num_bins);
	T one_over_bin = 1.0 / bin_size;
	int i;
	vector<int> counts;
	vector<float> percents;
	counts.resize(num_bins,0);

	for (i=0; i<vals.size(); i++)
	{
		if (vals[i]<min_val)
		{
			counts[0]++;
			continue;
		}

		int bin_idx = num_bins-1;
		if (vals[i]<max_val)
			bin_idx = (int)(one_over_bin*(vals[i]-min_val));
		
		counts[bin_idx]++;
	}

	T v = min_val;
	int tc=0;
	for (i=0; i<num_bins; i++)
	{
		os << setw(4) <<  setprecision(2) << left << v << " - ";
		v+= bin_size;
		os <<  setw(4) << left << v  << "  " <<  setw(6) << right << counts[i] << "  ";
		os << setw(6) << left << setprecision(4) << (float)counts[i]/(float)vals.size() << endl;
		tc+= counts[i];
	}

	os << "Total:       " << setw(6) << right << tc << "  " << setw(6) << setprecision(4) << left << (float)tc/(float)vals.size() << endl;
}


template<class T>
void createHistogram(const vector<T>& vals, const vector<T>& separator_vals, 
					  vector<size_t>& counts, ostream& os =cout)
{

	size_t i;
	
	vector<float> percents;
	counts.resize(separator_vals.size()+1,0);

	for (i=0; i<vals.size(); i++)
	{
		size_t j;
		for (j=0; j<separator_vals.size(); j++)
			if (vals[i]<separator_vals[j])
				break;
		
		counts[j]++;
	}
}



template<class T>
void calcMeanSd(const vector<T>& v, T *mean, T *sd)
{
	T m=0,var=0;
	
	if (v.size() == 0)
	{
		*mean=0;
		*sd=0;
		return;
	}

	if (v.size() == 1)
	{
		*mean=v[0];
		*sd=0;
	}

	size_t i;
	for (i=0; i<v.size(); i++)
		m+=v[i];

	m/=v.size();

	for (i=0; i<v.size(); i++)
		var+=(v[i]-m)*(v[i]-m);

	var /= v.size();

	*mean=m;
	*sd = sqrt(var);
}


template<class T>
void calcMeanSdFromCounts(const vector<T>& vals, const vector<int>& counts,
							  double *mean, double *sd)
{
	double m=0,var=0;
	
	if (vals.size() == 0)
	{
		*mean=0;
		*sd=0;
		return;
	}

	if (vals.size() != counts.size())
	{
		cout << "Error: values and counts should have same dimension!" << endl;
		exit(1);
	}

	size_t total_counts=0;
	size_t i;
	for (i=0; i<vals.size(); i++)
	{
		m+=vals[i]*counts[i];
		total_counts+=counts[i];
	}

	m/=(double)total_counts;

	for (i=0; i<vals.size(); i++)
		var+=(vals[i]-m)*(vals[i]-m)*counts[i];

	var /= (double)total_counts;

	*mean=m;
	*sd = sqrt(var);
}





#endif 
