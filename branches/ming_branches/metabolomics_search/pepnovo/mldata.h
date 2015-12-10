#ifndef __MLDATA_H__
#define __MLDATA_H__

#ifndef __INCLUDES_H__
#include "includes.h"
#endif

struct IdxVal {
	IdxVal() : index(0), value(0.0) {}
	IdxVal(size_t i, value_t v) : index(i), value(v) {}

	size_t index;
	value_t value;
};

struct Sample {
	Sample() : label(-1), weight(1.0) {}

	void addPair(size_t index, value_t value) { pairs.push_back(IdxVal(index,value)); }
	void clear() { label=-1; weight =0; pairs.clear(); }
	void clearPairs() { pairs.clear(); }
	void print(ostream& os = cout) const;

	bool checkConsistency() const
	{
		cout << "PAIRS: " << pairs.size() << endl;
		size_t i;
		for (i=1; i<pairs.size(); i++)
		{
			cout << i << "\t" << pairs[i].index << "\t" << pairs[i].value << endl;
			if (pairs[i-1].index>=pairs[i].index)
				return false;
		}
		return true;
	}

	int		 label;
	weight_t weight;
	vector<IdxVal> pairs;
};


struct MlDataSet {
public:
	MlDataSet() : totalWeight_(0), numClasses_(0), numPairs_(0), maxFeatureIndex_(0),
				indStatsValid_(false), outStream_(0), buffer_(0), bufferSize_(0), bufferPos_(0) {};
	~MlDataSet();

	size_t				  getNumSamples()  const { return samples_.size(); }
	
	void  addSamples(const vector<Sample>& samples);
	void  addSamples(const vector<Sample>& samples, int label); // adds the samples and sets the labels to the new value

	void  addSample(const Sample& sample) { samples_.push_back(sample); indStatsValid_=false; }
	const vector<Sample>& getSamples()     const { return samples_; }
	const Sample&		  getSample(size_t i) const { return samples_[i]; }
	Sample&				  getSample(size_t i) { return samples_[i]; }
	
	weight_t			  getTotalWeight() const {  return totalWeight_; }

	size_t				  getNumClasess()  const { return numClasses_; }

	size_t				  getNumPairs()	   const  { return numPairs_; }

	const vector<weight_t>& getClassWeights() const { return classWeights_; }

	void	writeWholeDataFile(const char* dataFile);

	void	readDataFile(const char* dataFile, 
						 size_t aa = 99999999, 
						 bool randomSelection=false,
						 bool verbose=false);

	void	initializeOutputFile(const char* path);
	void	writeSampleToBuffer(const Sample& sample);
	void	closeOutputFile();

	void	printDatasetStatistics(ostream& os = cout) const;
	
	void	reserve(size_t numSamples) { samples_.reserve(numSamples); }

	void	tallyClassStatistics();
	
protected:

	vector<weight_t> classWeights_;
	vector<Sample> samples_;

	weight_t totalWeight_;
	size_t	 numClasses_;

	size_t numPairs_;
	size_t maxFeatureIndex_;

	bool   indStatsValid_;

	// used for on the fly writing of the dataset
	FILE*	outStream_;
	char*	buffer_;
	size_t	bufferSize_;
	size_t	bufferPos_;

	void	flushOutputBuffer();
	
};


#endif

