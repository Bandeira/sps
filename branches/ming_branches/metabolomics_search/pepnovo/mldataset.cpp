#include "mldata.h"
#include "mlauxfun.h"
#include "auxfun.h"


MlDataSet::~MlDataSet()
{
	if (outStream_ && bufferPos_>0)
		closeOutputFile();

	if (buffer_)
		delete [] buffer_;
}


void MlDataSet::addSamples(const vector<Sample>& otherSamples)
{
	samples_.resize(samples_.size() + otherSamples.size());
	size_t i;
	for (i=0; i<otherSamples.size(); i++)
		samples_.push_back(otherSamples[i]);

	indStatsValid_ = false;
}


// adds the samples and sets the labels to the new value
void  MlDataSet::addSamples(const vector<Sample>& otherSamples, int label) 
{
	const size_t originalSize = samples_.size();
	samples_.resize(originalSize + otherSamples.size());
	size_t i;
	for (i=0; i<otherSamples.size(); i++)
	{
		samples_.push_back(otherSamples[i]);
		samples_[originalSize+i].label = label;
	}

	indStatsValid_ = false;
}



void MlDataSet::initializeOutputFile(const char* path)
{
	if (outStream_)
		closeOutputFile();

	outStream_=fopen(path,"w");
	if (! outStream_)
	{
		cout << "Error: couldn't open out stream: " << path << endl;
		exit(1);
	}

	if (buffer_)
	{
		bufferPos_=0;
		return;
	}


	bufferSize_=1048576;
	buffer_ = new char[bufferSize_];
	bufferPos_ =0;
}

void MlDataSet::writeSampleToBuffer(const Sample& sample)
{
	const size_t numPairs = sample.pairs.size();
	const size_t maxSampleSize = numPairs*sizeof(IdxVal) + 100;

	if (bufferSize_ - bufferPos_ < maxSampleSize)
		flushOutputBuffer();

	char* pos = buffer_ + bufferPos_;

	if (sample.weight != 1.0)
	{
		pos+=sprintf(pos,"%d $$$weight %.3f",sample.label, sample.weight);
	} else
		pos+=sprintf(pos,"%d",sample.label);

	for (size_t i=0; i<numPairs; i++)
	{
		const IdxVal& pair = sample.pairs[i];
		if (pair.value != 1.0)
		{
			pos+=sprintf(pos," F%lu 1.0",pair.index);
		}
		else
			pos+=sprintf(pos," F%lu %e",pair.index, pair.value);
	}
	*pos++='\n';
	bufferPos_ = (pos-buffer_);
}

void MlDataSet::closeOutputFile()
{
	if (bufferPos_>0)
		flushOutputBuffer();

	bufferPos_ = 0;
	fclose(outStream_);
}

void MlDataSet::flushOutputBuffer()
{
	if (! outStream_)
		return;

	size_t numBytesWritten = fwrite(buffer_,1,bufferPos_,outStream_);
	if (numBytesWritten != bufferPos_)
	{
		cout << "Error: could not empty data buffer to stream!" << endl;
		exit(1);
	}
	bufferPos_ = 0;
}


// performs a two pass scan:
// 1. first scan read number of lines, tokens, resolve feature names.
// 2. allocate space and read data.
void MlDataSet::readDataFile(const char *dataFile, 
						   size_t maxNumSamplesToRead,
						   bool randomSelection,
						   bool verbose)
{
	vector<size_t> numFeatures;

	ifstream stream(dataFile);
	if (! stream.is_open())
		error("couldn't open input file for reading: ",dataFile);

	char* buffer= new char[131072];
	
	// first pass
	while (! stream.eof())
	{
		stream.getline(buffer,131072);
		const size_t len = stream.gcount();
		if (len<=1 || buffer[0] == '#' || buffer[0] == '\n' || buffer[0] == '\r' )
			continue;

		// assume every feature name is prefixed with F (e.g., F0, F112, F1456, etc.)
		// count number of 'F'
		size_t n=0;
		size_t i;
		for (i=0; i<len; i++)
			if (buffer[i] == 'F')
				n++;
			
		numFeatures.push_back(n);
	}
	stream.close();

	// check if we need to be selective in choice of samples
	vector<bool> sampleFlags;
	sampleFlags.clear();
	if (numFeatures.size() > maxNumSamplesToRead)
	{
		sampleFlags.resize(numFeatures.size(),false);
		if (randomSelection)
		{
			if (verbose)
				cout << "Randomly reducing input file being read to " << 
						maxNumSamplesToRead << "/" << numFeatures.size() << endl;
			vector<size_t> idxs;
			chooseKFromN(maxNumSamplesToRead, numFeatures.size(), idxs);
			size_t i;
			for (i=0; i<idxs.size(); i++)
				sampleFlags[idxs[i]]=true;
		}
		else
		{
			if (verbose)
				cout << "Reducing input file being read to first " << 
						maxNumSamplesToRead << " samples (from a total of " << numFeatures.size() << ")" << endl;
			size_t i;
			for (i=0; i<maxNumSamplesToRead; i++)
				sampleFlags[i]=true;
		}
	}

	// second pass
	numPairs_ = 0;
	maxFeatureIndex_ = 0;

	stream.open(dataFile);
	if (! stream.is_open())
		error("couldn't open input file for reading: ",dataFile);

	samples_.resize(numFeatures.size());
	size_t samCounter=0;
	while (! stream.eof())
	{
		stream.getline(buffer,131072);
		const size_t len = stream.gcount();
		
		if (len<=1 || buffer[0] == '#' || buffer[0] == '\n' || buffer[0] == '\r' )
			continue;

		size_t i;
		for (i=0; i<len; i++)
			if (buffer[i]=='F')
				buffer[i]=' ';

		if (sampleFlags.size() > 0 && ! sampleFlags[samCounter])
		{
			samCounter++;
			continue;
		}

		const size_t n = numFeatures[samCounter];
		Sample& sample = samples_[samCounter];
		samCounter++;

		if (n == 0)
			continue;

		numPairs_ += n;
		sample.pairs.resize(n);
		istringstream iss(buffer);

		iss >> sample.label;
		if (iss.gcount()<1)
			error("error parsing input file, sample ",samCounter);

		size_t start=0;
		size_t index;
		iss >> index;
		if (iss.gcount()<1) // is this a weight string
		{	
			string wStr;
			iss >> wStr;

#ifdef MY_DEBUG_CHECKS
			if (wStr != "$$$weight")
				error("error parsing input file at first feature ($$$weight), sample ",samCounter);
#endif

		}
		else
		{
			sample.pairs[0].index = index;
			iss >> sample.pairs[0].value;
#ifdef MY_DEBUG_CHECKS
			if (iss.gcount()<1)
				error("error parsing input file at first feature, sample ",samCounter);
#endif
			start=1;
		}

		for (i=start; i<n; i++)
		{
			iss >> sample.pairs[i].index >> sample.pairs[i].value;
#ifdef MY_DEBUG_CHECKS
			if (iss.gcount()<1)
				error("error parsing input file, sample ",samCounter, " feature ",i);
#endif
		}

		if (sample.pairs[n-1].index > maxFeatureIndex_)
			maxFeatureIndex_ = sample.pairs[n-1].index;
	}

	delete [] buffer;

	tallyClassStatistics();
	if (verbose)
		printDatasetStatistics();
}


void MlDataSet::tallyClassStatistics()
{
	numPairs_=0;
	maxFeatureIndex_=0;

	vector<double> weights;
	double total=0;
	size_t i;
	for (i=0; i<samples_.size(); i++)
	{
		size_t c = samples_[i].label;
		double w = static_cast<double>(samples_[i].weight);

		total+=w;
		if (c>= weights.size())
			weights.resize(c+1,0);
		weights[c]+=w;

		const vector<IdxVal>& pairs = samples_[i].pairs;
		numPairs_ += pairs.size();
		if (pairs[pairs.size()-1].index > maxFeatureIndex_)
			maxFeatureIndex_ = pairs[pairs.size()-1].index;
	}

	numClasses_ = weights.size();
	totalWeight_ = static_cast<weight_t>(total);
	classWeights_.resize(numClasses_);
	for (i=0; i<numClasses_; i++)
		classWeights_[i]=static_cast<weight_t>(weights[i]);

	indStatsValid_ = true;
}


void	MlDataSet::writeWholeDataFile(const char* dataFile)
{
	initializeOutputFile(dataFile);
	
	size_t i;
	for (i=0; i<samples_.size(); i++)
		writeSampleToBuffer(samples_[i]);

	closeOutputFile();
}


void	MlDataSet::printDatasetStatistics(ostream& os) const
{
	os << "DataSet contains:" << endl;
	os << samples_.size() << "\tsamples" << endl;
	os << maxFeatureIndex_+1 << "\tfeatures (maximal index is " << maxFeatureIndex_ << ")" << endl;
	os << numPairs_ << " index-value pairs, avg of " << setprecision(1) << numPairs_ / static_cast<float>(samples_.size()) <<
		" per sample." << endl;
	os << "Samples belong to " << numClasses_ << " classes:" << endl;
	size_t i;
	for (i=0; i<numClasses_; i++)
		cout << "Class " << i <<"\t weight " << classWeights_[i] << " (" << classWeights_[i]/totalWeight_ << ")" << endl;
}


void Sample::print(ostream& os) const
{
	os << "LABEL  " << label << endl;
	os << "WEIGHT " << weight << endl;
	os << "INDEX-VALUE PAIRS: "<< pairs.size() << endl;
	int i;
	for (i=0; i<pairs.size(); i++)
		os << i << "\t" << pairs[i].index << "\t" << pairs[i].value << endl;
}

