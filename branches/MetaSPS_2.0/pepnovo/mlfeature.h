#ifndef __MLFEATURE_H__
#define __MLFEATURE_H__

#ifndef __INCLUDES_H__
#include "includes.h"
#endif



class MlFeature {
	friend class MlFeatureSet;
public:
	MlFeature() : index_(MAX_INT), name_("Uninitialized"), flagConsiderForUnary_(true), 
		flagConsiderForConditional_(true) {}

	MlFeature(size_t i, const string& n, bool u, bool b) :
		index_(i), name_(n), flagConsiderForUnary_(u), flagConsiderForConditional_(b) {}

	bool operator< (const MlFeature& rhs) const
	{
		return (index_ < rhs.index_);
	}


	size_t getIndex() const { return index_; }
	string getName() const { return name_; }
	bool   getFlagConsiderUnary() const { return flagConsiderForUnary_; }
	bool   getFlagConsiderConditional() const { return flagConsiderForConditional_; }

	const vector<size_t>& getGroups() const { return groups_; }

private:
	
	size_t index_;
	string name_;
	
	// internal statistics that can be read from the file
	bool flagConsiderForUnary_;
	bool flagConsiderForConditional_;

	vector<size_t> groups_;
};



class MlFeatureSet {
public:
	void readFeatureFile(const char* path, bool verbose=false);

	void writeFeatureFile(const char* path, bool verbose=false) const;

	size_t	getNumFeatures() const { return features_.size(); }

	const vector<MlFeature>& getFeatures() const { return features_; }
	const MlFeature& getFeature(size_t i) const { return features_[i]; }

	void initializeWithDefaultFeatures(size_t numFeatures);

private:
	vector<MlFeature> features_;

	map<string, size_t> groupMapping_;
	vector<string>	   groupNames_;
	vector< vector<size_t> > groupMembers_;

};



/*
xxx changes that need to be made:
1. Feature is only an index and a name
2. Original data is expressed as a vector of feature values
3. The processing of the data is done through operations.
4. There are two types of operations:
	- unary, that operate on a single feature, convert type BOOL, LOG, EXP, ZVAL, SPLIT
5. Binary/Multinary that operate on more than one feature, e.g. AND OR NOT_AND, etc.
6. Multinary operations can also express the result as different functions: ONE, SUM, PROD, MAX, MIN, AVG, MEDIAN, FEATURE, etc.
7. Operations will be applied to the data in order. Each operation can only access features that have already
   been manipulated in a unary fashion (those features cannot recieve additional unary operations).
8. All additional features ar appended (in order of creation) to the end of the feature list
9. Most operations will add a conditional feature in addition to the real value (if it is not a bool itself)
*/




#endif

