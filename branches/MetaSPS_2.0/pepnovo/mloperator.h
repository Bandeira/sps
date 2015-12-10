#ifndef __MLOPERATOR_H__
#define __MLOPERATOR_H__

#include "mlfeature.h"
#include "mldata.h"



typedef enum UnaryOperatorType  { // unary operators, need only one feature to work with
							 UOT_BOOL, // output f(x)=1 if the feature was found in the sample
							 UOT_LOG,  // output f(x)=log(1+x) for x>0, otherwise 0
							 UOT_EXP,  // output f(x)=e^x
							 UOT_NEG,  // output f(x)=-x
							 UOT_ABS,  // output f(x)=|x|
							 UOT_ZVAL, // normalize x so according to mean and variance
							 UOT_SPLIT,// split x into several features according to values
									   // if add bool flag is set, each bin also gets a boolean indicator
							 UOT_NUM_OPERATORS,
} UnaryOperatorType;

const char* const unaryOperatorLabels[]={"BOOL",  "LOG",  "EXP", "NEG", "ABS", "ZVAL",  "SPLIT"};
const size_t numUnaryOperatorLabels = sizeof(unaryOperatorLabels)/sizeof(char*);

typedef enum ConditionalOperatorType {
							 COT_AND,     // performs AND of a list of features
							 COT_OR,	  // performs OR of a list of features
							 COT_NOT_AND, // returns 1 unless all features are present
							 COT_NOT_OR,  // returns 1 only if all featureas are not present
							 COT_NUM_OPERATORS,
} ConditionalOperatorType;

const char* const conditionalOperatorLabels[]={"AND",  "OR",  "NOT_AND", "NOT_OR"};
const size_t numConditionalOperatorLabels = sizeof(conditionalOperatorLabels)/sizeof(char*);


// these are the values that get outputted in the result feature if the conditional
// operator evaluates to "true". Usually there will be a feature that ouputs 1 in addition
// to the type of value(s) listed below.
typedef enum ConditionalValueType {
									CVT_BOOL,	// output 1.0
									CVT_SUM,	// output the sum of values being operated on
									CVT_PROD,   // output the product of the values being operated on
									CVT_MAX,    // output the max val
									CVT_MIN,	// output tht minimal val
									CVT_AVG,	//
									CVT_COUNT,	// count how many features had a value (good only for OR)
									CVT_MEDIAN,
									CVT_FEATURE, // output the value of a specific feature
									CVT_NUM_VALUES,
} ConditionalOperatorValueType;

const char* const conditionalValueLabels[]={"BOOL","SUM","PROD","MAX","MIN","AVG","COUNT","MEDIAN","FEATURE"};
const size_t numConditionalValueLabels = sizeof(conditionalValueLabels)/sizeof(char*);


class MlOperator {
public:

	MlOperator() : indInitialized_(false), indUnaryOperator_(true), indAddIndicator_(true),
					type_(MAX_INT), valueType_(MAX_INT), targetFeatureIndex_(MAX_UINT),
					indexOfAdditionalFeatures_(MAX_INT), numAdditionalFeatures_(MAX_INT) {}
	
	void readOperator(const char* line);
	void writeOperator(char* buffer) const;
	inline void applyOperator(Sample& sam, vector<value_t>& expandedVector) const;

	static void initializeLabelNames();

	bool	getIndUnaryOperator() const { return indUnaryOperator_; }
	bool	getIndAddIndicator() const { return indAddIndicator_; }
	size_t  getType() const { return type_; }
	size_t  getValueType() const { return valueType_; }
	size_t	getTargetFeatureIndex() const { return targetFeatureIndex_; }
	size_t  getIndexOfAdditionalFeatures() const { return indexOfAdditionalFeatures_; }
	size_t	getNumAdditionalFeatures() const { return numAdditionalFeatures_; }
	

private:
	
	bool   indInitialized_;
	bool   indUnaryOperator_;
	bool   indAddIndicator_;

	size_t type_;
	size_t valueType_;
	size_t targetFeatureIndex_; // for unary operator the feature value we are looking at, 
								// for conditional operator with FEATURE value type, the feature index that should be used
	size_t indexOfAdditionalFeatures_;
	size_t numAdditionalFeatures_;

	vector<value_t> valueParameters_;
	vector<size_t>  indexParameters_;

	static map<string,size_t> unaryTypeMap_;
	static map<string,size_t> conditionalTypeMap_;
	static map<string,size_t> valueTypeMap_;

	inline void applyUnaryOperator(Sample& sam, vector<value_t>& expandedVector) const;

	inline bool evaluateCondition(const vector<value_t>& featureValues, 
								  size_t numWithValues,
								  size_t numFeatures) const;

	inline void createConditionResults(Sample& sam,
									   vector<value_t>& expandedVector,
									   const vector<value_t>& featureValues, 
									   size_t numWithValues) const;
	
};

class MlOperatorList {
public:

	void readOperatorList(const char* path);

	void writeOperatorList(const char* path) const;

	void applyOperatorList(Sample& sam) const;

	size_t getNumOperators() const { return operators_.size(); }

	size_t getNumFeatures() const { return featureSet_.getNumFeatures(); }


	void createNewOperators(const MlDataSet& ds,
							const vector<size_t>& designatedIndexes,
							size_t maxNumOperators);


private:

	MlFeatureSet		featureSet_;

	vector<MlOperator>	operators_;
	
	mutable vector<value_t>		expandedVector_;
};





#endif

