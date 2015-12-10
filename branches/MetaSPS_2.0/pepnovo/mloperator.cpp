#include "mloperator.h"
#include "mlauxfun.h"
#include "auxfun.h"

map<string,size_t> MlOperator::unaryTypeMap_;
map<string,size_t> MlOperator::conditionalTypeMap_;
map<string,size_t> MlOperator::valueTypeMap_;


void MlOperator::initializeLabelNames()
{
	unaryTypeMap_.clear();
	conditionalTypeMap_.clear();
	valueTypeMap_.clear();

	size_t i;

	for (i=0; i<numUnaryOperatorLabels; i++)
		unaryTypeMap_[unaryOperatorLabels[i]]=i;

	for (i=0; i<numConditionalOperatorLabels; i++)
		conditionalTypeMap_[conditionalOperatorLabels[i]]=i;

	for (i=0; i<numConditionalValueLabels; i++)
		valueTypeMap_[conditionalValueLabels[i]]=i;
}



/*
Syntax of operator line Unary
U   Type   +/-   #add idx_add #params <params>

Syntax of operator line binary
C   Type   ValueType +/- #add idx_add #params <params>

*/
void MlOperator::readOperator(const char *line)
{
	indInitialized_ = false;
	
	istringstream iss(line);
	char c;
	iss >> c;
	if (c != 'U' && c != 'C')
		error("operator line should start with U or B: ",line);
	indUnaryOperator_ = (c == 'C');

	if (indUnaryOperator_)
	{
		string typeStr;
		iss >> typeStr;
		map<string,size_t>::const_iterator it = unaryTypeMap_.find(typeStr);
		if (it == unaryTypeMap_.end())
			error("unrecognized unary type: ", typeStr.c_str());

		type_ = it->second;
	}
	else
	{
		string typeStr, valueStr;
		iss >> typeStr;
		map<string,size_t>::const_iterator it = conditionalTypeMap_.find(typeStr);
		if (it == conditionalTypeMap_.end())
			error("unrecognized conditional type: ", typeStr.c_str());
		type_ = it->second;

		it = valueTypeMap_.find(valueStr);
		if (it == valueTypeMap_.end())
			error("unrecognized value type: ", valueStr.c_str());
		valueType_ = it->second;
	}
	targetFeatureIndex_ = MAX_UINT;
	iss >> targetFeatureIndex_;
	if (targetFeatureIndex_ == MAX_UINT)
		error("need to supply target feature index (use 0 if not applicable to conditional operator)");

	iss >> c;
	if (c != '-' && c != '+')
		error("token should be + or - (for add conditional variable) : ",line);
	indAddIndicator_ = (c == '+');

	numAdditionalFeatures_ = MAX_UINT;
	indexOfAdditionalFeatures_ = MAX_UINT;
	iss >> numAdditionalFeatures_ >> indexOfAdditionalFeatures_;
	if (numAdditionalFeatures_ == MAX_UINT || indexOfAdditionalFeatures_ == MAX_UINT)
		error("bad operator line (additional features or index of additional) : ",line);

	size_t numParams=0;
	iss >> numParams;
	
	if (indUnaryOperator_)
	{
		valueParameters_.resize(numParams,0.0);
		size_t i;
		for (i=0; i<numParams; i++)
			iss >> valueParameters_[i];
	}
	else
	{
		indexParameters_.resize(numParams,0);
		size_t i;
		for (i=0; i<numParams; i++)
			iss >> indexParameters_[i];
	}
	indInitialized_ = true;
}

void MlOperator::writeOperator(char *buffer) const
{
	char* pos = buffer;
	if (indUnaryOperator_)
	{
		pos+= sprintf(pos,"U %s %d %c %d %d %d ",unaryOperatorLabels[type_], targetFeatureIndex_,
			(indAddIndicator_ ? '+' : '-'), numAdditionalFeatures_, 
			indexOfAdditionalFeatures_, valueParameters_.size());
		size_t i;
		for (i=0; i<valueParameters_.size(); i++)
			pos += sprintf(pos," %f",valueParameters_[i]);
		
	}
	else
	{
		pos+= sprintf(pos,"C %s %s %d %c %d %d %d ",conditionalOperatorLabels[type_], 
			conditionalValueLabels[valueType_], targetFeatureIndex_, 
			(indAddIndicator_ ? '+' : '-'), numAdditionalFeatures_, 
			indexOfAdditionalFeatures_, indexParameters_.size());
		size_t i;
		for (i=0; i<indexParameters_.size(); i++)
			pos += sprintf(pos," %d",indexParameters_[i]);
		
	}
}



inline void MlOperator::applyUnaryOperator(Sample& sam, vector<value_t>& expandedVector) const
{
	const value_t featureVal = expandedVector[targetFeatureIndex_];
	if (featureVal == NON_FLOAT) // no feature val
		return;

	if (type_ == UOT_SPLIT)
	{
		const size_t numBins = valueParameters_.size();
		size_t i;
		for (i=0; i<numBins; i++)
			if (featureVal <= valueParameters_[i])
				break;

		const size_t newFeatureIndex = indexOfAdditionalFeatures_+i;
		expandedVector[newFeatureIndex]=featureVal;
		sam.addPair(newFeatureIndex,featureVal);

		if (indAddIndicator_)
		{
			const size_t newBoolIndex = newFeatureIndex + numBins;
			expandedVector[newBoolIndex]=1.0;
			sam.addPair(newBoolIndex,1.0);
		}
		return;
	}

	if (type_ == UOT_ZVAL)
	{
		expandedVector[targetFeatureIndex_] = (featureVal - valueParameters_[0])/valueParameters_[1];
		if (indAddIndicator_)
		{
			expandedVector[indexOfAdditionalFeatures_]=1.0;
			sam.addPair(indexOfAdditionalFeatures_,1.0);
		}
		return;
	}

	if (type_ == UOT_BOOL)
	{
		expandedVector[targetFeatureIndex_]=1.0;
		return;
	}
	
	if (type_ == UOT_EXP)
	{
		expandedVector[targetFeatureIndex_]=exp(featureVal);
		if (indAddIndicator_)
		{
			expandedVector[indexOfAdditionalFeatures_]=1.0;
			sam.addPair(indexOfAdditionalFeatures_,1.0);
		}
		return;
	}

	if (type_ == UOT_NEG)
	{
		expandedVector[targetFeatureIndex_]=(-featureVal);
		if (indAddIndicator_)
		{
			expandedVector[indexOfAdditionalFeatures_]=1.0;
			sam.addPair(indexOfAdditionalFeatures_,1.0);
		}
		return;
	}

	if (type_ == UOT_ABS)
	{
		expandedVector[targetFeatureIndex_]=exp(featureVal);
		if (indAddIndicator_)
		{
			expandedVector[indexOfAdditionalFeatures_]=1.0;
			sam.addPair(indexOfAdditionalFeatures_,1.0);
		}
		return;
	}
	
	if (type_ == UOT_LOG)
	{
		value_t x = featureVal + static_cast<value_t>(1.0);
		expandedVector[targetFeatureIndex_]=static_cast<value_t>( x> 1.0 ? log(x) : 0.0);
		if (indAddIndicator_)
		{
			expandedVector[indexOfAdditionalFeatures_]=1.0;
			sam.addPair(indexOfAdditionalFeatures_,1.0);
		}
		return;
	}
}


inline bool MlOperator::evaluateCondition(const vector<value_t>& featureValues, 
										  size_t numWithValues,
										  size_t numFeatures) const
{
	if (type_ == COT_OR)
		return (numWithValues>0);

	if (type_ == COT_NOT_OR)
		return (numWithValues == 0);

	if (type_ == COT_AND)
		return (numWithValues == numFeatures);

	if (type_ == COT_NOT_AND)
		return (numWithValues < numFeatures);

	return false;
}



inline void MlOperator::createConditionResults(Sample& sam,
											   vector<value_t>& expandedVector,
											   const vector<value_t>& featureValues, 
											   size_t numWithValues) const
{
	
	if (valueType_ == CVT_BOOL)
	{
		expandedVector[indexOfAdditionalFeatures_ ]=1.0;
		sam.addPair(indexOfAdditionalFeatures_,1.0);
		return;
	}
	
	// single result variable
	if (valueType_ <= CVT_FEATURE)
	{
		value_t result=0;
		if (valueType_ == CVT_SUM)
		{
			size_t i;
			for (i=0; i<numWithValues; i++)
				result+=featureValues[i];
		}
		else if (valueType_ == CVT_PROD)
		{
			result =1.0;
			size_t i;
			for (i=0; i<numWithValues; i++)
				result*=featureValues[i];
		}
		else if (valueType_ == CVT_MAX)
		{
			result = MIN_FLOAT;
			size_t i;
			for (i=0; i<numWithValues; i++)
				if (featureValues[i]>result)
					result = featureValues[i];
		}
		else if (valueType_ == CVT_MIN)
		{
			result = MAX_FLOAT;
			size_t i;
			for (i=0; i<numWithValues; i++)
				if (featureValues[i]<result)
					result = featureValues[i];
		}
		else if (valueType_ == CVT_AVG)
		{
			size_t i;
			for (i=0; i<numWithValues; i++)
				result+=featureValues[i];

			result /= numWithValues;
		}
		else if (valueType_ == CVT_COUNT)
		{
			result = static_cast<value_t>(numWithValues);
		}
		else if (valueType_ == CVT_MEDIAN)
		{
			result = featureValues[numWithValues/2];
		}
		else if (valueType_ == CVT_FEATURE)
		{
			result = expandedVector[targetFeatureIndex_];
		}

		// add result 
		if (result != NON_FLOAT && (result == result)) // test for NaN
		{
			expandedVector[indexOfAdditionalFeatures_ ]=result;
			sam.addPair(indexOfAdditionalFeatures_,result);

			if (indAddIndicator_)
			{
				expandedVector[indexOfAdditionalFeatures_+1 ]=1.0;
				sam.addPair(indexOfAdditionalFeatures_+1,1.0);
			}
		}

		return;
	}

	// will some conditions create multiple results?
}


/**************************************************************************
***************************************************************************/
inline void MlOperator::applyOperator(Sample& sam, vector<value_t>& expandedVector) const
{
	if (indUnaryOperator_)
	{
		applyUnaryOperator(sam, expandedVector);
		return;
	}
	else
	{
		// collect values
		static vector<value_t> featureValues;
		const size_t numFeatures = indexParameters_.size();

		if (featureValues.size() < numFeatures)
			featureValues.resize(numFeatures);
	
		size_t numWithValues=0;
		size_t i;
		for (i=0; i<numFeatures; i++)
		{
			const float& v=expandedVector[indexParameters_[i]];
			if (v != NON_FLOAT)
				featureValues[numWithValues++]=v;
		}

		if (evaluateCondition(featureValues, numWithValues, numFeatures))
		{
			createConditionResults(sam, expandedVector, featureValues, numWithValues);
		}
		return;
	}
	
	
}


/*************************************************************************
This function sequentially applies the operators to the feature values in
the sample. Intermidate results are stored in an expanded array (as opposed
to the sparse feature representation in the sample).

An invariant that is maintained:
- Every entry in the expanded array that is > MIN_FLOAT has a corresponding
  feature-value pair in the sample (so if the operator sees such a value in
  the array, it does not add a pair to the sample).
**************************************************************************/
void MlOperatorList::applyOperatorList(Sample& sam) const
{
	if (sam.pairs.size() == 0)
		return;

	// enlarge expanded vector
	if (expandedVector_.size() < getNumFeatures())
	{
		expandedVector_.clear();
		expandedVector_.resize(getNumFeatures(), NON_FLOAT); // NON_FLOAT represents NO VALUE
	}

	// copy original feature values
	size_t i;
	for (i=0; i<sam.pairs.size(); i++)
		expandedVector_[sam.pairs[i].index]=sam.pairs[i].value;

	// apply operators

	for (i=0; i<operators_.size(); i++)
	{
		const MlOperator& mlo = operators_[i];
		// quick test to see if the relevant feature is there. This works for all Unary operators
		// and conditional operators that return the type feature
		if (expandedVector_[mlo.getTargetFeatureIndex()]>MIN_FLOAT) 
			mlo.applyOperator(sam, expandedVector_);
	}

	// copy updated values and revert cells to MIN_FLOAT
	for (i=0; i<sam.pairs.size(); i++)
	{
		const size_t index = sam.pairs[i].index;
		sam.pairs[i].value = expandedVector_[index];
		expandedVector_[index] = NON_FLOAT;
	}	
}





void MlOperatorList::createNewOperators(const MlDataSet& ds,
							const vector<size_t>& missIndexes,
							size_t maxNumOperators)
{

}
