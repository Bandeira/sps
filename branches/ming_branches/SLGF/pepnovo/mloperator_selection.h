#ifndef __MLOPERATOR_SELECTION_H__
#define __MLOPERATOR_SELECTION_H__

/************************************************************************

This class contains data structures and functions that are used in the
process of operator selection (during the initial training).

*************************************************************************/

#include "mloperator.h"

// pair of two features
// assumed but not inforced: idx1<=idx2
struct IdxPair {
	IdxPair() : idx1(0), idx2(0) {}
	bool operator< (const IdxPair& rhs) const
	{
		return (idx1<rhs.idx1 || (idx1 == rhs.idx1 && idx2<rhs.idx2));
	}
	size_t idx1, idx2;
};


// holds statistics about the co-occurrence of two features
// Let Si= #samples that have feature i
// and T = total weight of sample
struct CoOccurrenceRecord {
	CoOccurrenceRecord() : ratioCoOccurrence(0), relRatioCoOccurrence(0), appearCorrelation(0), vlaueCorrelation(0) {}

	float ratioCoOccurrence;     // = S(i,j)/T
	float relRatioCoOccurrence;  // = S(i,j)/ min{Si,Sj}
	float appearCorrelation;     // = Corr(I[Xi],I[Xj])
	float vlaueCorrelation;      // = Corr(Xi,Xj)
};


const size_t listHeaderStartMarker = 4040404040U; // marks the start of a list (for debug and memory corruption tests)

// Is a pointer to a (usually sorted) list of feature/sample indexes
struct ListHeader {
	bool operator< (const ListHeader& rhs) const
	{
		return (listNumber < rhs.listNumber);
	}

	bool isValid() const { return (*(list-2)== listHeaderStartMarker &&
								   *(list-1)== listNumber); }

	size_t  listNumber;
	size_t  length;
	size_t* list;
};


// Holds a large memory block where lists can be stored.
// Lists are identified by the pair index. Once a list is added
// it might not stay for ever, it might get over written, in which
// case we will nedd to compute it again.

class ListStorage {
public:
	ListStorage() : numWrapAround_(0), numAdditions_(0), storageBegin_(0), 
		storageEnd_(0),  nextEmpty_(0) {}
	ListStorage(size_t k);
	~ListStorage();
								
	bool initialize(size_t k); // size in k=1024 bytes

	// will copy the header info if a valid record is found, otherwise
	// returns NULL. Note that the header returned is not guaranteed to
	// be valid when used if between the time it was given and the time it
	// is accessed, there have been several additional insertions (if the
	// storage area is large, this should not be a problem, use judgement).
	bool getList(ListHeader& header) const;

	// doesn't check if the IdxPair is already there, if it is, it will
	// overwrite the list. Might return false if the list could not be written
	// (maybe too large for storage)
	bool addList(const ListHeader& header);

	void clear() { headerMap_.clear(); numWrapAround_=0; numAdditions_=0; nextEmpty_ = storageBegin_; }

private:
	
	size_t  numWrapAround_;
	size_t	numAdditions_;
	size_t* storageBegin_;
	size_t*	storageEnd_;
	size_t* nextEmpty_;
	map<ListHeader,size_t> headerMap_;
};



/********************************************************************************
Data structure that contains all temporary info from samples that is needed
in the operator generation phase.
*********************************************************************************/
struct MlOperatorSearchData {

	MlOperatorSearchData() {}

	vector<value_t> featureMeans_; // of values (only when the feature is in a sample's list)
	vector<value_t> featureSds_;   

	vector<bool> flagsFeatureInAll_;	// indicates that the feature appears in all the samples
										// (in that case, the occurence list will remain empty)

	vector< vector<size_t> > featureOccurenceLists_; // for each feature holds a sorted list of the indexes
													// of the samples in which it occurs

	map<IdxPair,CoOccurrenceRecord> coOccurrenceRecords_;

	ListStorage featureIntersectinLists_;
};




#endif
