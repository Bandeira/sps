#include "mloperator_selection.h"
#include "mlauxfun.h"
#include "auxfun.h"

ListStorage::ListStorage(size_t k)
{
	if (! initialize(k) && ! initialize(k/4))
		error("Couldn't initialize a storage area (or even a quarter of it), size in K = ",k);
}

ListStorage::~ListStorage()
{
	if (storageBegin_ && storageEnd_)
	{
		size_t size = storageEnd_ - storageBegin_;
		memset(storageBegin_,0,size*sizeof(size_t)); // invalidate the lists in case some one is still looking
		delete [] storageBegin_;
	}
}

bool ListStorage::initialize(size_t k)
{
	size_t size = k *1024;

	storageBegin_ = new size_t[size];
	if (! storageBegin_)
		return false;

	storageEnd_ = storageBegin_ + size;
	nextEmpty_	= storageBegin_;

	numWrapAround_ =0;
	numAdditions_ =0;
	headerMap_.clear();

	return true;
}


bool ListStorage::addList(const ListHeader& header)
{
	size_t neededSpace = header.length + 3;

	if (nextEmpty_ + neededSpace >= storageEnd_)
	{
		numWrapAround_++;
		nextEmpty_=0;
	}

	if (nextEmpty_ + neededSpace >= storageEnd_)
		return false;

	// write list
	*nextEmpty_++ = listHeaderStartMarker;
	*nextEmpty_++ = numAdditions_;
	

	memcpy(nextEmpty_, header.list, header.length * sizeof(size_t));

	ListHeader newHeader;
	newHeader.listNumber = numAdditions_;
	newHeader.length  = header.length;
	newHeader.list	  = nextEmpty_;
	headerMap_[newHeader]=numAdditions_++;

	nextEmpty_ += header.length;

	return true;
}

bool ListStorage::getList(ListHeader& header) const
{
	map<ListHeader,size_t>::const_iterator it = headerMap_.find(header);
	if (it == headerMap_.end() || ! it->first.isValid())
		return false;

	header = it->first;

	return true;
}


