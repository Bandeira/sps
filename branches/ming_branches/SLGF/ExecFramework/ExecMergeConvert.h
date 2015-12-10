/*
 * ExecMergeConvert.h
 *
 *  Created on: Oct 5, 2011
 *      Author: aguthals
 */

#ifndef EXECMERGECONVERT_H_
#define EXECMERGECONVERT_H_

#include "ExecBase.h"

#include "mzxml.h"
#include "utils.h"
#include "abruijn.h"
#include "spectrum.h"

// System Includes
#include <string>
#include <vector>

using namespace std;

namespace specnets {
class ExecMergeConvert: public ExecBase {
public:

	static bool loadSpecsetMultiple(const string& specFilePaths,
			SpecSet* spectra, const char* separator = ":");

	static bool loadSpecset(const string& specFilePath, SpecSet* spectra);

	static bool saveSpecset(const string& specFilePath, SpecSet* spectra);

	static pair<bool, list<string> > getFileList(const string& listFilePath);

	vector<pair<string, int> > m_recordedInput;
	vector<pair<string, int> > m_recordedOutput;
	stringstream m_recordedProcessing;
	string m_separator;

	ExecMergeConvert(void);

	ExecMergeConvert(const ParameterList & inputParams);

	ExecMergeConvert(const ParameterList & inputParams,
			vector<pair<int, int> >* inputSpecsetIdx);

	ExecMergeConvert(const ParameterList & inputParams,
			vector<pair<int, int> >* inputSpecsetIdx, SpecSet * outputSpecset);

	virtual ~ExecMergeConvert(void);

	virtual ExecBase * clone(const ParameterList & input_params) const;

	virtual bool invoke(void);

	virtual bool loadInputData(void);

	virtual bool saveOutputData(void);

	virtual bool saveInputData(std::vector<std::string> & filenames);

	virtual bool loadOutputData(void);

	virtual std::vector<ExecBase *> const & split(int numSplit);

	virtual bool merge(void);

	virtual bool validateParams(std::string & error);

	virtual bool saveActivity(const string& infoFilePath);

private:
	vector<pair<int, int> >* m_inputSpecsetIdx; // Set of input specsets
	bool ownInput; //! Does this object "own" the input data structures (and hence have to free them)

	SpecSet * m_outputSpecset; //! Output contig/contig alignments
	bool ownOutput; //! Does this object "own" the output data pointers (and hence have to free them)
};
}

#endif /* EXECMERGECONVERT_H_ */
