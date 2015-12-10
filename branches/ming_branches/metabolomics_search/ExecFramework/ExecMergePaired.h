/*
 * ExecMergePaired.h
 *
 *  Created on: Dec 12, 2011
 *      Author: aguthals
 */

#ifndef EXECMERGEPAIRED_H_
#define EXECMERGEPAIRED_H_

#include "ExecBase.h"
#include "SpecSet.h"

// System Includes
#include <string>
#include <vector>

using namespace std;

namespace specnets {
class ExecMergePaired: public ExecBase {
public:

	ExecMergePaired(void);

	ExecMergePaired(const ParameterList & inputParams);

	ExecMergePaired(const ParameterList & inputParams, SpecSet * pairedSpectra);

	ExecMergePaired(const ParameterList & inputParams, SpecSet * pairedSpectra,
			SpecSet * mergedSpectra);

	virtual ~ExecMergePaired(void);

	virtual ExecBase * clone(const ParameterList & input_params) const;

	virtual bool invoke(void);

	virtual bool loadInputData(void);

	virtual bool saveOutputData(void);

	virtual bool saveInputData(std::vector<std::string> & filenames);

	virtual bool loadOutputData(void);

	virtual std::vector<ExecBase *> const & split(int numSplit);

	virtual bool merge(void);

	virtual bool validateParams(std::string & error);

private:
	bool ownInput;
	bool ownOutput;
	SpecSet * m_mergedSpecset;
	SpecSet * m_inputSpecset;
};
}

#endif /* EXECMERGEPAIRED_H_ */
