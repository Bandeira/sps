/**
 * Parsing for spectrum annotation statistics
 */

#ifndef SPECTRUM_ANNOT_STATISTICS_PARSE_H_
#define SPECTRUM_ANNOT_STATISTICS_PARSE_H_

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class SpectrumAnnotParam {
public:
	/** Spectrum annotation parameter  parsing
	 *	@param header - string vector of column names for input
	 *	@param fields - string vector of current line split by "\t"
	 */
	bool parseFromFields(const vector<string>& header, const vector<string>& fields);

	string	ionNames;
	string	statistic;
};

class SpectrumAnnotParams {
public:
	vector<SpectrumAnnotParam> params;

	bool loadSpectrumAnnotFile(const char *spectrum_annot_file);

	SpectrumAnnotParam &operator[](unsigned int i) {
        return params[i];
    }
    const SpectrumAnnotParam &operator[](unsigned int i) const {
        return params[i];
    }

    SpectrumAnnotParams &operator=(SpectrumAnnotParams &other) {
    	params.resize(other.params.size());
    	params = other.params;
    }

    unsigned int size() {
        return params.size();
    }

    unsigned int resize(unsigned int newSize) {
        params.resize(newSize);
        return params.size();
    }
};

#endif /* SPECTRUM_STATISTICS_H_ */
