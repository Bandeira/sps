#ifndef MZXML_H
#define MZXML_H

#include <xercesc/sax2/DefaultHandler.hpp>
#include "spectrum.h"
#include "SpecSet.h"
#include <string>

using namespace std;
/** TODO: add comment */
//#define XERCES_CPP_NAMESPACE_USE

/**
 * TODO: add description
 *
 *@param filename
 *@param specs
 *@param msLevels
 *@param minMsLevel
 *@return
 */

//unsigned int LoadMzxml(char *filename, vector<Spectrum> &specs, vector<unsigned int> &scanNums,
//                 vector<short> *msLevels = 0, short minMsLevel = 2 );

template <typename T>
  unsigned int LoadMzxml(const T & xmlFile, specnets::SpecSet &specs, vector<short> *msLevels = 0, short minMsLevel = 2);

unsigned int LoadMzxml(const char *filename, specnets::SpecSet &specs, const string & sscan, vector<short> *msLevels = 0, short minMsLevel = 2);


/**
 * TODO: add description
 *
 */
class mzxmlSAX2Handler: public XERCES_CPP_NAMESPACE_QUALIFIER DefaultHandler {
protected:
	/**
	 * TODO: add description
	 *
	 */
	static short STATE_SKIPPING;

	/**
	 * TODO: add description
	 *
	 */
	static short STATE_SCAN;

	/**
	 * TODO: add description
	 *
	 */
	static short STATE_PRECURSORMZ;

	/**
	 * TODO: add description
	 *
	 */
	static short STATE_PEAKS;
    

public:

	/**
	 * Set of retained spectra.
	 */
	specnets::SpecSet *specs;

	//	vector<Spectrum> *specs;    // Set of retained spectra
	//	vector<unsigned int> *scanNums; // Scan number for each retained spectrum

	/**
	 * MS level for each retained spectrum.
	 */
	vector<short> *msLevels;

	/**
	 * Spectrum currently being read (list is necessary because of nested scans).
	 */
	list<int> curSpec;

	/**
	 * Indicates current state (see STATE_* above).
	 */
	short state;

	/**
	 * Updated count of retained spectra.
	 */
	unsigned int numSpecs;

	/**
	 * Minimum msLevel to read a spectrum into specs.
	 */
	short minMsLevel;

	/**
	 * TODO: add description
	 *
	 *@param in_specs
	 *@param in_msLevels
	 *@param in_minMsLevel
	 */
	mzxmlSAX2Handler(specnets::SpecSet &in_specs, vector<short> *in_msLevels = 0,
			short in_minMsLevel = 2) {
		numSpecs = 0;
		state = STATE_SKIPPING;
		specs = &in_specs;
		msLevels = in_msLevels;
		minMsLevel = in_minMsLevel;
		curSpec.clear();
	}

	/**
	 * TODO: add description
	 *
	 *@param uri
	 *@param localname
	 *@param qname
	 *@param attrs
	 */
	void startElement(const XMLCh* const uri, const XMLCh* const localname,
			const XMLCh* const qname, const XERCES_CPP_NAMESPACE_QUALIFIER Attributes& attrs);

	/**
	 * TODO: add description
	 *
	 *@param chars
	 *@param length
	 */
	void characters(const XMLCh* const chars, const XMLSize_t length);

	/**
	 * TODO: add description
	 *
	 *@param uri
	 *@param localname
	 *@param qname
	 */
	void endElement(const XMLCh* const uri, const XMLCh* const localname,
			const XMLCh* const qname);

	/**
	 * TODO: add description
	 *
	 */
	void endDocument();

	/**
	 * TODO: add description
	 *
	 */
	void fatalError(const XERCES_CPP_NAMESPACE_QUALIFIER SAXParseException&);
};

#endif
