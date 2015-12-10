/*
 * PeptideSpectrumMatchSet.h
 *
 *  Created on: Apr 27, 2011
 *      Author: jsnedecor
 */

#ifndef PEPTIDESPECTRUMMATCHSET_H_
#define PEPTIDESPECTRUMMATCHSET_H_

//Module includes
#include "DelimitedTextReader.h"
#include "PeptideSpectrumMatch.h"
#include "SpecSet.h"

//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/memory>
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <memory>
#  include <unordered_map>
#endif

//System includes
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>

namespace specnets
{

  /**
   * @see PeptideSpectrumMatch.h
   */
  class PeptideSpectrumMatch;

  /**
   * @see SpecSet
   */
  class SpecSet;

  typedef std::tr1::shared_ptr<PeptideSpectrumMatch> psmPtr;

  class PeptideSpectrumMatchSet
  {
  public:

    /*! \brief Return associated  PeptideSpectrumMatch for
     that vector position.
     */
    psmPtr & operator[](unsigned int i);

    /*! \brief Return associated  PeptideSpectrumMatch
     for that vector position.
     */
    const psmPtr & operator[](unsigned int i) const;

    /*! \brief Set value of one PeptideSpectrumMatchSet to another

     */
    PeptideSpectrumMatchSet & operator=(PeptideSpectrumMatchSet &other);

    /*! \brief Adds to psmSet vector
     *
     */
    PeptideSpectrumMatchSet & push_back(psmPtr &other);

    /*! \brief Returns size of m_psmSet vector

     */
    unsigned int size();

    /*! \brief Resizes m_psmSet vector

     @param newSize the new size of the parameter vector
     */
    unsigned int resize(unsigned int newSize);


    bool loadFromFile(const char * filename);
    bool saveToFile(const char * filename);

    /** Inspect result file parsing
     *
     *@param resultsFile = input inspect results
     *@param zeroIndexed = is file zero indexed or not (yes for MGF, no for mzXML)
     */
    bool loadInspectResultsFile(const char * resultsFile, bool zeroIndexed = 1);
    /** Load multiple inspect results files.
     *
     *@param resultsFileList = tab delimited file which contains columns Path and (optionally) isZeroIndexed.
     */
    bool loadInspectResultsFiles(const char * resultsFileList);

    /** MSGFDB result file parsing
     *
     *@param resultsFile = input MSGFDB results
     *@param mergeConsecutive = if true, merge consecutive PSMs with the same precursor and different fragmentation methods
     *@param updateScanIdxs = if loading from mgf and scan numbers are indices and
     *  mergeConsecutive == true, set this to true in order to match the scan nums to indices of the merged SpecSet
     *@param zeroIndexed = is file zero indexed or not (should always be false for MSGFDB)
     */
    bool loadMSGFDBResultsFile(const char * resultsFile,
    		                   bool mergeConsecutive = true,
    		                   bool updateScanIdxs = true,
    		                   bool zeroIndexed = false);

    /** Specnets result file parsing
     *
     *@param resultsFile = input inspect results
     *@param zeroIndexed = is file zero indexed or not (yes for MGF, no for mzXML)
     */
    bool
    loadSpecnetsResultsFile(const char * resultsFile, bool zeroIndexed = 0);

    /**
     * Match up annotations with spectra
     * @param spectra spectra to add to psms
     *  @param addToSpectra = choose whether to update psmList on
     * spectra
     */
    void addSpectra(SpecSet * spectra, bool addToSpectra = true);

    /**
     * Match up annotations with spectra
     *
     * @param spectra spectra to add to psms
     * @param filename = original file name
     * @param addToSpectra = choose whether to update psmList on
     * spectra
     */
    void addSpectra(SpecSet * spectra, string filename, bool addToSpectra =
        true);

    /** Return PSM for a particular scan number
     *
     */
    bool getResultsByScan(int scan,
                          vector<PeptideSpectrumMatchSet>& output_results);

    /** Build PSMSet from input SpecSet
     *  Note: this will overwrite any existing PSMs.
     */
    void getPSMSet(SpecSet * spectra);


    void removePsmSetItem(psmPtr);

    vector<psmPtr> m_psmSet;
  protected:
    bool parseLine(psmPtr currMatch,
                   vector<string> &line,
                   int zeroIndexed,
                   int scanIndex,
                   int specIdxIndex,
                   int spectrumFileIndex,
                   int annotationIndex,
                   int chargeIndex,
                   int scoreIndex,
                   int pvalueIndex,
                   int proteinIndex,
                   bool isInspect);

  };

  class PeptideSpectrumMatchSetSpectralLibraryLoader: public PeptideSpectrumMatchSet
  {
  public:
    bool
    loadSpecnetsResultsFile(const char * resultsFile, bool zeroIndexed = 0);

  protected:
    bool parseLine(psmPtr currMatch,
                   vector<string> &line,
                   int zeroIndexed,
                   int scanIndex,
                   int spectrumFileIndex,
                   int annotationIndex,
                   int chargeIndex,
                   int scoreIndex,
                   int pvalueIndex,
                   int proteinIndex,
                   int strictenvIndex,
                   int unstrictenvIndex);

  };
}

#endif /* PEPTIDESPECTRUMMATCHSET_H_ */
