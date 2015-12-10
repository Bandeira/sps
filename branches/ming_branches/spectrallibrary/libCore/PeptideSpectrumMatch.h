/*
 * PeptideSpectrumMatch.h
 *
 *  Created on: Apr 5, 2011
 *      Author: jsnedecor
 */

#ifndef PEPTIDESPECTRUMMATCH_H_
#define PEPTIDESPECTRUMMATCH_H_

//Module includes
#include "utils.h"
#include "spectrum.h"
#include "spectrum_scoring.h"
#include "aminoacid.h"

//System includes
#include <iostream>
#include <fstream>
#include <vector>
//TR1 includes. GCC 4.0 and above only!
#ifdef __GLIBCXX__
#  include <tr1/unordered_map>
#else
#  ifdef __IBMCPP__
#    define __IBMCPP_TR1__
#  endif
#  include <unordered_map>
#endif

namespace specnets
{
  /**
   * @see spectrum.h
   */
  class Spectrum;

  /*! \brief Peptide spectrum match class. A single spectrum is matched to a single annotation
   *
   * It is possible to have more than one annotation per spectrum, however, they will show up
   * as different PSMs.
   */
  class PeptideSpectrumMatch
  {
  public:
    //! \name CONSTRUCTORS
    //@{
    /*! \brief constructor for PSM class
     *
     */
    PeptideSpectrumMatch();
    //@}

    //! \name DESTRUCTOR
    //@{
    virtual ~PeptideSpectrumMatch();
    //@}

    /*! \brief set PSM to another PSM, copies peak annotations, spectrum pointers, etc.
     *
     */
    PeptideSpectrumMatch(const PeptideSpectrumMatch &other);

    /*! \brief set PSM to another PSM, copies peak annotations, spectrum pointers, etc.
     *
     */
    virtual PeptideSpectrumMatch &operator=(const PeptideSpectrumMatch &other);

    /** helper function for annotate. Runs through match vector and
     *sets annotation vector
     *@param matches - vector of ion indices (from 0) to set matches for
     *@param annotation - vector of fragments to set
     *@param ionIdx - index of ion 0..N-1 where N is length of peptide. i.e. for b1, ionIdx = 1
     *@param currIonFrag - pointer to current ftIonFragment we're currently annotating
     */
    virtual void
        setAnnotationToMatches(vector<int> &matches,
                               vector<pair<const ftIonFragment*, short> > &annotation,
                               int ionIdx,
                               const ftIonFragment* currIonFrag);

    /** helper function for annotate. Runs through match vector and
     *sets annotation vector, taking only the higher intensity match in case of duplicates
     *@param matches - vector of ion indices (from 0) to set matches for
     *@param annotation - vector of fragments to set
     *@param ionIdx - index of ion 0..N-1 where N is length of peptide. i.e. for b1, ionIdx = 1
     *@param currIonFrag - pointer to current ftIonFragment we're currently annotating
     */
    virtual void
        setAnnotationToMatchesNoDups(vector<int> &matches,
                                     vector<pair<const ftIonFragment*, short> > &annotation,
                                     int ionIdx,
                                     const ftIonFragment* currIonFrag);

    void setDbMatch(string & protein, int dbIndex, float startMass);

    void getDbMatch(string & protein, int & dbIndex, float & startMass);

    /**
	 * add annotations for matched peaks from MS2ScoringModel (using per-peak tolerances)
	 * @param peptide amino acid sequence used to determine peak annotations
	 * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
	 * then just include all ions in MS2Model. ex. "y,b,y++,b++"
	 * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
	 * is copied into ionTypes
	 * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
	 * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
	 * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
	 * @param jumps - AAJumps indicating amino acid masses
	 */
	virtual bool annotate(string &peptide,
			              string &ionNamesInclude,
			              MS2ScoringModel &inputIonTypes,
			              float prmOffset,
			              float srmOffset,
			              AAJumps &jumps,
			              bool removeDuplicates = true,
			              bool ignoreParentCharge = false);

    /**
     * add annotations for matched peaks from MS2ScoringModel
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     * @param jumps - AAJumps indicating amino acid masses
     */
    virtual bool annotate(string &peptide,
                          string &ionNamesInclude,
                          MS2ScoringModel &inputIonTypes,
                          float prmOffset,
                          float srmOffset,
                          float peakTol,
                          AAJumps &jumps,
                          bool removeDuplicates = true,
                          bool ignoreParentCharge = false);

    /**
     * add annotations for matched peaks from MS2ScoringModel
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param removeDuplicates - if true, no duplicate annotations are allowed for the same peak
     */
    virtual bool annotate(string &peptide,
                          string &ionNamesInclude,
                          MS2ScoringModel &inputIonTypes,
                          float prmOffset,
                          float srmOffset,
                          float peakTol,
                          bool removeDuplicates = true,
                          bool ignoreParentCharge = false);
    /**
     * Count the number of annotated peaks and the total explained intensity (must call annotate() to annotate peaks first)
     * @param ionNamesInclude comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     *   then just include all ions in MS2Model. ex. "y,b,y++,b++".
     * @param outputIonCounts optional output map pairing each ion type's name w/ annotated peaks to a pair of: pair.first = number of peaks
     *   annotated to ion type, pair.second = total intensity of matched peaks
     * @return pair.first = total number of matching peaks, pair.second = total matched intensity

    virtual pair<int, float> countAnnotatedPeaks(string &ionNamesInclude,
    		                         map<string, pair<int, float> >* outputIonCounts = 0);*/

    /** Maps peakList indices to ion names.
     *
     *
     * @param outputMatches output for matched indices
     * @param ionMap mapping between ion name + ion index keyed to
     * index of outputMatches
     */
    virtual void mapIons(vector<vector<int> > &outputMatches,
                         std::tr1::unordered_map<string, int> &ionMap);

    /** Inspect annotation parsing
     * Changes Inspect annotations to SpecNets (should add converse function)
     * @param inspect = inspect string
     * @param specnets = output string for specnets
     */
    static void inspectToSpecNets(string &inspect, string &specnets);

    /** Return whether an annotation contains a modification
     *
     *  This simply looks for the first "(" character in the string
     *  if it's there, we assume the peptide is modified
     */
    bool isModified(void);

    /** Returns all associated modifications for annotation
     *
     * This might be worth replacing with something that parses modifications
     * when m_annotation is changed.
     * @param modifications = output vector for modifications
     * @return = returns false when modification is unable to be parsed
     */
    bool getModifications(vector<float> &modifications);

    /**
     * Returns unmodified version of the peptide
     *
     * @param inputPeptide = modified input peptide sequence
     * @param outputPeptide = output unmodified annotation
     */
    static void getUnmodifiedPeptide(string &inputPeptide,
                                     string &outputPeptide);

    /**
     * Returns unmodified version of the peptide.
     *
     * Note that this does not remove gaps, only modifications. (denoted by ())
     * @param outputPeptide = output unmodified annotation.
     */
    void getUnmodifiedPeptide(string &outputPeptide);

    /**
     * Inserts modifications into peptide
     *
     * @param modifications = array of modifications to put into peptide
     * @param positions = locations to insert modifications into
     * @param outputPeptide = output annotation
     * @return = returns false when modifications aren't able to be inserted.
     */
    bool insertModifications(vector<float> &modifications,
                             vector<unsigned int> &positions,
                             string &outputPeptide);
    /**
     * Inserts modifications into peptide
     *
     * @param unmodifiedPeptide = input string
     * @param modifications = array of modifications to put into peptide
     * @param positions = locations to insert modifications into
     * @param outputPeptide = output annotation
     * @return = returns false when modifications aren't able to be inserted.
     */
    static bool insertModifications(string &unmodifiedPeptide,
                                    vector<float> &modifications,
                                    vector<unsigned int> &positions,
                                    string &outputPeptide);

    /**
     * Generates all theoretical masses from model and peptide sequence
     *
     * @param peptide amino acid sequence used to determine peak annotations
     * @param ionNamesInclude  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
     * then just include all ions in MS2Model. ex. "y,b,y++,b++"
     * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
     * is copied into ionTypes
     * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
     * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
     * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
     * @param jumps - AAJumps indicating amino acid masses
     * @param ionNames - output ion names, i.e. b2, y1 etc.
     * @param theoreticalMasses - output theoretical masses for ionNames
     */
    static bool generateTheoreticalMasses(string &peptide,
                                          string &ionNamesInclude,
                                          MS2ScoringModel &inputIonTypes,
                                          float prmOffset,
                                          float srmOffset,
                                          AAJumps &jumps,
                                          vector<string> &ionNames,
                                          vector<float> &theoreticalMasses);


    static bool generateTheoreticalMasses(string &peptide,
                                              string &ionNamesInclude,
                                              MS2ScoringModel &inputIonTypes,
                                              float prmOffset,
                                              float srmOffset,
                                              vector<string> &ionNames,
                                              vector<float> &theoreticalMasses);

    bool loadFromFile(std::ifstream & ifs);
    bool saveToFile(std::ofstream & ofs);

    string m_spectrumFile; //!original spectrum file from which annotation was loaded. Optional
    int m_scanNum; //! Scan number of spectrum.
    string m_annotation; //! Peptide annotation. Uses SpecNets format regardless of original input format.
    string m_origAnnotation; //! Original annotation, if original format has been changed.
    /**
     * Annotation for each peak in peakList.
     * annotation[i] - contains the annotation for i-th spectrum peak
     * annotation[i].first is a pointer to a structure defining the type of ion (ftIonFragment *)
     * annotation[i].second is the index of the ion in that corresponding ion series (e.g., first b-ion has index 1)
     */
    vector<pair<const ftIonFragment*, short> > m_peakAnnotations;

    vector<TwoValues<int> > m_matchedPeaks;

    /**
     *  vector loaded in from MS2Model for annotation;
     */
    vector<ftIonFragment> m_ionTypes;
    string m_protein; //! Protien annotation.
    int m_dbIndex; //! Protein index in DB (if needed)
    int m_numMods; //! Number of modifications
    int m_matchOrientation; //! Forward (0) or reverse (1) tag match orientation
    float m_startMass; //! Starting mass of protein match
    int m_charge; //! Peptide charge
    float m_score; //! associated score
    double m_pValue; //! associated p-value.
    bool m_isDecoy; //! indicates whether hit is decoy or not
    float m_strict_envelope_score; //! strict envelope score of MS1
    float m_unstrict_envelope_score; //! strict envelope score of MS2
    Spectrum * m_spectrum; //! associated spectrum match
    
    /**
     *  Below is meta data regarding the spectrum used for loading spectra
     *  into CCMS for creating spectral libraries
     */
    string m_submission_user;
    string m_submission_id;
    string m_submission_date;
    string m_organism;
    float  m_molecule_mass;
    string m_compound_name;

  private:
    void internalCopy(const PeptideSpectrumMatch &other);
  };
}
#endif /* PEPTIDESPECTRUMMATCH_H_ */
