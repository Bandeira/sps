#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <cstring>

#include "vector.h"
#include "aminoacid.h"
#include "spectrum_scoring.h"
#include "utils.h"

using namespace std;

/**
 * @see label.h
 */
class SpectrumPeakLabels;

/**
 * @see ms1.h
 */
class IsoEnvelope;

/**
 * TODO: add description
 */

class Spectrum
{

  /**
   * TODO: add description
   *
   *@param startIdx
   *@param targetMass
   *@param tolerance
   *@param minIdx
   *@param maxIdx
   */
  void locateWithinTolerance(int startIdx,
                             float targetMass,
                             float tolerance,
                             int &minIdx,
                             int &maxIdx);

  /**
   * TODO: add description
   *
   *@param minIdx
   *@param maxIdx
   *@param avgMass
   *@param totScore
   */
  void weightedAverage(int minIdx, int maxIdx, float &avgMass, float &totScore);
public:

  static const unsigned short BIN_VERSION;

  static const unsigned short BIN_SUBVERSION;

  static const string BIN_VERSION_ID;

  static const string BIN_SUBVERSION_ID;

  enum FragType
  {
    FragType_CID = 0,
    FragType_ETD = 1,
    FragType_HCD = 2
  };

  /**
   * TODO: add description
   */
  sps::vector<TwoValues<float> > peakList;

  /**
   * Annotation for each peak in peakList.
   * annotation[i] - contains the annotation for i-th spectrum peak
   * annotation[i].first is a pointer to a structure defining the type of ion (ftIonFragment *)
   * annotation[i].second is the index of the ion in that corresponding ion series (e.g., first b-ion has index 1)
   */
  sps::vector<pair<const ftIonFragment*, short> > annotation;

  /**
   *	vector loaded in from MS2Model for annotation;
   */
  vector<ftIonFragment> ionTypes;

  string annotation_peptide;

  /**
   * Parent M/Z
   */
  double parentMZ;

  //sps::vector<float> annotation;
  /**
   * Monoisotopic parent mass (sum of AA masses + 19).
   */
  float parentMass;

  /**
   * Precursor charge.
   */
  short parentCharge;

  /**
   * Scan number for this spectrum.
   */
  unsigned int scan;

  /**
   * Level of MS
   */
  short msLevel;

  /**
   * Type of MS fragmentation - see enum FragType above.
   */
  FragType msFragType;

  /**
   *  Miscelaneous information (e.g. spectrum filename)
   */
  char *info;

  /**
   * How many daltons between peaks that are 1 integer unit apart.
   */
  float resolution;

  /**
   * Inter-Dalton distance: number of integer units between peaks that
   * are 1 Da apart (idDist=1/resolution), e.g. 1 for regular spectra and
   * 1/RESOLUTION for spectra with masses=round(mass/RESOLUTION),
   * e.g. idDist=10 for RESOLUTION=0.1
   */
  float idDist;

  Spectrum();
  ~Spectrum();

  /**
   * TODO: add description
   *
   *@param i
   *@return
   */
  TwoValues<float> &operator[](unsigned int i)
  {
    return peakList[i];
  }

  const TwoValues<float> &operator[](unsigned int i) const
  {
    return peakList[i];
  }

  /**
   * TODO: add description
   *
   *@return
   */
  unsigned int size() const
  {
    return peakList.size();
  }

  /**
   * Resizes current peak list to another given size.
   *
   *@param newSize the new size for the peakList.
   *@return the new size after the size has been set.
   */
  unsigned int resize(unsigned int newSize)
  {
    annotation.resize(newSize);
    peakList.resize(newSize);
    return peakList.size();
  }

  /**
   * Spectrum assignment
   *
   *@param other
   *@return
   */
  Spectrum &operator=(const Spectrum &other);

  Spectrum(const Spectrum &other);
  /**
   * Copy all member fields except peakList.
   *
   *@param other
   *@return
   */
  Spectrum &copyNP(const Spectrum &other);

  /**
   * TODO: add description
   *
   *@param newResolution
   *@param enforceRounding
   */
  void changeResolution(float newResolution, bool enforceRounding);

  /**
   * TODO: add description
   *
   *@param newResolution
   */
  void setResolution(float newResolution)
  {
    resolution = newResolution;
    idDist = 1 / resolution;
  }

  /**
   * Add peaks for b0/bk/y0/yk
   *
   *@param tolerance
   *@param ionOffset
   *@param includeY0k
   *@param ctermH2O
   *@param labels
   */
  void addZPMpeaks(float tolerance,
                   float ionOffset,
                   bool includeY0k,
                   bool ctermH2O = true,
                   SpectrumPeakLabels *labels = 0);

  /**
   * TODO: add description
   *
   *@param tolerance
   *@param ionOffset
   *@param includeY0k
   */
  void maximizeZPMpeaks(float tolerance, float ionOffset, bool includeY0k =
      false);

  /**
   * Normalizes total intensity to 100.
   *
   *@param newTotalInt
   *@param removeNegatives
   */
  void normalize(float newTotalInt = 100, bool removeNegatives = true);

  /**
   * TODO: add description
   *
   *@param parentMS1
   *@param peakTol
   *@param maxZ
   *@param isoEnvs
   *@param strictMode
   */
  void guessPrecursorZPM(Spectrum &parentMS1,
                         float peakTol,
                         short maxZ,
                         IsoEnvelope &isoEnvs,
                         bool strictMode = false);

  /**
   * Find all peaks within tolerance of baseMass.
   *
   *@param baseMass
   *@param peakTol
   *@param matchesIdx
   *@return
   */
  short findMatches(float baseMass,
                    float peakTol,
                    vector<int> &matchesIdx,
                    int startIdx = -1);

  /**
   * TODO: add description
   *
   *@param mass
   *@return
   */
  int findClosest(float mass);

  /**
   * TODO: add description
   *
   *@param newPeaks
   *@param putHere
   */
  void mergePeakList(vector<TwoValues<float> > &newPeaks, Spectrum *putHere);

  /**
   * TODO: add description
   *
   *@param newPeaks
   *@param putHere
   */
  void mergePeakListRev(vector<TwoValues<float> > &newPeaks, Spectrum *putHere);

  /**
   * TODO: add description
   *
   *@param pmOffset
   *@param putHere
   */
  void reverse(float pmOffset, Spectrum *putHere = 0);

  /**
   * TODO: add description
   *
   *@param offset
   *@param cterm
   */
  void rotate(float offset, float cterm);

  /**
   * TODO: add description
   *
   *@param topK
   *@param putHere
   */
  void selectTopK(unsigned int topK, Spectrum *putHere = 0);

  /**
   * TODO: add description
   *
   *@param topK
   *@param putHere
   */
  void selectTopK(unsigned int topK, TwoValues<float> w, Spectrum *putHere = 0);

  /**
   * TODO: add description
   *
   *@param resolution
   */
  void roundMasses(float resolution = 0.1);

  /**
   * Converts a list of masses to a list of corresponding spectrum peak
   * indices (closest mass).
   *
   *@param masses
   *@param indices
   *@param peakTol
   */
  void massesToIndices(vector<TwoValues<float> > &masses,
                       vector<int> &indices,
                       float peakTol);

  /**
   * Retain only the peaks with indices in idx (idx MUST be sorted).
   *
   *@param idx
   */
  void selectIndices(vector<int> &idx);

  /**
   * Removes peaks with intensity smaller than minIntensity
   *
   *@param minIntensity
   *@return
   */
  unsigned int filterLowInt(float minIntensity = 0);

  /**
   * Retains only peaks with masses in [minMass,maxMass] (inclusive)
   *
   *@param minMass Minimum peak mass
   *@param maxMass Maximum peak mass; set to -1 to disable maxMass filter
   *@param peakTol Peak mass tolerance
   */
  void filterMasses(float minMass, float maxMass = -1, float peakTol = 0.5);

  /**
   * TODO: add description
   *
   *@param output
   */
  void output(ostream &output);

  /**
   * TODO: add description
   *
   *@param output
   */
  void output_ms2(ostream &output);

  /**
   * TODO: add description
   *
   *@param pmOffset
   *@param tolerance
   *@param pairs
   *@param pairsIdx
   */
  void getPairs(float pmOffset,
                float tolerance,
                vector<vector<float> > &pairs,
                vector<vector<int> > &pairsIdx);

  /**
   * TODO: add description
   *
   *@param pmOffset
   *@param tolerance
   */
  void makeSymmetric(float pmOffset, float tolerance);

  // Generic functions but initially defined for alignment of consensus to PRM spectra (batch.cpp::getPairAlignsPAext_aux)
  /**
   * TODO: add description
   *
   *@param toSpec
   *@return
   */
  bool compare(Spectrum &toSpec);

  /**
   * TODO: add description
   *
   *@param withPeaks
   *@param toSpec
   *@param scores1
   *@param scores2
   *@param scoresMerged
   */
  void merge(vector<TwoValues<float> > &withPeaks, Spectrum &toSpec, vector<
      float> &scores1, vector<float> &scores2, vector<float> &scoresMerged);

  /**
   * TODO: add description
   *
   *@param withSpectrum
   *@param toSpec
   *@param shift
   *@param peakTol
   *@param mergeType
   */
  void mergeCommon(Spectrum &withSpectrum,
                   Spectrum *toSpec,
                   float shift,
                   float peakTol,
                   short mergeType = 0);

  /**
   * TODO: add description
   *
   *@param other
   *@param shift
   *@param tolerance
   *@param output
   *@param otherScores
   *@param otherEPScores
   */
  void filterAdd(Spectrum &other,
                 float shift,
                 float tolerance,
                 Spectrum &output,
                 vector<float> &otherScores,
                 float &otherEPScores);

  /**
   * De-novo intepretation functions.
   *
   *@param jumps
   *@param peakTol
   *@param matchedIdx
   *@return
   */
  float denovoLR(AAJumps &jumps, float peakTol, vector<int> &matchedIdx);

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
   */
  void annotate(string &peptide,
                string &ionNamesInclude,
                MS2ScoringModel &inputIonTypes,
                float prmOffset,
                float srmOffset,
                float peakTol);

  bool saveDTA(const char* outfile);

  bool loadFromBinaryStream(FILE* fp, map<string, unsigned short>& versions);
};

/**
 * TODO: add description
 */
class SpecSet
{
public:

  /**
   * TODO: add description
   */
  sps::vector<Spectrum> specs;

  /**
   * TODO: add description
   *
   *@param sz
   */
  SpecSet(unsigned int sz = 0)
  {
    specs.resize(sz);
  }

  /**
   * TODO: add description
   *
   *@param filename
   */
  SpecSet(char *filename)
  {
    LoadSpecSet_pkl(filename);
  }

  /**
   * TODO: add description
   *
   *@param i
   *@return
   */
  Spectrum &operator[](unsigned int i)
  {
    return specs[i];
  }

  const Spectrum & operator[](unsigned int i) const
  {
    return specs[i];
  }

  /**
   * Get scan at scan number
   * @param scan_num - scan number (indexed from 1!)
   */
  Spectrum* getScan(int scan_num);

  /**
   * TODO: add description
   *
   *@return
   */
  unsigned int size()
  {
    return specs.size();
  }

  /**
   * TODO: add description
   */
  unsigned int resize(unsigned int newSize)
  {
    specs.resize(newSize);
    return specs.size();
  }

  /**
   * TODO: add description
   *
   *@param results
   *@param resolution
   *@return
   */
  vector<list<int> > &getMassesHistogram(vector<list<int> > &results,
                                         float resolution = 1.0);

  /**
   * TODO: add description
   *
   *@param newResolution
   *@param enforceRounding
   */
  void changeResolution(float newResolution, bool enforceRounding)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].changeResolution(newResolution, enforceRounding);
  }

  /**
   * TODO: add description
   *
   *@param tolerance
   *@param ionOffset
   *@param includeY0k
   */
  void addZPMpeaks(float tolerance, float ionOffset, bool includeY0k)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].addZPMpeaks(tolerance, ionOffset, includeY0k);
  }

  /**
   * TODO: add description
   *
   *@param newTotalInt
   */
  void normalize(float newTotalInt = 1)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].normalize(newTotalInt);
  }

  /**
   * TODO: add description
   *
   *@param features
   *@param featureValue
   *@param output
   *@return
   */
  template<class T> unsigned int extract(vector<T> &features,
                                         T featureValue,
                                         SpecSet &output);

  bool LoadSpecSet_pkl_mic(const char* filename);
  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  unsigned int LoadSpecSet_pkl(const char *filename);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  unsigned int LoadSpecSet_ms2(const char *filename);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  unsigned int LoadSpecSet_mgf(const char *filename);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  unsigned int LoadSpecSet_prms(const char *filename);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  unsigned int LoadSpecSet_prmsv3(const char *filename);

  /**
   * VERSION 1: Loads all data fields of each spectrum from an open pklbin file
   * @param fp Open file pointer where the number of spectra and version of
   *    already been read from. This will be closed upon return
   * @param numSpecs Number of spectra to read. This has already been read from the file
   * @param oldVersion true if reading in the original pklbin format (pre version 1).
   * @param subversion Sub-version of this format. Currently not used (at version 1, subversion 0)
   * @return true if the SpecSet was successfully loaded, false if not. In either case the file will
   *    be closed
   */
  bool loadPklBin_1(FILE* fp, int numSpecs, bool oldVersion, char subversion);

  /**
   * VERSION 2: Loads all data fields of each spectrum from an open pklbin file
   * @param fp Open file pointer where the number of spectra and version of
   *    already been read from. This will be closed upon return
   * @param numSpecs Number of spectra to read. This has already been read from the file
   * @param subversion Sub-version of this format. Currently not used (at version 2, subversion 0)
   * @return true if the SpecSet was successfully loaded, false if not. In either case the file will
   *    be closed
   */
  bool loadPklBin_2(FILE* fp, int numSpecs, char subversion);

  /**
   * VERSION 3: Loads all data fields of each spectrum from an open pklbin file
   * @param fp Open file pointer where the number of spectra and version of
   *    already been read from. This will be closed upon return
   * @param numSpecs Number of spectra to read. This has already been read from the file
   * @param subversion Sub-version of this format. Currently not used (at version 2, subversion 0)
   * @return true if the SpecSet was successfully loaded, false if not. In either case the file will
   *    be closed
   */
  bool loadPklBin_3(FILE* fp, int numSpecs, char subversion);

  /**
   * Loads all data fields of each spectrum in a pklbin file. Backwards compatible with all old versions.
   *
   *@param filename
   *@param psmFileame
   *@param peaksFilename
   *@return
   */
  int loadPklBin(const char * filename,
                 const char * psmFileame = 0x0,
                 const char * peaksFilename = 0x0);

  /**
   * LoadSpecSet_pklbin - loads a set of spectra in binary format. File format
   * 1 int - number of spectra in the file
   * numSpecs shorts - number of peaks per spectrum in the file.
   * arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
   *
   *@param filename Name of the pklbin file to load
   *@param countSpectraOnly If then only counts number of spectra in the file without loading the spectra (defaults to false)
   *@return
   */
  unsigned int LoadSpecSet_pklbin(const char *filename, bool countSpectraOnly =
      false);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  short SaveSpecSet_info(const char *filename);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  short SaveSpecSet_mgf(const char *filename, const char* activation = 0);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  short SaveSpecSet_pkl(const char *filename);

  /**
   * Saves the annotation_peptide field of each spectrum to the indexed line in output file
   *@param filename output filename, if a spectrum has an annotation_peptide it will be written to the same line # as its index
   *@return true if file was written successfully, false otherwise
   */
  bool SaveAnnotations(const char* filename);

  short SaveSpecSet_ms2(const char* filename);

  /**
   * TODO: add description
   *
   *@param filename
   *@param filename
   *@return
   */
  short
  SaveSpecSet_pklbin(const char *filename, const char *binFilename = NULL);

  /**
   * TODO: Current pklbin format. Includes scan # and msLevel.
   *
   *@param filename
   *@return
   */
  int savePklBin(const char *filename,
                 const char * psmFileame = 0x0,
                 const char * peaksFilename = 0x0);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  short SaveSpecSet_pklbin2(const char *filename);

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  void SaveScanNums(const char *filename);
};

TwoValues<int> countConsecutiveMP(Spectrum& contig1,
                                  Spectrum& contig2,
                                  float shiftFor,
                                  float shiftRev,
                                  float pmTol,
                                  float minOverlapArea,
                                  bool debug = false);

TwoValues<float> getContigOverlapScores(Spectrum& contig1,
                                        Spectrum& contig2,
                                        float shiftFor,
                                        float shiftRev,
                                        float pmTol,
                                        float minOverlapArea,
                                        bool debug = false);

TwoValues<int> getContigOverlapPeaks(Spectrum& contig1,
                                     Spectrum& contig2,
                                     float shiftFor,
                                     float shiftRev,
                                     float pmTol,
                                     float minOverlapArea,
                                     bool debug = false);

/**
 * TODO: add description
 *
 *@param filename
 *@param specs
 *@return
 */
short SaveSpecSet_pklbin(char *filename, vector<Spectrum *> specs);

/**
 * TODO: add description
 *
 *@param inSpec
 *@param parentMass
 *@param peakTol
 *@param outSpec
 */
void MakeSymmetric(vector<TwoValues<float> > *inSpec,
                   float parentMass,
                   float peakTol,
                   vector<TwoValues<float> > *outSpec = 0,
                   vector<int> *indices = 0);

/**
 * TODO: add description
 *
 *@param spec1
 *@param spec2
 *@param spec3
 *@param shift12
 *@param shift13
 *@param peakTol
 *@param minPeakCount
 *@param mergedSpec
 *@param peaksIndices
 *@return
 */
unsigned int MergeSpectra3(Spectrum &spec1,
                           Spectrum &spec2,
                           Spectrum &spec3,
                           float shift12,
                           float shift13,
                           float peakTol,
                           int minPeakCount,
                           Spectrum &mergedSpec,
                           vector<vector<int> > &peaksIndices);

#endif
