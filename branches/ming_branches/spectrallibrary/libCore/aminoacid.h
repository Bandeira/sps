#ifndef AMINOACID_H
#define AMINOACID_H

#include "utils.h"
#include "inputParams.h"
#include "mzrange.h"

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>

namespace specnets
{
  class Spectrum; // defined in spectrum.h

  //deprecated, use AAJumps::getPRMMasses
  void getMassesCummulative(const char *sequence,
                            vector<float> &masses,
                            float offset);

  //deprecated, use AAJumps::getPRMMasses
  void getMassesCummulativeNoShift(char *sequence, vector<float> &masses);

  /**
   * Transform zero terminated string into vector of amino acid masses.
   *
   *@param sequence
   *@param masses
   */
  void getMasses(char *sequence, vector<float> &masses);

  /**
   * Transform zero terminated string into vector of amino acid masses.
   *
   *@param sequence
   *@param masses
   *
   */
  void getMasses(vector<char> &sequence, vector<float> &masses);

  /**
   * From a c string get all characters that correspond to the mass of an amino acid
   *@param sequence char* to input
   *@param destination char* to output
   *@return
   */
  void getPepSeq(const char* sequence, char* destination);

  /**
   * Transform zero terminated string into vector of amino acid masses.
   *
   *@param sequence
   *@param masses
   *
   */
  void getMasses(vector<char> &sequence, vector<float> &masses);

  /**
   * Define sets of valid mass jumps based on amino acid and modifications masses.
   */
  class AAJumps
  {
  protected:

    /**
     * Global reference amino acid masses, initialized from const AAmasses
     * (see .cpp).
     *
     */
    static vector<float> glbMasses;

    /**
     * Global reference amino acid letters, initialized from const AAletters
     * (see .cpp).
     *
     */
    static vector<char> glbLetters;

    /**
     * Reference amino acid masses, initialized from const AAmasses (see .cpp)
     */
    vector<float> refMasses;

    /**
     * TODO: add description
     *
     */
    vector<char> refLetters;
  public:

    /*
     static const short NO_MODS=0, USE_MODS=1;
     static const double massHion =  1.0072763,
     massH2O  = 18.010564686,
     minAAmass= 57.0214637230,
     massMH = 18.010564686+1.0072763;
     */

    /**
     * TODO: add description
     *
     */
    static const short NO_MODS;

    /**
     * TODO: add description
     *
     */
    static const short USE_MODS;

    /**
     * TODO: add description
     *
     */
    static const double massHion;

    /**
     * TODO: add description
     *
     */
    static const double minAAmass;

    /**
     * TODO: add description
     *
     */
    static const double massMH;

    /**
     * TODO: add description
     *
     */
    static const double massH2O;

    /**
     * TODO: add description
     *
     */
    static const double massNH3;

    /**
     * TODO: add description
     *
     */
    static const double massCO;

    static const double massNH;

    /**
     * TODO: add description
     *
     */
    short modsUsed;

    /**
     * List of valid amino acid jumps' masses.
     */
    vector<float> masses;

    //    vector< vector<float> > massesMods; // List of modification masses that can be added to entries in masses
    //    vector<list<string> > letters;

    /**
     * Amino acid letters for the single-residue jumps.
     */
    vector<char> aaLetters;

    /**
     * Position i contains the start/end indices of jumps within tolerance of mass
     * i*resolution (see computeIndex() below)
     *
     * Note: index is not valid for negative-mass jumps
     */
    vector<TwoValues<unsigned short> > index;

    /**
     * Default constructor, equivalent to calling other constructor as AAJumps(0)
     */
    AAJumps();

    /**
     * TODO: add description
     *
     *@param maxJumpSize
     *@param resolution
     *@param peakTol
     *@param useMods
     */
    AAJumps(short maxJumpSize,
            float resolution = 0.01,
            float peakTol = -1,
            short useMods = NO_MODS,
            bool uniqueMasses = true);

    /**
     * TODO: add description
     *
     *@param maxJumpSize
     *@param resolution
     *@param peakTol
     *@param useMods
     */
    void getjumps(short maxJumpSize, float resolution = 0.01, float peakTol =
        -1, short useMods = NO_MODS, bool uniqueMasses = true);

    /**
     * Quick test function to check AA refmass
     */
    bool getAAref(const char aa, float &mass);

    /**
     * TODO: add description
     *
     *@param maxJumpMass
     *@param resolution
     *@param peakTol
     *@param useMods
     */
    void alljumps(float maxJumpMass, float resolution = 0.01, float peakTol =
        -1, short useMods = NO_MODS);

    /**
     * TODO: add description
     *
     *@param newJumps
     *@param newNames
     */
    void addJumps(vector<float> &newJumps, vector<char> *newNames = 0);

    /**
     * TODO: add description
     *
     *@param filename
     *@param setGlobal
     *@return
     */
    bool loadJumps(const char *filename, bool setGlobal = false);

    /**
     * TODO: add description
     *
     *@param i
     *@return
     */
    float &operator[](const int i)
    {
      return masses[i];
    }

    /**
     * TODO: add description
     *
     *@return
     */
    const unsigned int size()
    {
      return masses.size();
    }

    /**
     * TODO: add description
     *
     *@param resolution
     */
    void forceUnique(float resolution = 0.01)
    {
      masses = Utils::unique(masses, resolution);
      index.resize(0);
    }

    /**
     * Makes masses = [-masses; masses];
     */
    void forceDoubleSided();

    /**
     * Adds a given mass to the set of valid jumps.
     *
     *@param mass
     */
    void forceJump(float mass)
    {
      masses.push_back(mass);
      sort(masses.begin(), masses.end());
      index.resize(0);
    }

    /**
     * Adds a set of given masses to the set of valid jumps.
     *
     *@param newMasses
     */
    void forceJumps(vector<float> newMasses);

    /**
     * Every mass m in masses is replaced by a set of masses
     * m+[-tolerance:resolution:tolerance].
     *
     *@param tolerance
     *@param resolution
     */
    void forceTolerance(float tolerance, float resolution = 0.01);

    /**
     * TODO: add description
     *
     *@param largestJump
     */
    void removeHigherJumps(float largestJump);

    /**
     * Test whether a given mass is a valid jump in the current set.
     *
     *@param mass
     *@param tolerance
     *@return
     */
    bool isValid(float mass, float tolerance);

    /**
     * TODO: add description
     *
     *@param mass
     *@param tolerance
     *@param idxBounds
     *@return
     */
    unsigned int find(float mass,
                      float tolerance,
                      TwoValues<unsigned int> &idxBounds);

    /**
     * TODO: add description
     *
     *@param peakTol
     *@param resolution
     *@param maxMass
     */
    void computeIndex(float peakTol, float resolution, float maxMass = -1);

    /** AAJumps::getPRMMasses
     * Calculates PRM masses for input annotation in SpecNets format
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param masses: 	Vector of masses to be modified.
     * @param offset: offset for prefix masses
     * @param addZeroMass: if true add offset mass to beginning of masses
     */
    void getPRMMasses(string sequence, vector<float> &masses, float offset =
        0.0, bool addZeroMass = false);

    /** AAJumps::getPRMMasses
     * Calculates PRM masses for input annotation in SpecNets format
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param spec:   Spectrum object with vector of masses to be modified.
     * @param offset: offset for prefix masses
     * @param addZeroMass: if true add offset mass to beginning of masses
     */
    void getPRMMasses(string sequence, Spectrum &spec, float offset=0.0, bool addZeroMass = false);

    /** AAJumps::getSRMMasses
     * Calculates SRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param masses: 	Vector of masses to be modified.
     * @param offset: offset for suffix masses
     * @param addZeroMass: if true add offset mass to beginning of masses
     */
    void getSRMMasses(string sequence, vector<float> &masses, float offset =
        0.0, bool addZeroMass = false);

    /** AAJumps::getPRMandSRMMasses
     * Calculates SRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     * @param prm_masses: 	Vector of prefix masses to be modified.
     * @param srm_masses: 	Vector of suffix masses to be modified.
     * @param addZeroMass: if true add offset mass to beginning of masses
     *
     */
    void getPRMandSRMMasses(string &sequence,
                            vector<float> &prm_masses,
                            vector<float> &srm_masses,
                            float &peptide_mass,
                            bool addZeroMass = false);

    /** AAJumps::getPeptideFromSpectrum
     * Converts a spectrum to a peptide sequence by instantiating consecutive mass differences
     *   to amino acids; mass differences are reported as "[mass]" whenever indistinguishable between amino acids (e.g., I/L)
     * @param spec: spectrum with the sequence represented as a series of cumulative prefix masses
     * @param sequence: output sequence
     * @param peakTol: mass error tolerance when deciding when a mass difference matches an amino acid (defaults to 0.45 Da)
     * @param offset: offset for prefix masses from theoretical masses (set to 0.0 for PRMs)
     */
    void getPeptideFromSpectrum(const Spectrum &spec,
                                string &sequence,
                                float peakTol = 0.45,
                                float offset = 0.0);

    /** AAJumps::getPeptideMass
     * Calculates sum of all aa masses for input annotation in SpecNets format *.SEQ[-17]UENCE.* NOT PARENT MASS
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     *
     */
    double getPeptideMass(string &sequence);

    /** AAJumps::getPeptideLength
     * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
     */
    int getPeptideLength(string &sequence);

  private:

    void init(short maxJumpSize,
              float resolution = 0.01,
              float peakTol = -1,
              short useMods = NO_MODS,
              bool uniqueMasses = true);
  };

  /**
   * Finds the longest amino acid tag contained in a sequence of masses
   * @param masses sequence of cummulative masses
   * @param jumps AAJumps class
   * return longest amino acid found by checking consecutive mass differences
   */
  string getLongestTag(vector<MZRange>& masses, AAJumps& jumps);
}
#endif
