#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <cstring>

#include "vector.h"
#include "aminoacid.h"
#include "utils.h"
#include "PeptideSpectrumMatch.h"
#include "PeptideSpectrumMatchSet.h"

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

namespace specnets {
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
 * @see PeptideSpectrumMatch.h
 */
class PeptideSpectrumMatch;

/**
 * TODO: add description
 */
typedef std::tr1::shared_ptr<PeptideSpectrumMatch> psmPtr;

class Spectrum {

	/**
	 * TODO: add description
	 *
	 *@param minIdx
	 *@param maxIdx
	 *@param avgMass
	 *@param totScore
	 */
	void weightedAverage(int minIdx, int maxIdx, float &avgMass,
			float &totScore);

	/**
	 * The internal copy method used by copy ctor and operator=
	 *
	 *@param other  spectrum to copy
	 */
	void internalCopy(const Spectrum &other);

public:

	enum FragType {
		FragType_CID = 0, FragType_ETD = 1, FragType_HCD = 2
	};

	/**
	 * Associated pointers to annotations having to do with this spectrum
	 */
	list<psmPtr> psmList;

	/**
	 * Parent mass tolerance
	 */
	float parentMassTol;

	float oldParentMassTol;

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
	 * Below includes fields specifically required for Uploading
	 * Spectral Libraries onto CCMS
	 */

	/**
	 * Instrument Name
	 */
	string instrument_name;

	/**
	 * Ion Mass Tolerance
	 */
	float ITOL;

	/**
	 * Ion Mass Tolerance Units
	 */
	string ITOLU;

	/**
	 * Precursor Tolerance
	 */
	float TOL;

	/**
	 * Precusor Mass Tolerance Units
	 */
	string TOLU;

	/**
	 * Spectrum Quality, 1-5, 1 = low, 5 = high
	 */
	int spectrum_quality;

	/**
	 * Inter-Dalton distance: number of integer units between peaks that
	 * are 1 Da apart (idDist=1/resolution), e.g. 1 for regular spectra and
	 * 1/RESOLUTION for spectra with masses=round(mass/RESOLUTION),
	 * e.g. idDist=10 for RESOLUTION=0.1
	 */
	float idDist;

	/**
	 * Default constructor
	 */
	Spectrum();

	/**
	 * Copy constructor
	 */
	Spectrum(const Spectrum &other);

	/**
	 * Deconstructor
	 */
	~Spectrum();

	/**
	 * Initializes all class data structures
	 */
	void initialize();

	/**
	 * TODO: add description
	 *
	 *@param i
	 *@return
	 */
	TwoValues<float> &operator[](unsigned int i) {
		return peakList[i];
	}

	const TwoValues<float> &operator[](unsigned int i) const {
		return peakList[i];
	}

	/**
	 * TODO: add description
	 *
	 *@return
	 */
	unsigned int size() const {
		return peakList.size();
	}

	/**
	 * Resizes current peak list to another given size.
	 *
	 *@param newSize the new size for the peakList.
	 *@return the new size after the size has been set.
	 */
	unsigned int resize(unsigned int newSize) {
		peakList.resize(newSize);
		peakTols.resize(newSize);
		return peakList.size();
	}

	/**
	 * Calls reserve on peakList and peakTols
	 */
	void reserve(int newSize) {
		peakList.reserve(newSize);
		peakTols.reserve(newSize);
	}

	/**
	 * Spectrum assignment
	 *
	 * This will copy any associated annotations with the spectrum
	 *
	 *@param other
	 *@return
	 */
	Spectrum &operator=(const Spectrum &other);

	/**
	 * Copy all member fields except peakList and psmList.
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
	void setResolution(float newResolution, bool enforceRounding);

	/**
	 * Sets the tolerance of each peak in peakTols
	 */
	void setPeakTolerance(float peakTol);

	/**
	 * Sets the ppm tolerance of each peak in peakTols
	 */
	void setPeakTolerancePPM(float peakTolPPM);

	/**
	 * Sets the peak tolerance of the parent mass
	 */
	void setParentMassTol(float pmTol);

	/**
	 * Sets the ppm tolerance of the parent mass
	 */
	void setParentMassTolPPM(float pmTolPPM);

	/**
	 * Copies all per-peak tolerances into an internal vector
	 */
	void rememberTolerances();

	/**
	 * Restores peak tolerances to as they were before last call to rememberTolerances()
	 */
	void revertTolerances();

	/**
	 * Returns a pointer to the first peak in the Spectrum
	 */
	TwoValues<float>* front();

	/**
	 * Returns a pointer to the last peak in the Spectrum
	 */
	TwoValues<float>* back();

	/**
	 * Returns the peak tolerance for the peak at index
	 */
	float getTolerance(int index);

	/**
	 * Sets the peak tolerance for the peak at index
	 */
	void setTolerance(int index, float tolerance);

	/**
	 * Sorts peaks by ascending order of mass without invalidating peak tolerance indices
	 */
	void sortPeaks();

	/**
	 * Inserts a new peak into the spectrum.
	 * @param mass new peak's mass
	 * @param intensity new peak's intensity
	 * @param tolerance new peak's tolerance
	 * @param atIndex if specified, insert the new peak at this index. Otherwise
	 *  find an appropriate index such that the peak list stays in sorted order.
	 * @return index the of the new peak
	 */
	int insertPeak(float mass, float intensity, float tolerance, int atIndex =
			-1);

	/**
	 * Inserts a new peak into the spectrum.
	 * @param peak
	 * @param atIndex if specified, insert the new peak at this index. Otherwise
	 *  find an appropriate index such that the peak list stays in sorted order.
	 * @return index the of the new peak
	 */
	int insertPeak(MZRange* peak, int atIndex = -1);

	/**
	 * Inserts multiple peaks into the spectrum and keeps peaks in sorted order
	 *   (faster than multiple calls to insertPeak)
	 * @param newPeaks peaks to insert (each has mass, intensity, and tolerance)
	 * @return the new size of the spectrum
	 */
	int insertPeaks(list<MZRange>& newPeaks);

	/**
	 * Determines if a peak overlaps with peaks in this spectrum and finds the closest match along with all matches
	 * @param mass query peak's mass
	 * @param tolerance query peak's tolerance. mass is evaluated as mass +/- tolerance
	 * @param matches pointer to output list. Ignored if 0. Otherwise, this will hold the indicies of all peaks that
	 *   overlap with the query peak
	 * @param startIdx if != NULL, start looking here and avoid binary search to find closest mass
	 * @return index of closest peak to mass +/- tolerance if the ranges overlap. Otherwise, returns -1.
	 */
	int findPeaks(float mass, float tolerance = 0, list<int>* matches = 0, int* startIdx = 0);

	/**
	 * Determines if a peak overlaps with peaks in this spectrum and finds all matches
	 * @param range
	 * @param matches pointer to output list. Ignored if 0. Otherwise, this will hold the indices of all peaks that
	 *   overlap with the query peak
	 * @param startIdx if != NULL, start looking here and avoid binary search to find closest mass
	 * @return index of closest peak to range->getMass() if the ranges overlap. Otherwise, returns -1.
	 */
	int findPeaks(const MZRange& range, list<int>* matches = 0, int* startIdx = 0);

	/**
	 * Removes a peak from the spectrum. This fails if peakIndex is out of range
	 * @param peakIndex valid peak index of spectrum.
	 * @return true if peak was removed, false if not
	 */
	bool removePeak(int peakIndex);

	/**
	 * Removes multiple peaks from the spectrum (faster than multiple calls to removePeak()).
	 *   This fails if any peak indices are out of range
	 * @param peakIndices list of valid peak indices to be removed.
	 * @return true if all peaks were removed, false if not
	 */
	bool removePeaks(list<int>& peakIndices);

	/**
	 * Returns the spectrum's peak density as num peaks / parent mass
	 */
	float getPeakDensity();

	/**
	 * Add peaks for b0/bk/y0/yk
	 *
	 *@param tolerance
	 *@param ionOffset
	 *@param includeY0k
	 *@param ctermH2O
	 *@param labels
	 */
	void addZPMpeaks(float tolerance, float ionOffset, bool includeY0k,
			bool ctermH2O = true, SpectrumPeakLabels *labels = 0);

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
	 * Normalizes to Euclidian norm (norm 2)
	 *
	 *@param newNorm2
	 */
	void normalize2(float newNorm2 = 1.0);

	/**
	 * TODO: add description
	 *
	 *@param parentMS1
	 *@param peakTol
	 *@param maxZ
	 *@param isoEnvs
	 *@param strictMode
	 */
	void guessPrecursorZPM(Spectrum &parentMS1, float peakTol, short maxZ,
			IsoEnvelope &isoEnvs, bool strictMode = false);

	/**
	 * Find all peaks within tolerance of baseMass.
	 *
	 * If startIdx is >=0 then the search starts at startIdx. Function returns the number of peaks within tolerance.
	 *
	 *@param baseMass
	 *@param peakTol
	 *@param matchesIdx
	 *@return
	 */
	short findMatches(float baseMass, float peakTol, vector<int> &matchesIdx,
			int startIdx = -1);

	/**
	 * Finds the spectrum peak with peak mass closest to mass
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
	 * Computes the reversed spectrum and returns it in putHere (or reverses current if putHere==NULL)
	 *   The reverse of a peak with mass m is parentMass-1+pmOffset-m :
	 *   use pmOffset=0 for PRM spectra and pmOffset=2*AAJumps::massHion for MS/MS spectra.
	 *
	 *@param pmOffset
	 *@param putHere
	 */
	void reverse(float pmOffset, Spectrum *putHere = 0);

	/**
	 * Merges this spectrum's peakList with newSpec. Peaks in newSpec are merged with
	 * their closest mass in this spectrum if both peak masses are within peak tolerance of each other.
	 *
	 *@param newSpec Spectrum to merge with with this spectrum
	 *@param mergeType
	 *           0: the mass of each merged peak is the weighted average of merged masses
	 *           1: the mass of each merged peak is that of the merged peak with the highest intensity
	 *           2: the mass of each merged peak is that of its merged peak in this spectrum
	 *           3: the mass of each merged peak is that of its merged peak with the lowest tolerance
	 */
	void mergeClosestPeaks(Spectrum& newSpec, short mergeType);

	/**
	 * Adds offset to every peak mass and rotates negative masses to parentMass-cterm+mass
	 *
	 *@param offset
	 *@param cterm
	 */
	void rotate(float offset, float cterm);

	/**
	 * Selects the top k peaks in the spectrum.
	 *
	 *@param topK
	 *@param putHere
	 */
	void selectTopK(unsigned int topK, Spectrum *putHere = 0);

	/**
	 * A peak must be in the top k in a window of [w(0),w(1)] around its mass to be retained in the spectrum.
	 *
	 *@param topK
	 *@param putHere
	 */
	void selectTopK(unsigned int topK, TwoValues<float> w, Spectrum *putHere =
			0);

	/**
	 * Converts a list of masses to a list of corresponding spectrum peak
	 * indices (closest mass). Masses is assumed to be sorted by increasing mass values.
	 *
	 *@param masses
	 *@param indices
	 *@param peakTol
	 */
	void massesToIndices(Spectrum &masses, vector<int> &indices, float peakTol);

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
	unsigned int filterLowIntensity(float minIntensity = 0);

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
	 * Computes the sets of pairs in the spectrum and returns (in pairs) a list of
	 * pairs sorted by minimum distance to the closest endpoint. Two peaks are
	 * considered a pair if their masses add up to parentMass-idDist+pmOffset ->
	 * use pmOffset=0 for PRM spectra and pmOffset=2 for MS/MS spectra.
	 *
	 *@param pmOffset
	 *@param tolerance
	 *@param pairs - returns the masses of the peaks in the pairs (per pair)
	 *@param pairsIdx - returns the indices of the peaks in the pairs (per pair)
	 */
	void getSymmetricPeakPairs(float pmOffset, float tolerance,
			vector<vector<float> > &pairs, vector<vector<int> > &pairsIdx);

	/**
	 * Removes peaks that have an intensity ranked worse than maxRank
	 *   compared to all neighboring peaks +/- windowRadius
	 * @param maxRank maximum allowable rank of each peak
	 * @param windowRadius radius of peak comparison in the spectrum
	 * @return
	 */
	void rankFilterPeaks(int maxRank, float windowRadius = AAJumps::minAAmass
			- 1.0);
	/**
	 * Forces a spectrum to be symmetric by adding symmetric peaks
	 * whenever missing from the spectrum.
	 *
	 *@param pmOffset - use 0 for PRM spectra, 2 for MS/MS spectra
	 *@param tolerance
	 *@param indices if indices!=NULL then each position contains the original index of the entry (or -1 if it was added by this function)
	 */
	void makeSymmetric(float pmOffset, float tolerance, vector<int>* indicies = 0);

	// Generic functions but initially defined for alignment of consensus to PRM spectra (batch.cpp::getPairAlignsPAext_aux)
	/**
	 * TODO: add description
	 *
	 *@param toSpec
	 *@return
	 */
	bool compare(Spectrum &toSpec);

	/**
	 * mergeCommon - merge the common peaks in the current spectrum with those in
	 * withSpec and output the result to toSpec. Multiple peak matches (within tolerance)
	 * are represented by different peaks in toSpec. Each peak in toSpec has its mass
	 * given by the weighted average mass of the matched peaks. The intensity of
	 * a merged peak is determined by the intensities of the paired peaks acording
	 * to mergeType: 0 (sum), 1 (max), 2 (min), 3 (cur spec score only), 4 (withSpectrum score only)
	 *
	 *@param withSpectrum - 0 (sum), 1 (max), 2 (min), 3 (cur spec score only), 4 (withSpectrum score only)
	 *@param toSpec
	 *@param shift
	 *@param peakTol
	 *@param mergeType
	 */
	void mergeCommon(Spectrum &withSpectrum, Spectrum *toSpec, float shift,
			float peakTol, short mergeType = 0);

	bool saveDTA(const char* outfile);

	/**
	 * Get total ion current (total intensity)
	 */
	float getTotalIonCurrent(void);

protected:
	/**
	 * TODO: add description
	 */
	vector<TwoValues<float> > peakList;

	/**
	 * vector of peak tolerances, parallel to peakList
	 */
	vector<float> peakTols;

	/**
	 * backup vector of tolerances. This is used to revert peak tolerances as older functions are called with the peakTol parameter
	 */
	vector<float> oldPeakTols;
};

}
#endif
