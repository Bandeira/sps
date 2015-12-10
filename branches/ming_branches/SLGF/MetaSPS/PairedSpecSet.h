/*
 * PairedSpecSet.h
 *
 *  Created on: Nov 11, 2011
 *      Author: aguthals
 *
 */

#ifndef PAIREDSPECSET_H_
#define PAIREDSPECSET_H_

#include "spectrum.h"
#include "SpecSet.h"

using namespace std;

namespace specnets {

class PairedSpecSet {
public:

	/**
	 * All allocated heap memory for this class is automatically freed upon deallocation of this class.
	 */

	// Set of merged spectra. This is parallel to CID, HCD, and/or ETD spectra.
	SpecSet* mergedSet;

	// Labels for peaks in mergedSet that were derived from multiple complementary peaks in CID, HCD, and ETD.
	// 0 - Unlabeled, 1 - PRM, 2 - SRM, 3 - End point
	vector<vector<short> >* mergedLabels;

	// CID, HCD, and ETD PRM spectra
	SpecSet* spectraSetCID;
	SpecSet* spectraSetHCD;
	SpecSet* spectraSetETD;

	// Labels for peaks in CID, HCD, and ETD spectra that were merged. Keeping track of these makes sure peak
	//  intensities are not added to merged spectra more than once.
	// 0 - Unlabeled, 1 - PRM, 2 - SRM, 3 - End point
	vector<vector<short> >* usedPeaksCID;
	vector<vector<short> >* usedPeaksHCD;
	vector<vector<short> >* usedPeaksETD;

	// number of spectra in spectraSetCID, spectraSetHCD, spectraSetETD, and mergedSet
	int numSpecs;

	/**
	 * Outputs a spectrum and its labels to the output stream.
	 */
	static void outputLabels(ostream &output, SpecSet* spectra, vector<vector<
			short> >* usedPeaks, int specIdx);

	/**
	 * Default constructor. If this is called, the user must allocate memory for
	 *   spectraSetCID, spectraSetHCD and/or spectraSetETD before calling initialize()
	 */
	PairedSpecSet();

	PairedSpecSet(SpecSet* sequentialPairs);

	/**
	 * Constructor used for CID/ETD or HCD/ETD input pairs.
	 * @param input1 pointer to CID, HCD, or ETD input spectra, which will be copied into this class
	 * @param input2 pointer to CID, HCD, or ETD input spectra, which will be copied into this class
	 * @param frag1 CID, HCD, or ETD label for input1
	 * @param frag2 CID, HCD, or ETD label for input2
	 */
	PairedSpecSet(SpecSet* input1, SpecSet* input2, Spectrum::FragType& frag1,
			Spectrum::FragType& frag2);

	/**
	 * Constructor used for CID/HCD/ETD triples.
	 * @param inputCID
	 * @param inputHCD
	 * @param inputETD
	 */
	PairedSpecSet(SpecSet* inputCID, SpecSet* inputHCD, SpecSet* inputETD);

	~PairedSpecSet(void);

	/**
	 * Called in each constructor to initialize data structures. Call this to
	 *   re-initialize if CID, HCD, or ETD spectra have been modified.
	 */
	void initialize();

	/**
	 * Finds peaks with same mass in CID/ETD pairs or HCD/ETD pairs and translates
	 *   them and their complementary SRMs to the PRM mass in their merged spectrum with summed PRM+SRM
	 *   intensities (at a singular PRM peak mass)
	 */
	void mergePRMsStage1(float maxPeakDensity = -1.0);

	/**
	 * Finds peaks with mass difference 15+18 in CID/ETD pairs or HCD/ETD pairs AND where
	 *   at least one spectrum (CID, HCD or ETD) contains a PRM for these SRMs, translates them and
	 *   their complementary PRMs to the PRM mass in their merged spectrum with summed PRM+SRM
	 *   intensities
	 */
	void mergePRMsStage2(float maxPeakDensity = -1.0);

	/**
	 * Finds PRM / SRM pairs in CID/ETD pairs or HCD/ETD pairs and translates them to the PRM mass
	 *   in their merged spectrum with summed PRM+SRM intensities
	 */
	void mergePRMsStage3(float maxPeakDensity = -1.0);

	/**
	 * Finds peaks with mass difference 15+18 in CID/ETD pairs or HCD/ETD pairs and translates them and
	 *   their complementary PRMs to the PRM mass in their merged spectrum with summed PRM+SRM
	 *   intensities
	 */
	void mergePRMsStage4(float maxPeakDensity = -1.0);

	/**
	 * Copies all un-labeled left-over peaks from CID, HCD, and/or ETD spectra to their merged spectra
	 */
	void mergePRMsStage5(float maxPeakDensity = -1.0);

protected:

	void
	initializeUsedPeaks(vector<vector<short> >* usedPeaks, SpecSet* spectra);

	void moveSamePRMs(float maxPeakDensity, Spectrum::FragType fragCIDHCD);

	void moveSameSRMs(float maxPeakDensity, bool checkPRM,
			Spectrum::FragType fragCIDHCD);

	void movePRMsSRMs(float maxPeakDensity, Spectrum::FragType fragCIDHCD);

	void insertMergedLabel(vector<short>* peakLabels, int index, short label);

	void moveSamePRMsPair(float maxPeakDensity, int mergedIdx, int CIDIdx,
			int ETDIdx, Spectrum::FragType fragCIDHCD);

	void moveSameSRMsPair(float maxPeakDensity, int mergedIdx, int CIDIdx,
			int ETDIdx, bool checkPRM, Spectrum::FragType fragCIDHCD,
			list<int>* CIDscheckPRM = 0, list<int>* HCDscheckPRM = 0,
			list<int>* ETDscheckPRM = 0);

	void movePRMsSRMsPair(float maxPeakDensity, int mergedIdx, int CIDIdx,
			int ETDIdx, Spectrum::FragType fragCIDHCD);

	void moveLeftoversPair(float maxPeakDensity, int mergedIdx, int childIdx,
			Spectrum::FragType fragCIDHCDETD);

	void addNewPeaks(float maxPeakDensity, list<MZRange>* newPeaks,
			int mergedIdx, short label);

};
}

#endif /* PAIREDSPECSET_H_ */
