/*
 * PairedSpecSet.cpp
 *
 *  Created on: Nov 11, 2011
 *      Author: aguthals
 */

#include "PairedSpecSet.h"
#include "math.h"
#include "Logger.h"

using namespace std;

namespace specnets {

struct SortPeakIntensity: public std::binary_function<MZRange, MZRange, bool> {
	bool operator()(MZRange left, MZRange right) const {
		return left.getIntensity() > right.getIntensity();
	}
	;
};

void PairedSpecSet::outputLabels(ostream &output, SpecSet* spectra, vector<
		vector<short> >* usedPeaks, int specIdx) {

	if (spectra->size() != usedPeaks->size()) {
		ERROR_MSG("Used peaks are different size than SpecSet!");
		return;
	} else if (specIdx < 0 || specIdx >= spectra->size()) {
		ERROR_MSG("Invalid spectrum index " << specIdx << " for SpecSet of size " << spectra->size());
		return;
	} else if ((*spectra)[specIdx].size() != (*usedPeaks)[specIdx].size()) {
		ERROR_MSG("Used peaks are different size than spectrum at index " << specIdx);
		DEBUG_VAR((*spectra)[specIdx].size());
		DEBUG_VAR((*usedPeaks)[specIdx].size());
		return;
	}
	output << (*spectra)[specIdx].parentMass << " "
			<< (*spectra)[specIdx].parentCharge << endl;
	for (unsigned int i = 0; i < (*spectra)[specIdx].size(); i++) {
		output << (*spectra)[specIdx][i][0] << " " << (*spectra)[specIdx][i][1]
				<< " " << (*spectra)[specIdx].getTolerance(i) << " "
				<< (*usedPeaks)[specIdx][i] << endl;
	}
	output << endl;
}

PairedSpecSet::PairedSpecSet() :
	mergedSet(0x0), mergedLabels(0x0), spectraSetCID(0x0), spectraSetHCD(0x0),
			spectraSetETD(0x0), usedPeaksCID(0x0), usedPeaksHCD(0x0),
			usedPeaksETD(0x0), numSpecs(0) {
	initialize();
}

PairedSpecSet::PairedSpecSet(SpecSet* sequentialPairs) :
	mergedSet(0x0), mergedLabels(0x0), spectraSetCID(0x0), spectraSetHCD(0x0),
			spectraSetETD(0x0), usedPeaksCID(0x0), usedPeaksHCD(0x0),
			usedPeaksETD(0x0), numSpecs(0) {

	bool cid = false, hcd = false, etd = false;
	int numPaired = 0;
	for (int i = 0; i < 3; i++) {
		if ((*sequentialPairs)[i].msFragType == Spectrum::FragType_CID && !cid) {
			cid = true;
			++numPaired;
		}
		if ((*sequentialPairs)[i].msFragType == Spectrum::FragType_HCD && !hcd) {
			hcd = true;
			++numPaired;
		}
		if ((*sequentialPairs)[i].msFragType == Spectrum::FragType_ETD && !etd) {
			etd = true;
			++numPaired;
		}
	}
	DEBUG_VAR(numPaired);
	int numParallel = sequentialPairs->size() / numPaired;
	if (cid) {
		spectraSetCID = new SpecSet(numParallel);
	}
	if (hcd) {
		spectraSetHCD = new SpecSet(numParallel);
	}
	if (etd) {
		spectraSetETD = new SpecSet(numParallel);
	}

	for (int i = 0; i < sequentialPairs->size(); i++) {
		if ((*sequentialPairs)[i].msFragType == Spectrum::FragType_CID) {
			(*spectraSetCID)[i / numPaired] = (*sequentialPairs)[i];
		} else if ((*sequentialPairs)[i].msFragType == Spectrum::FragType_HCD) {
			(*spectraSetHCD)[i / numPaired] = (*sequentialPairs)[i];
		} else {
			(*spectraSetETD)[i / numPaired] = (*sequentialPairs)[i];
		}
	}

	initialize();
}

PairedSpecSet::PairedSpecSet(SpecSet* input1, SpecSet* input2,
		Spectrum::FragType& frag1, Spectrum::FragType& frag2) :
	mergedSet(0x0), mergedLabels(0x0), spectraSetCID(0x0), spectraSetHCD(0x0),
			spectraSetETD(0x0), usedPeaksCID(0x0), usedPeaksHCD(0x0),
			usedPeaksETD(0x0), numSpecs(0) {
	if (frag1 == frag2) {
		ERROR_MSG("Must have separate fragmentation methods for input spectra");
		throw -1;
	}

	if (frag1 == Spectrum::FragType_CID) {
		spectraSetCID = new SpecSet(input1->size());
		*spectraSetCID = *input1;
	} else if (frag1 == Spectrum::FragType_HCD) {
		spectraSetHCD = new SpecSet(input1->size());
		*spectraSetHCD = *input1;
	} else if (frag1 == Spectrum::FragType_ETD) {
		spectraSetETD = new SpecSet(input1->size());
		*spectraSetETD = *input1;
	} else {
		ERROR_MSG("Found unsupported fragmentation identifier " << frag1 << " for first set of spectra");
		throw -1;
	}

	if (frag2 == Spectrum::FragType_CID) {
		spectraSetCID = new SpecSet(input2->size());
		*spectraSetCID = *input2;
	} else if (frag2 == Spectrum::FragType_HCD) {
		spectraSetHCD = new SpecSet(input2->size());
		*spectraSetHCD = *input2;
	} else if (frag2 == Spectrum::FragType_ETD) {
		spectraSetETD = new SpecSet(input2->size());
		*spectraSetETD = *input2;
	} else {
		ERROR_MSG("Found unsupported fragmentation identifier " << frag2 << " for second set of spectra");
		throw -1;
	}

	initialize();
}

PairedSpecSet::PairedSpecSet(SpecSet* inputCID, SpecSet* inputHCD,
		SpecSet* inputETD) :
	mergedSet(0x0), mergedLabels(0x0), spectraSetCID(0x0), spectraSetHCD(0x0),
			spectraSetETD(0x0), usedPeaksCID(0x0), usedPeaksHCD(0x0),
			usedPeaksETD(0x0), numSpecs(0) {
	spectraSetCID = new SpecSet(inputCID->size());
	*spectraSetCID = *inputCID;

	spectraSetHCD = new SpecSet(inputHCD->size());
	*spectraSetHCD = *inputHCD;

	spectraSetETD = new SpecSet(inputETD->size());
	*spectraSetETD = *inputETD;
	initialize();
}

PairedSpecSet::~PairedSpecSet() {
	if (mergedSet != 0) {
		delete mergedSet;
	}

	if (mergedLabels != 0) {
		delete mergedLabels;
	}

	if (spectraSetCID != 0) {
		delete spectraSetCID;
	}

	if (spectraSetHCD != 0) {
		delete spectraSetHCD;
	}

	if (spectraSetETD != 0) {
		delete spectraSetETD;
	}

	if (usedPeaksCID != 0) {
		delete usedPeaksCID;
	}

	if (usedPeaksHCD != 0) {
		delete usedPeaksHCD;
	}

	if (usedPeaksETD != 0) {
		delete usedPeaksETD;
	}
}

void PairedSpecSet::initialize() {
	numSpecs = -1;
	int numCID = -1;
	int numHCD = -1;
	int numETD = -1;
	SpecSet* templateSpecs;
	int numPaired = 0;

	if (spectraSetCID != 0) {
		numSpecs = spectraSetCID->size();
		numCID = spectraSetCID->size();
		templateSpecs = spectraSetCID;
		numPaired++;
	}
	if (spectraSetHCD != 0) {
		numSpecs = spectraSetHCD->size();
		numHCD = spectraSetHCD->size();
		templateSpecs = spectraSetHCD;
		numPaired++;
	}
	if (spectraSetETD != 0) {
		numSpecs = spectraSetETD->size();
		numETD = spectraSetETD->size();
		templateSpecs = spectraSetETD;
		numPaired++;
	}
	if (numSpecs < 0) {
		ERROR_MSG("No spectra have been loaded, initialization of PairedSpecSet failed!");
		throw -1;
	}

	if (numCID >= 0 && numHCD >= 0 && numCID != numHCD) {
		ERROR_MSG("Cannot pair " << numCID << " CID spectra with " << numHCD << " HCD spectra, both must be equal");
		throw -1;
	}
	if (numCID >= 0 && numETD >= 0 && numCID != numETD) {
		ERROR_MSG("Cannot pair " << numCID << " CID spectra with " << numETD << " ETD spectra, both must be equal");
		throw -1;
	}

	if (numHCD >= 0 && numETD >= 0 && numHCD != numETD) {
		ERROR_MSG("Cannot pair " << numHCD << " HCD spectra with " << numETD << " ETD spectra, both must be equal");
		throw -1;
	}

	if (mergedSet == 0) {
		mergedSet = new SpecSet(numSpecs);
	} else {
		mergedSet->resize(numSpecs);
	}

	if (mergedLabels == 0) {
		mergedLabels = new vector<vector<short> > (numSpecs);
	} else {
		mergedLabels->resize(numSpecs);
	}

	if (spectraSetCID != 0) {
		if (usedPeaksCID == 0) {
			usedPeaksCID = new vector<vector<short> > (numSpecs);
		}
		initializeUsedPeaks(usedPeaksCID, spectraSetCID);
	} else {
		if (usedPeaksCID != 0) {
			delete usedPeaksCID;
		}
	}

	if (spectraSetHCD != 0) {
		if (usedPeaksHCD == 0) {
			usedPeaksHCD = new vector<vector<short> > (numSpecs);
		}
		initializeUsedPeaks(usedPeaksHCD, spectraSetHCD);
	} else {
		if (usedPeaksHCD != 0) {
			delete usedPeaksHCD;
		}
	}

	if (spectraSetETD != 0) {
		if (usedPeaksETD == 0) {
			usedPeaksETD = new vector<vector<short> > (numSpecs);
		}
		initializeUsedPeaks(usedPeaksETD, spectraSetETD);
	} else {
		if (usedPeaksETD != 0) {
			delete usedPeaksETD;
		}
	}

	float totalScore = 0, numScores = 0, averageScore = 0, standardDev = 0,
			endPtScore = 0;
	list<float> scores;
	for (int i = 0; i < numSpecs; i++) {
		MZRange massB0(0, 0, 0.0001);
		MZRange massBk((*templateSpecs)[i].parentMass - AAJumps::massMH, 0,
				(*templateSpecs)[i].parentMassTol);
		MZRange massY0(AAJumps::massH2O, 0, 0.0001);
		MZRange massYk((*templateSpecs)[i].parentMass - AAJumps::massHion, 0,
				(*templateSpecs)[i].parentMassTol);
		int idxFound;

		if (spectraSetCID != 0) {
			idxFound = (*spectraSetCID)[i].findPeaks(massB0);
			if (idxFound >= 0) {
				(*usedPeaksCID)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetCID)[i].findPeaks(massBk);
			if (idxFound >= 0) {
				(*usedPeaksCID)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetCID)[i].findPeaks(massY0);
			if (idxFound >= 0) {
				(*usedPeaksCID)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetCID)[i].findPeaks(massYk);
			if (idxFound >= 0) {
				(*usedPeaksCID)[i][idxFound] = 3;
			}
			for (int j = 0; j < (*spectraSetCID)[i].size(); j++) {
				totalScore += (*spectraSetCID)[i][j][1];
				scores.push_back((*spectraSetCID)[i][j][1]);
				numScores += 1.0;
			}
		}

		if (spectraSetHCD != 0) {
			idxFound = (*spectraSetHCD)[i].findPeaks(massB0);
			if (idxFound >= 0) {
				(*usedPeaksHCD)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetHCD)[i].findPeaks(massBk);
			if (idxFound >= 0) {
				(*usedPeaksHCD)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetHCD)[i].findPeaks(massY0);
			if (idxFound >= 0) {
				(*usedPeaksHCD)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetHCD)[i].findPeaks(massYk);
			if (idxFound >= 0) {
				(*usedPeaksHCD)[i][idxFound] = 3;
			}
			for (int j = 0; j < (*spectraSetHCD)[i].size(); j++) {
				totalScore += (*spectraSetHCD)[i][j][1];
				scores.push_back((*spectraSetHCD)[i][j][1]);
				numScores += 1.0;
			}
		}

		MZRange massZk((*templateSpecs)[i].parentMass - AAJumps::massMH
				- AAJumps::massNH, 0, (*templateSpecs)[i].parentMassTol);

		if (spectraSetETD != 0) {
			idxFound = (*spectraSetETD)[i].findPeaks(massB0);
			if (idxFound >= 0) {
				(*usedPeaksETD)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetETD)[i].findPeaks(massBk);
			if (idxFound >= 0) {
				(*usedPeaksETD)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetETD)[i].findPeaks(massY0);
			if (idxFound >= 0) {
				(*usedPeaksETD)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetETD)[i].findPeaks(massYk);
			if (idxFound >= 0) {
				(*usedPeaksETD)[i][idxFound] = 3;
			}
			idxFound = (*spectraSetETD)[i].findPeaks(massZk);
			if (idxFound >= 0) {
				(*usedPeaksETD)[i][idxFound] = 3;
			}
			for (int j = 0; j < (*spectraSetETD)[i].size(); j++) {
				totalScore += (*spectraSetETD)[i][j][1];
				scores.push_back((*spectraSetETD)[i][j][1]);
				numScores += 1.0;
			}
		}
	}
	averageScore = (totalScore / numScores);
	for (list<float>::iterator scoreIt = scores.begin(); scoreIt
			!= scores.end(); scoreIt++) {
		float diffScore = (*scoreIt) - averageScore;
		standardDev += diffScore * diffScore;
	}

	standardDev = sqrt(standardDev / numScores);

	endPtScore = (averageScore + (2.0 * standardDev)) * ((float) numPaired);
	int insertIdx;
	for (int i = 0; i < numSpecs; i++) {
		(*mergedSet)[i].copyNP((*templateSpecs)[i]);
		(*mergedSet)[i].resize(0);
		float massB0 = 0;
		float massBk = (*mergedSet)[i].parentMass - AAJumps::massMH;
		float massY0 = AAJumps::massH2O;
		float massYk = (*mergedSet)[i].parentMass - AAJumps::massHion;
		insertIdx = (*mergedSet)[i].insertPeak(massB0, endPtScore, 0.0001);
		insertIdx = (*mergedSet)[i].insertPeak(massBk, endPtScore,
				(*mergedSet)[i].parentMassTol);
		insertIdx = (*mergedSet)[i].insertPeak(massY0, endPtScore, 0.0001);
		insertIdx = (*mergedSet)[i].insertPeak(massYk, endPtScore,
				(*mergedSet)[i].parentMassTol);
		(*mergedLabels)[i].resize(0);
		(*mergedLabels)[i].resize(4, (short) 3);
	}
}

void PairedSpecSet::mergePRMsStage1(float maxPeakDensity) {

	if (spectraSetETD == 0) {
		WARN_MSG("Found no ETD spectra, skipping Stage 1 merge");
		return;
	}
	if (spectraSetCID != 0) {
		moveSamePRMs(maxPeakDensity, Spectrum::FragType_CID);
	} else {
		DEBUG_MSG("No CID spectra, skipping CID/ETD Stage 1 merge");
	}
	if (spectraSetHCD != 0) {
		moveSamePRMs(maxPeakDensity, Spectrum::FragType_HCD);
	} else {
		DEBUG_MSG("No HCD spectra, skipping HCD/ETD Stage 1 merge");
	}
}

void PairedSpecSet::mergePRMsStage2(float maxPeakDensity) {
	if (spectraSetETD == 0) {
		WARN_MSG("Found no ETD spectra, skipping Stage 2 merge");
		return;
	}
	if (spectraSetCID != 0) {
		moveSameSRMs(maxPeakDensity, true, Spectrum::FragType_CID);
	} else {
		DEBUG_MSG("No CID spectra, skipping CID/ETD Stage 2 merge");
	}
	if (spectraSetHCD != 0) {
		moveSameSRMs(maxPeakDensity, true, Spectrum::FragType_HCD);
	} else {
		DEBUG_MSG("No HCD spectra, skipping HCD/ETD Stage 2 merge");
	}
}

void PairedSpecSet::mergePRMsStage3(float maxPeakDensity) {
	if (spectraSetETD == 0) {
		WARN_MSG("Found no ETD spectra, skipping Stage 4 merge");
		return;
	}
	if (spectraSetCID != 0) {
		movePRMsSRMs(maxPeakDensity, Spectrum::FragType_CID);
	} else {
		DEBUG_MSG("No CID spectra, skipping CID/ETD Stage 4 merge");
	}
	if (spectraSetHCD != 0) {
		movePRMsSRMs(maxPeakDensity, Spectrum::FragType_HCD);
	} else {
		DEBUG_MSG("No HCD spectra, skipping HCD/ETD Stage 4 merge");
	}
}

void PairedSpecSet::mergePRMsStage4(float maxPeakDensity) {

	if (spectraSetETD == 0) {
		WARN_MSG("Found no ETD spectra, skipping Stage 3 merge");
		return;
	}
	if (spectraSetCID != 0) {
		moveSameSRMs(maxPeakDensity, false, Spectrum::FragType_CID);
	} else {
		DEBUG_MSG("No CID spectra, skipping CID/ETD Stage 3 merge");
	}
	if (spectraSetHCD != 0) {
		moveSameSRMs(maxPeakDensity, false, Spectrum::FragType_HCD);
	} else {
		DEBUG_MSG("No HCD spectra, skipping HCD/ETD Stage 3 merge");
	}
}

void PairedSpecSet::mergePRMsStage5(float maxPeakDensity) {

	bool addedETD, addedHCD;

	if (spectraSetCID == 0) {
		DEBUG_MSG("No CID spectra, skipping CID Stage 5 merge");
	}

	if (spectraSetHCD == 0) {
		DEBUG_MSG("No HCD spectra, skipping HCD Stage 5 merge");
	}

	if (spectraSetETD == 0) {
		DEBUG_MSG("No ETD spectra, skipping ETD Stage 5 merge");
	}

	for (int i = 0; i < numSpecs; i++) {
		addedETD = false;
		addedHCD = false;

		// Add ETD and HCD before CID for spectra with highly charged precursors
		if (spectraSetETD != 0 && (*mergedSet)[i].parentCharge > 3) {
			moveLeftoversPair(maxPeakDensity, i, i, Spectrum::FragType_ETD);
			addedETD = true;
		}

		if (spectraSetHCD != 0 && (*mergedSet)[i].parentCharge > 3) {
			moveLeftoversPair(maxPeakDensity, i, i, Spectrum::FragType_HCD);
			addedHCD = true;
		}

		if (spectraSetCID != 0) {
			moveLeftoversPair(maxPeakDensity, i, i, Spectrum::FragType_CID);
		}
		if (spectraSetHCD != 0 && !addedHCD) {
			moveLeftoversPair(maxPeakDensity, i, i, Spectrum::FragType_HCD);
		}
		if (spectraSetETD != 0 && !addedETD) {
			moveLeftoversPair(maxPeakDensity, i, i, Spectrum::FragType_ETD);
		}
	}
}

void PairedSpecSet::initializeUsedPeaks(vector<vector<short> >* usedPeaks,
		SpecSet* spectra) {
	for (int i = 0; i < spectra->size(); i++) {
		(*usedPeaks)[i].resize((*spectra)[i].size());
		for (int j = 0; j < (*usedPeaks)[i].size(); j++) {
			(*usedPeaks)[i][j] = 0;
		}
	}
}

void PairedSpecSet::moveSamePRMs(float maxPeakDensity,
		Spectrum::FragType fragCIDHCD) {
	for (int i = 0; i < numSpecs; i++) {
		moveSamePRMsPair(maxPeakDensity, i, i, i, fragCIDHCD);
	}
}

void PairedSpecSet::moveSameSRMs(float maxPeakDensity, bool checkPRM,
		Spectrum::FragType fragCIDHCD) {
	list<int> idxCheck;
	for (int i = 0; i < numSpecs; i++) {
		idxCheck.clear();
		idxCheck.push_back(i);
		moveSameSRMsPair(maxPeakDensity, i, i, i, checkPRM, fragCIDHCD,
				&idxCheck, &idxCheck, &idxCheck);
	}
}

void PairedSpecSet::movePRMsSRMs(float maxPeakDensity,
		Spectrum::FragType fragCIDHCD) {
	for (int i = 0; i < numSpecs; i++) {
		movePRMsSRMsPair(maxPeakDensity, i, i, i, fragCIDHCD);
	}
}

void PairedSpecSet::insertMergedLabel(vector<short>* peakLabels, int index,
		short label) {
	vector<short>::iterator labelIt = peakLabels->begin();
	labelIt += index;
	peakLabels->insert(labelIt, label);
}

void PairedSpecSet::moveSamePRMsPair(float maxPeakDensity, int mergedIdx,
		int CIDIdx, int ETDIdx, Spectrum::FragType fragCIDHCD) {
	SpecSet* CIDspecs;
	vector<vector<short> >* CIDused;

	if (fragCIDHCD == Spectrum::FragType_CID) {
		CIDspecs = spectraSetCID;
		CIDused = usedPeaksCID;
	} else {
		CIDspecs = spectraSetHCD;
		CIDused = usedPeaksHCD;
	}

	Spectrum* CIDspec = &(*CIDspecs)[CIDIdx];
	Spectrum* ETDspec = &(*spectraSetETD)[ETDIdx];
	Spectrum* mergedSpec = &(*mergedSet)[mergedIdx];

	if (mergedSpec->getPeakDensity() > maxPeakDensity) {
		return;
	}

	if (CIDspec->size() == 0 || ETDspec->size() == 0) {
		return;
	}
	/*
	 bool debug = false;
	 if (mergedIdx == 2207) {
	 debug = true;
	 }
	 */
	vector<short>* CIDlabels = &(*CIDused)[CIDIdx];
	vector<short>* ETDlabels = &(*usedPeaksETD)[ETDIdx];
	vector<short>* mergedLocLabels = &(*mergedLabels)[mergedIdx];

	MZRange rangeCID;
	MZRange rangeETD;
	MZRange cPeakCID;
	MZRange cPeakETD;
	MZRange mergeRange;
	list<MZRange> peaksToAdd;

	int numPeaks = CIDspec->size();

	float pmTol = mergedSpec->parentMassTol;

	float massBk = mergedSpec->parentMass - AAJumps::massMH;

	for (int j = 0; j < numPeaks; j++) {

		if ((*CIDlabels)[j] == 2 || (*CIDlabels)[j] == 3) {
			continue;
		}

		float massCID = (*CIDspec)[j][0];
		float intenCID = ((*CIDlabels)[j] == 0) ? (*CIDspec)[j][1] : 0;
		float massTolCID = CIDspec->getTolerance(j);
		rangeCID.set(massCID, intenCID, massTolCID);

		int closestETD = ETDspec->findClosest(massCID);
		float massETD = (*ETDspec)[closestETD][0];

		float intenETD =
				((*ETDlabels)[closestETD] == 0) ? (*ETDspec)[closestETD][1] : 0;
		float massTolETD = ETDspec->getTolerance(closestETD);
		rangeETD.set(massETD, intenETD, massTolETD);

		if (rangeCID != rangeETD || (*ETDlabels)[closestETD] == 2
				|| (*ETDlabels)[closestETD] == 3) {
			continue;
		}

		if ((*CIDlabels)[j] == 1 && (*ETDlabels)[closestETD] == 1) {
			continue;
		} else {
			(*CIDlabels)[j] = 1;
			(*ETDlabels)[closestETD] = 1;
		}

		mergeRange = (massTolCID < massTolETD) ? rangeCID : rangeETD;
		mergeRange.setIntensity(intenCID + intenETD);

		cPeakCID.set(massBk - mergeRange.getMass() + AAJumps::massH2O, 0, pmTol
				+ mergeRange.getTolerance());
		cPeakETD.set(massBk - mergeRange.getMass() - AAJumps::massNH, 0, pmTol
				+ mergeRange.getTolerance());

		int cPeakIdxCID = CIDspec->findPeaks(cPeakCID);
		int cPeakIdxETD = ETDspec->findPeaks(cPeakETD);

		if (cPeakIdxCID >= 0 && (*CIDlabels)[cPeakIdxCID] == 0) {
			mergeRange.setIntensity(mergeRange.getIntensity()
					+ (*CIDspec)[cPeakIdxCID][1]);
			(*CIDlabels)[cPeakIdxCID] = 2;
		}

		if (cPeakIdxETD >= 0 && (*ETDlabels)[cPeakIdxETD] == 0) {
			mergeRange.setIntensity(mergeRange.getIntensity()
					+ (*ETDspec)[cPeakIdxETD][1]);
			(*ETDlabels)[cPeakIdxETD] = 2;
		}

		peaksToAdd.push_back(mergeRange);
	}
	addNewPeaks(maxPeakDensity, &peaksToAdd, mergedIdx, 1);
	/*
	 if (debug) {
	 DEBUG_VAR(mergedIdx);
	 outputLabels(cout, mergedSet, mergedLabels, mergedIdx);
	 cout << "\nCID spectrum (" << fragCIDHCD << ")\n";
	 outputLabels(cout, CIDspecs, CIDused, CIDIdx);
	 cout << "\nETD spectrum:";
	 outputLabels(cout, spectraSetETD, usedPeaksETD, ETDIdx);
	 }
	 */
}

void PairedSpecSet::moveSameSRMsPair(float maxPeakDensity, int mergedIdx,
		int CIDIdx, int ETDIdx, bool checkPRM, Spectrum::FragType fragCIDHCD,
		list<int>* CIDscheckPRM, list<int>* HCDscheckPRM,
		list<int>* ETDscheckPRM) {
	SpecSet* CIDspecs;
	vector<vector<short> >* CIDused;

	if (fragCIDHCD == Spectrum::FragType_CID) {
		CIDspecs = spectraSetCID;
		CIDused = usedPeaksCID;
	} else {
		CIDspecs = spectraSetHCD;
		CIDused = usedPeaksHCD;
	}

	Spectrum* CIDspec = &(*CIDspecs)[CIDIdx];
	Spectrum* ETDspec = &(*spectraSetETD)[ETDIdx];
	Spectrum* mergedSpec = &(*mergedSet)[mergedIdx];

	if (mergedSpec->getPeakDensity() > maxPeakDensity) {
		return;
	}
	if (CIDspec->size() == 0 || ETDspec->size() == 0) {
		return;
	}

	vector<short>* CIDlabels = &(*CIDused)[CIDIdx];
	vector<short>* ETDlabels = &(*usedPeaksETD)[ETDIdx];
	vector<short>* mergedLocLabels = &(*mergedLabels)[mergedIdx];

	int numPeaks = CIDspec->size();

	float pmTol = mergedSpec->parentMassTol;

	float massBk = mergedSpec->parentMass - AAJumps::massMH;

	MZRange rangeETD;
	MZRange rangeCID;
	MZRange mergeRange;
	list<MZRange> peaksToAdd;

	for (int j = 0; j < numPeaks; j++) {

		if ((*CIDlabels)[j] == 1 || (*CIDlabels)[j] == 3) {
			continue;
		}

		float massCID = (*CIDspec)[j][0];
		float intenCID = ((*CIDlabels)[j] == 0) ? (*CIDspec)[j][1] : 0;
		float massTolCID = CIDspec->getTolerance(j);
		rangeCID.set(massCID - AAJumps::massH2O - AAJumps::massNH, intenCID,
				massTolCID);

		int closestETD = ETDspec->findClosest(rangeCID.getMass());

		float massETD = (*ETDspec)[closestETD][0];
		float intenETD =
				((*ETDlabels)[closestETD] == 0) ? (*ETDspec)[closestETD][1] : 0;
		float massTolETD = ETDspec->getTolerance(closestETD);
		rangeETD.set(massETD, intenETD, massTolETD);

		if (rangeCID != rangeETD || (*ETDlabels)[closestETD] == 1
				|| (*ETDlabels)[closestETD] == 3) {
			continue;
		}

		if ((*CIDlabels)[j] == 2 && (*ETDlabels)[closestETD] == 2) {
			continue;
		}

		mergeRange = (massTolCID < massTolETD) ? rangeCID : rangeETD;
		mergeRange.setMass(massBk - (mergeRange.getMass() + AAJumps::massNH));
		mergeRange.setTolerance(pmTol + mergeRange.getTolerance());
		mergeRange.setIntensity(intenCID + intenETD);

		int prmIdx = -1;
		if (checkPRM) {
			if (CIDscheckPRM != 0 && spectraSetCID != 0) {
				for (list<int>::iterator idxIt = CIDscheckPRM->begin(); idxIt
						!= CIDscheckPRM->end(); idxIt++) {
					int i = *idxIt;
					int prmIdxTemp = (*spectraSetCID)[i].findPeaks(mergeRange);
					if (prmIdxTemp >= 0 && (*usedPeaksCID)[i][prmIdxTemp] != 2
							&& (*usedPeaksCID)[i][prmIdxTemp] != 3) {
						prmIdx = prmIdxTemp;
						if ((*usedPeaksCID)[i][prmIdx] == 0) {
							mergeRange.setIntensity(mergeRange.getIntensity()
									+ (*spectraSetCID)[i][prmIdx][1]);
						}
						(*usedPeaksCID)[i][prmIdx] = 1;
						if ((*spectraSetCID)[i].getTolerance(prmIdx)
								< mergeRange.getTolerance()) {
							mergeRange.setMass((*spectraSetCID)[i][prmIdx][0]);
							mergeRange.setTolerance(
									(*spectraSetCID)[i].getTolerance(prmIdx));
						}
					}
				}
			}
			if (HCDscheckPRM != 0 && spectraSetHCD != 0) {
				for (list<int>::iterator idxIt = HCDscheckPRM->begin(); idxIt
						!= HCDscheckPRM->end(); idxIt++) {
					int i = *idxIt;
					int prmIdxTemp = (*spectraSetHCD)[i].findPeaks(mergeRange);
					if (prmIdxTemp >= 0 && (*usedPeaksHCD)[i][prmIdxTemp] != 2
							&& (*usedPeaksHCD)[i][prmIdxTemp] != 3) {
						prmIdx = prmIdxTemp;
						if ((*usedPeaksHCD)[i][prmIdx] == 0) {
							mergeRange.setIntensity(mergeRange.getIntensity()
									+ (*spectraSetHCD)[i][prmIdx][1]);
						}
						(*usedPeaksHCD)[i][prmIdx] = 1;
						if ((*spectraSetHCD)[i].getTolerance(prmIdx)
								< mergeRange.getTolerance()) {
							mergeRange.setMass((*spectraSetHCD)[i][prmIdx][0]);
							mergeRange.setTolerance(
									(*spectraSetHCD)[i].getTolerance(prmIdx));
						}
					}
				}
			}
			if (ETDscheckPRM != 0 && spectraSetETD != 0) {
				for (list<int>::iterator idxIt = ETDscheckPRM->begin(); idxIt
						!= ETDscheckPRM->end(); idxIt++) {
					int i = *idxIt;
					int prmIdxTemp = (*spectraSetETD)[i].findPeaks(mergeRange);
					if (prmIdxTemp >= 0 && (*usedPeaksETD)[i][prmIdxTemp] != 2
							&& (*usedPeaksETD)[i][prmIdxTemp] != 3) {
						prmIdx = prmIdxTemp;
						if ((*usedPeaksETD)[i][prmIdx] == 0) {
							mergeRange.setIntensity(mergeRange.getIntensity()
									+ (*spectraSetETD)[i][prmIdx][1]);
						}
						(*usedPeaksETD)[i][prmIdx] = 1;
						if ((*spectraSetETD)[i].getTolerance(prmIdx)
								< mergeRange.getTolerance()) {
							mergeRange.setMass((*spectraSetETD)[i][prmIdx][0]);
							mergeRange.setTolerance(
									(*spectraSetETD)[i].getTolerance(prmIdx));
						}
					}
				}
			}
			if (prmIdx < 0) {
				continue;
			}
		}

		(*CIDlabels)[j] = 2;
		(*ETDlabels)[closestETD] = 2;

		peaksToAdd.push_back(mergeRange);
	}
	addNewPeaks(maxPeakDensity, &peaksToAdd, mergedIdx, 1);
}

void PairedSpecSet::movePRMsSRMsPair(float maxPeakDensity, int mergedIdx,
		int CIDIdx, int ETDIdx, Spectrum::FragType fragCIDHCD) {
	SpecSet* CIDspecs;
	vector<vector<short> >* CIDused;

	if (fragCIDHCD == Spectrum::FragType_CID) {
		CIDspecs = spectraSetCID;
		CIDused = usedPeaksCID;
	} else {
		CIDspecs = spectraSetHCD;
		CIDused = usedPeaksHCD;
	}

	Spectrum* CIDspec = &(*CIDspecs)[CIDIdx];
	Spectrum* ETDspec = &(*spectraSetETD)[ETDIdx];
	Spectrum* mergedSpec = &(*mergedSet)[mergedIdx];

	if (mergedSpec->getPeakDensity() > maxPeakDensity) {
		return;
	}

	if (CIDspec->size() == 0 || ETDspec->size() == 0) {
		return;
	}

	vector<short>* CIDlabels = &(*CIDused)[CIDIdx];
	vector<short>* ETDlabels = &(*usedPeaksETD)[ETDIdx];
	vector<short>* mergedLocLabels = &(*mergedLabels)[mergedIdx];

	MZRange rangeCID;
	MZRange rangeETD;
	MZRange cPeakETD;
	MZRange mergeRange;
	list<MZRange> peaksToAdd;

	int numPeaks = CIDspec->size();

	float pmTol = mergedSpec->parentMassTol;

	float massBk = mergedSpec->parentMass - AAJumps::massMH;

	for (int j = 0; j < numPeaks; j++) {

		if ((*CIDlabels)[j] == 2 || (*CIDlabels)[j] == 3) {
			continue;
		}

		float massCID = (*CIDspec)[j][0];
		float intenCID = ((*CIDlabels)[j] == 0) ? (*CIDspec)[j][1] : 0;
		float massTolCID = CIDspec->getTolerance(j);
		rangeCID.set(massCID, intenCID, massTolCID);
		cPeakETD.set(massBk - massCID - AAJumps::massNH, intenCID, massTolCID
				+ pmTol);

		int closestETD = ETDspec->findClosest(cPeakETD.getMass());

		float massETD = (*ETDspec)[closestETD][0];
		float intenETD =
				((*ETDlabels)[closestETD] == 0) ? (*ETDspec)[closestETD][1] : 0;
		float massTolETD = ETDspec->getTolerance(closestETD);
		rangeETD.set(massETD, intenETD, massTolETD);

		if (rangeETD != cPeakETD || (*ETDlabels)[closestETD] == 1
				|| (*ETDlabels)[closestETD] == 3) {
			continue;
		}

		if ((*CIDlabels)[j] == 1 && (*ETDlabels)[closestETD] == 2) {
			continue;
		} else {
			(*CIDlabels)[j] = 1;
			(*ETDlabels)[closestETD] = 2;
		}

		mergeRange = rangeCID;
		mergeRange.setIntensity(intenCID + intenETD);

		peaksToAdd.push_back(mergeRange);
	}

	MZRange cPeakCID;
	MZRange closestETD;
	numPeaks = ETDspec->size();

	for (int j = 0; j < numPeaks; j++) {

		if ((*ETDlabels)[j] == 2 || (*ETDlabels)[j] == 3) {
			continue;
		}

		float massETD = (*ETDspec)[j][0];
		float intenETD = ((*ETDlabels)[j] == 0) ? (*ETDspec)[j][1] : 0;
		float massTolETD = ETDspec->getTolerance(j);
		rangeETD.set(massETD, intenETD, massTolETD);
		cPeakCID.set(massBk - massETD + AAJumps::massH2O, intenETD, massTolETD
				+ pmTol);

		int closestCID = CIDspec->findClosest(cPeakCID.getMass());

		float massCID = (*CIDspec)[closestCID][0];
		float intenCID =
				((*CIDlabels)[closestCID] == 0) ? (*CIDspec)[closestCID][1] : 0;
		float massTolCID = CIDspec->getTolerance(closestCID);
		rangeCID.set(massCID, intenCID, massTolCID);

		if (rangeCID != cPeakCID || (*CIDlabels)[closestCID] == 1
				|| (*ETDlabels)[closestCID] == 3) {
			continue;
		}

		if ((*ETDlabels)[j] == 1 && (*CIDlabels)[closestCID] == 2) {
			continue;
		} else {
			(*CIDlabels)[closestCID] = 2;
			(*ETDlabels)[j] = 1;
		}

		mergeRange = rangeETD;
		mergeRange.setIntensity(intenCID + intenETD);

		peaksToAdd.push_back(mergeRange);
	}
	addNewPeaks(maxPeakDensity, &peaksToAdd, mergedIdx, 1);
}

void PairedSpecSet::moveLeftoversPair(float maxPeakDensity, int mergedIdx,
		int childIdx, Spectrum::FragType fragCIDHCDETD) {
	SpecSet* childSpecs;
	vector<vector<short> >* childUsed;

	if (fragCIDHCDETD == Spectrum::FragType_CID) {
		childSpecs = spectraSetCID;
		childUsed = usedPeaksCID;
	} else if (fragCIDHCDETD == Spectrum::FragType_HCD) {
		childSpecs = spectraSetHCD;
		childUsed = usedPeaksHCD;
	} else {
		childSpecs = spectraSetETD;
		childUsed = usedPeaksETD;
	}

	Spectrum* childSpec = &(*childSpecs)[childIdx];
	Spectrum* mergedSpec = &(*mergedSet)[mergedIdx];

	if (mergedSpec->getPeakDensity() > maxPeakDensity) {
		return;
	}

	vector<short>* childLabels = &(*childUsed)[childIdx];
	vector<short>* mergedLocLabels = &(*mergedLabels)[mergedIdx];

	MZRange mergeRange;
	list<MZRange> peaksToAdd;

	int numPeaks = childSpec->size();

	for (int j = 0; j < numPeaks; j++) {
		if ((*childLabels)[j] != 0) {
			continue;
		}

		mergeRange.set((*childSpec)[j][0], (*childSpec)[j][1],
				childSpec->getTolerance(j));

		peaksToAdd.push_back(mergeRange);
	}
	addNewPeaks(maxPeakDensity, &peaksToAdd, mergedIdx, 0);
}

void PairedSpecSet::addNewPeaks(float maxPeakDensity, list<MZRange>* newPeaks,
		int mergedIdx, short label) {
	Spectrum* mergedSpec = &(*mergedSet)[mergedIdx];
	vector<short>* mergedLocLabels = &(*mergedLabels)[mergedIdx];

	if (mergedSpec->getPeakDensity() > maxPeakDensity && maxPeakDensity >= 0) {
		return;
	}

	newPeaks->sort(SortPeakIntensity());

	for (list<MZRange>::iterator peakIt = newPeaks->begin(); (mergedSpec->getPeakDensity()
			<= maxPeakDensity || maxPeakDensity < 0) && peakIt != newPeaks->end(); peakIt++) {
		int mergeIdx = mergedSpec->findPeaks(*peakIt);
		if (mergeIdx >= 0) {
			(*mergedSpec)[mergeIdx][1] += peakIt->getIntensity();
		} else {
			int insertIdx = mergedSpec->insertPeak(&(*peakIt));
			insertMergedLabel(mergedLocLabels, insertIdx, label);
		}
	}
}
}

