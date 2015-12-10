#include "spectrum.h"
#include "label.h"
#include "Logger.h"
#include "IsoEnvelope.h"
#include "utils.h"

#include <errno.h>
#include <string>
#include <map>

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <list>

namespace specnets
{
  // -------------------------------------------------------------------------
  Spectrum::Spectrum()
  {
    initialize();
  }
  // -------------------------------------------------------------------------
  Spectrum::~Spectrum()
  {
    //DEBUG_MSG("destructor");
    //	std::cerr << "destructor" << endl;
    if (info != (char *)0)
    {
    }
    else
    {
      free(info);
      info = (char *)0;
    }

    list<psmPtr>::iterator it;
    //set Spectrum pointer to null for associated PSMs
    for (it = psmList.begin(); it != psmList.end(); it++)
    {
      (*it)->m_spectrum = (Spectrum *)NULL;
    }
  }
  // -------------------------------------------------------------------------
  Spectrum::Spectrum(const Spectrum &other)
  {
	initialize();
    internalCopy(other);
  }

  void Spectrum::initialize() {
	// DEBUG_MSG("constructor");
	//	std::cerr << "constructor" << endl;
	parentMass = 0;
	parentMassTol = 0;
	oldParentMassTol = 0;
	parentCharge = 0;
	parentMZ = 0;
	scan = 0;
	msLevel = 0;
	msFragType = FragType_CID;
	peakList.resize(0);
	psmList.resize(0);
	peakTols.resize(0);
	oldPeakTols.resize(0);
	resolution = 1;
	idDist = 1;
	info = (char *) 0;
	instrument_name = "";
	ITOL = -1.0;
	ITOLU = "";
	TOL = -1.0;
	TOLU = "";
	spectrum_quality = -1;
  }

  Spectrum &Spectrum::operator=(const Spectrum &other)
  {
    //	std::cerr << "assignment" << endl;

    if (this == &other)
      return *this;
    internalCopy(other);
    return (*this);
  }
  // -------------------------------------------------------------------------
  Spectrum &Spectrum::copyNP(const Spectrum &other)
  {
    //	std::cerr << "copynp" << endl;
    parentMass = other.parentMass;
    parentMassTol = other.parentMassTol;
    parentCharge = other.parentCharge;
    parentMZ = other.parentMZ;
    scan = other.scan;
    msLevel = other.msLevel;
    msFragType = other.msFragType;
    resolution = other.resolution;
    idDist = other.idDist;
    instrument_name = other.instrument_name;
    ITOL = other.ITOL;
    ITOLU = other.ITOLU;
    TOL = other.TOL;
    TOLU = other.TOLU;
    spectrum_quality = other.spectrum_quality;
    
    return (*this);
  }
  // -------------------------------------------------------------------------
  void Spectrum::internalCopy(const Spectrum &other)
  {
    copyNP(other);

    peakList.resize(other.peakList.size());
    for (unsigned int i = 0; i < other.peakList.size(); i++)
      peakList[i] = other.peakList[i];
    peakTols.resize(other.peakTols.size());
    for (unsigned int i = 0; i < other.peakTols.size(); i++) {
      peakTols[i] = other.peakTols[i];
    }
    oldPeakTols = other.oldPeakTols;

    //copy psms
    psmList.clear();
    for (list<psmPtr>::const_iterator it = other.psmList.begin(); it
        != other.psmList.end(); it++)
    {
      psmPtr newPsm(new PeptideSpectrumMatch);
      *newPsm = **it;
      psmList.push_back(newPsm);
    }
    return;
  }
    
  // -------------------------------------------------------------------------
  inline void Spectrum::weightedAverage(int minIdx,
                                        int maxIdx,
                                        float &avgMass,
                                        float &totScore)
  {
    int i = 0;
    totScore = 0;
    avgMass = 0;
    for (i = minIdx; i <= maxIdx; i++)
      totScore += peakList[i][1];
    if (totScore <= 0)
      for (i = minIdx; i <= maxIdx; i++)
        avgMass += peakList[i][0] / (maxIdx - minIdx + 1); // regular average
    else
      for (i = minIdx; i <= maxIdx; i++)
        avgMass += peakList[i][0] * peakList[i][1] / totScore; // weighted average
  }

  void Spectrum::setResolution(float newResolution, bool enforceRounding)
  {
    resolution = newResolution;
    idDist = 1 / resolution;
    parentMass = enforceRounding ? round(parentMass / resolution) : parentMass
        / resolution;
    if (!enforceRounding)
      for (unsigned int i = 0; i < peakList.size(); i++)
        peakList[i][0] = peakList[i][0] / resolution;
    else
    {
      unsigned int j = 0;
      for (unsigned int i = 0; i < peakList.size(); i++)
      {
        peakList[i][0] = round(peakList[i][0] / resolution);
        if (j < i)
        {
          if (abs(peakList[j][0] - peakList[i][0]) < 0.0001)
            peakList[j][1] += peakList[i][1];
          else
          {
            j++;
            peakList[j] = peakList[i];
            peakTols[j] = peakTols[i];
          }
        }
      }
      resize(j + 1);
    }
  }

  /**
 * Sets the tolerance of each peak in Da
 */
void Spectrum::setPeakTolerance(float peakTol) {
	for (int i = 0; i < size(); i++) {
		peakTols[i] = peakTol;
	}
}

/**
 * Sets the ppm tolerance of each peak
 */
void Spectrum::setPeakTolerancePPM(float peakTolPPM) {
	for (int i = 0; i < size(); i++) {
		peakTols[i] = peakList[i][0] * peakTolPPM * PPM_FACTOR;
	}
}

TwoValues<float>* Spectrum::front() {
	return &peakList.front();
}

TwoValues<float>* Spectrum::back() {
	return &peakList.back();
}

/**
 * Sets the peak tolerance of the parent mass
 */
void Spectrum::setParentMassTol(float pmTol) {
	parentMassTol = pmTol;
}

/**
 * Sets the ppm tolerance of the parent mass
 */
void Spectrum::setParentMassTolPPM(float pmTolPPM) {
	parentMassTol = parentMass * pmTolPPM * PPM_FACTOR;
}

/**
 * Copies all per-peak tolerances into an internal vector
 */
void Spectrum::rememberTolerances() {
	oldPeakTols = peakTols;
	oldParentMassTol = parentMassTol;
}

/**
 * Restores peak tolerances to as they were before last call to rememberTolerances()
 */
void Spectrum::revertTolerances() {
	if (oldPeakTols.size() == peakTols.size()) {
		peakTols = oldPeakTols;
		parentMassTol = oldParentMassTol;
	} else {
		WARN_MSG("Could not revert peak tolerances because peak list size has changed!");
		return;
	}
}

float Spectrum::getTolerance(int index) {
	return peakTols[index];
}

void Spectrum::setTolerance(int index, float tolerance) {
	peakTols[index] = tolerance;
}

void Spectrum::sortPeaks() {
	map<TwoValues<float>*, float> peakTolMap;

	for (int i = 0; i < peakList.size(); i++) {
		peakTolMap[&peakList[i]] = peakTols[i];
	}

	sort(peakList.begin(), peakList.end());

	for (int i = 0; i < peakList.size(); i++) {
		peakTols[i] = peakTolMap[&peakList[i]];
	}
}

int Spectrum::insertPeak(float mass, float intensity, float tolerance,
		int atIndex) {
	vector<TwoValues<float> >::iterator peakIt = peakList.begin();
	vector<float>::iterator tolIt = peakTols.begin();
	int insertIdx = 0;

	if (atIndex >= 0) {
		if (atIndex >= peakList.size() || atIndex >= peakTols.size()) {
			atIndex = min(peakList.size(), peakTols.size()) - 1;
		}
		peakIt += atIndex;
		tolIt += atIndex;
		insertIdx = atIndex;
	} else {
		if (size() > 0) {
			insertIdx = findClosest(mass);
			if (peakList[insertIdx][0] < mass) {
				insertIdx++;
			}
		} else {
			insertIdx = 0;
		}
		peakIt += insertIdx;
		tolIt += insertIdx;
	}

	peakList.insert(peakIt, TwoValues<float> (mass, intensity));
	peakTols.insert(tolIt, tolerance);

	return insertIdx;
}

int Spectrum::insertPeak(MZRange* peak, int atIndex) {
	return insertPeak(peak->getMass(), peak->getIntensity(),
			peak->getTolerance(), atIndex);
}

int Spectrum::insertPeaks(list<MZRange>& newPeaks) {

	int idxUse = size() - 1;
	resize(size() + newPeaks.size());
	for (list<MZRange>::iterator peakIt = newPeaks.begin(); peakIt != newPeaks.end(); peakIt ++) {
		peakList[idxUse][0] = peakIt->getMass();
		peakList[idxUse][1] = peakIt->getIntensity();
		peakTols[idxUse] = peakIt->getTolerance();
		++idxUse;
	}
	sortPeaks();
}

int Spectrum::findPeaks(float mass, float tolerance, list<int>* matches,
		int* startIdx) {
	if (size() == 0) {
		return -1;
	}
	int idx;

	if (startIdx == (int*) 0) {
		idx = findClosest(mass);
	} else {
		idx = *startIdx;
		if (idx < 0) {
			idx = 0;
		}
		while (idx < size() && peakList[idx][0] < mass) {
			++idx;
		}
		if (idx > 0 && peakList[idx][0] > mass && peakList[idx - 1][0] < mass
				&& peakList[idx][0] - mass > mass - peakList[idx - 1][0]) {
			--idx;
		}
	}
	MZRange peakRange(peakList[idx][0], 0, peakTols[idx]);
	if (peakRange == mass || peakRange == mass - tolerance || peakRange == mass
			+ tolerance) {
		if (matches != 0) {
			matches->clear();
			matches->push_back(idx);
			float lowerBound = mass - tolerance;
			for (int i = idx - 1; i >= 0 && peakList[i][0] + peakTols[i]
					>= lowerBound; i--) {
				matches->push_back(i);
			}
			float upperBound = mass + tolerance;
			for (int i = idx + 1; i < size() && peakList[i][0] - peakTols[i]
					<= upperBound; i++) {
				matches->push_back(i);
			}
			matches->sort();
		}
		return idx;
	} else {
		if (matches != 0) {
			matches->clear();
		}
		return -1;
	}
}

int Spectrum::findPeaks(const MZRange& range, list<int>* matches, int* startIdx) {
	return findPeaks(range.getMass(), range.getTolerance(), matches, startIdx);
}

bool Spectrum::removePeak(int peakIndex) {
	if (peakIndex < 0 || peakIndex >= peakList.size() || peakIndex
			>= peakTols.size()) {
		WARN_MSG("Invalid peak index " << peakIndex);
		return false;
	}
	vector<TwoValues<float> >::iterator peakIt = peakList.begin();
	vector<float>::iterator tolIt = peakTols.begin();

	peakList.erase(peakIt + peakIndex);
	peakTols.erase(tolIt + peakIndex);

	return true;
}

bool Spectrum::removePeaks(list<int>& peakIndices) {

	if (size() < peakIndices.size()) {
		return false;
	}
	vector<bool> removed(size(), false);

	for (list<int>::iterator peakIt = peakIndices.begin(); peakIt
			!= peakIndices.end(); peakIt++) {
		if (*peakIt < 0 || *peakIt >= size()) {
			return false;
		}
		removed[*peakIt] = true;
	}

	vector<TwoValues<float> > newPeakList(size());
	vector<float> newTols(size());
	int idxUse = 0;
	for (int i = 0; i < removed.size(); i++) {
		if (! removed[i]) {
			newPeakList[idxUse] = peakList[i];
			newTols[idxUse] = peakTols[i];
			++idxUse;
		}
	}
	newPeakList.resize(idxUse);
	newTols.resize(idxUse);
	resize(idxUse);

	peakList.assign(newPeakList.begin(), newPeakList.end());
	peakTols.assign(newTols.begin(), newTols.end());

	return true;
}

float Spectrum::getPeakDensity() {
	return ((float) size()) / parentMass;
}

// -------------------------------------------------------------------------
void Spectrum::addZPMpeaks(float tolerance, float ionOffset, bool includeY0k,
		bool ctermH2O, SpectrumPeakLabels *labels) {
	float massB0 = ionOffset * idDist, massBk = parentMass
			- (ctermH2O ? AAJumps::massMH : AAJumps::massHion - ionOffset)
					* idDist, massY0 = (ctermH2O ? AAJumps::massH2O : 0
			+ ionOffset) * idDist, massYk = parentMass - (AAJumps::massHion
			- ionOffset) * idDist;

	float endptIntensity = 0.1;

	bool modTols = false;

	if (tolerance >= 0) {
		modTols = true;
		rememberTolerances();
		setPeakTolerance(tolerance);
	}

	if (findPeaks(massB0) < 0) {
		int idx = insertPeak(massB0, endptIntensity, 0);
		if (modTols) {
			vector<float>::iterator peakIt = oldPeakTols.begin();
			peakIt += idx;
			oldPeakTols.insert(peakIt, 0);
		}
		if (labels != 0) {
			vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
			labelIt += idx;
			PeakLabel newLabel(0, 0, 0, 0, 1);
			labels->peakLabels.insert(labelIt, newLabel);
		}
	}

	if (findPeaks(massBk) < 0) {
		int idx = insertPeak(massBk, endptIntensity, 0);
		if (modTols) {
			vector<float>::iterator peakIt = oldPeakTols.begin();
			peakIt += idx;
			oldPeakTols.insert(peakIt, 0);
		}
		if (labels != 0) {
			vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
			labelIt += idx;
			PeakLabel newLabel(0, 0, 0, 0, 1);
			labels->peakLabels.insert(labelIt, newLabel);
		}
	}

	if (findPeaks(massY0) < 0) {
		int idx = insertPeak(massY0, endptIntensity, 0);
		if (modTols) {
			vector<float>::iterator peakIt = oldPeakTols.begin();
			peakIt += idx;
			oldPeakTols.insert(peakIt, 0);
		}
		if (labels != 0) {
			vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
			labelIt += idx;
			PeakLabel newLabel(0, 0, 0, 0, 1);
			labels->peakLabels.insert(labelIt, newLabel);
		}
	}

	if (includeY0k && findPeaks(massYk) < 0) {
		int idx = insertPeak(massYk, endptIntensity, 0);
		if (modTols) {
			vector<float>::iterator peakIt = oldPeakTols.begin();
			peakIt += idx;
			oldPeakTols.insert(peakIt, 0);
		}
		if (labels != 0) {
			vector<PeakLabel>::iterator labelIt = labels->peakLabels.begin();
			labelIt += idx;
			PeakLabel newLabel(0, 0, 0, 0, 1);
			labels->peakLabels.insert(labelIt, newLabel);
		}
	}

	if (modTols) {
		revertTolerances();
	}
}
// -------------------------------------------------------------------------
void Spectrum::maximizeZPMpeaks(float tolerance, float ionOffset,
		bool includeY0k) {

	bool modTols = false;

	if (tolerance >= 0) {
		modTols = true;
		rememberTolerances();
		setPeakTolerance(tolerance);
	}

	float massB0 = ionOffset * idDist, massBk = parentMass - (AAJumps::massMH
			- ionOffset) * idDist, massY0 = (AAJumps::massH2O + ionOffset)
			* idDist, massYk = parentMass - (AAJumps::massHion - ionOffset)
			* idDist;

	unsigned int pivot;
	float maxScore = 0;
	for (pivot = 0; pivot < peakList.size(); pivot++)
		if (peakList[pivot][1] > maxScore)
			maxScore = peakList[pivot][1];
	for (pivot = 0; pivot < peakList.size() and peakList[pivot][0] <= massB0
			+ getTolerance(pivot); pivot++)
		peakList[pivot][1] = maxScore;
	if (includeY0k) {
		for (; pivot < peakList.size() and peakList[pivot][0] < massY0
				- getTolerance(pivot); pivot++)
			;
		for (; pivot < peakList.size() and peakList[pivot][0] <= massY0
				+ getTolerance(pivot); pivot++)
			peakList[pivot][1] = maxScore;
	}
	for (; pivot < peakList.size() and peakList[pivot][0] < massBk; pivot++)
		;
	for (; pivot < peakList.size() and peakList[pivot][0] <= massBk; pivot++)
		peakList[pivot][1] = maxScore;
	if (includeY0k) {
		for (; pivot < peakList.size() and peakList[pivot][0] < massYk
				- getTolerance(pivot); pivot++)
			;
		for (; pivot < peakList.size() and peakList[pivot][0] <= massYk
				+ getTolerance(pivot); pivot++)
			peakList[pivot][1] = maxScore;
	}

	if (modTols) {
		revertTolerances();
	}
}
  // -------------------------------------------------------------------------
  void Spectrum::normalize(float newTotalInt, bool removeNegatives)
  { // Normalizes total intensity to 100
    float totalIntensity = 0;
    for (unsigned int i = 0; i < peakList.size(); i++)
      totalIntensity += peakList[i][1];
    for (unsigned int i = 0; i < peakList.size(); i++)
      peakList[i][1] = newTotalInt * peakList[i][1] / totalIntensity;
  }
  // -------------------------------------------------------------------------
  void Spectrum::normalize2(float newNorm2)
  { // Normalizes to Euclidian norm newNorm2
    float factor = 0;
    for (unsigned int i = 0; i < peakList.size(); i++)
      factor += peakList[i][1] * peakList[i][1];
    factor = sqrt(factor) / newNorm2;
    for (unsigned int i = 0; i < peakList.size(); i++)
      peakList[i][1] = peakList[i][1] / factor;
  }
  // -------------------------------------------------------------------------
  void Spectrum::guessPrecursorZPM(Spectrum &parentMS1,
                                   float peakTol,
                                   short maxZ,
                                   IsoEnvelope &isoEnvs,
                                   bool strictMode)
  {

	  WARN_MSG("guessPrecursorZPM has not been updated for per-peak tolerances");
	  cerr << "guessPrecursorZPM has not been updated for per-peak tolerances\n";

    // Look for the best isotopic envelope match over all possible charges and monoisotopic masses
    float curMonoMass, chargeIncrement, bestPM = 0, curMatch, bestMatch =
        1000000; // Match score of zero is optimal
    short bestZ = 0;
    vector<float> massEnvelope;
    //cerr<<"Processing scan "<<scan<<endl;
    for (short charge = 1; charge <= maxZ; charge++)
    {
      chargeIncrement = 1.0 / ((float)charge);
      for (short massOffset = -2; massOffset <= 0; massOffset++)
      {
        curMonoMass = parentMass + ((float)massOffset) * chargeIncrement;
        isoEnvs.ExtractEnvelope(curMonoMass,
                                charge,
                                parentMS1,
                                peakTol,
                                massEnvelope,
                                strictMode);
        //cerr<<"massEnvelope ["<<charge<<","<<massOffset<<"]: "; for(unsigned int i=0; i<massEnvelope.size(); i++) cerr<<massEnvelope[i]<<" "; cerr<<"\n"; cerr.flush();
        isoEnvs.normalize(massEnvelope);
        //cerr<<"normalized massEnvelope: "; for(unsigned int i=0; i<massEnvelope.size(); i++) cerr<<massEnvelope[i]<<" "; cerr<<"\n"; cerr.flush();
        curMatch = isoEnvs.ScoreEnvelope(curMonoMass, massEnvelope, strictMode);
        //cerr<<"score = "<<curMatch<<endl;
        if (curMatch < bestMatch)
        {
          bestMatch = curMatch;
          bestPM = curMonoMass * charge - charge + 1;
          bestZ = charge;
        }
      }
    }
    //cerr<<" --> Chose parentCharge = "<<bestZ<<", parentMass = "<<bestPM<<" (was "<<parentMass<<")\n";
    parentCharge = bestZ;
    parentMass = bestPM;
  }
  // -------------------------------------------------------------------------
  short Spectrum::findMatches(float baseMass,
                              float peakTol,
                              vector<int> &matchesIdx,
                              int startIdx)
  {
	vector<float> tempTols;
	bool modTols = false;

	if (peakTol >= 0) {
		modTols = true;
		rememberTolerances();
		setPeakTolerance(peakTol);
	}
	int* startPtr = &startIdx;
	if (startIdx < 0) {
		startPtr = (int*) 0;
	}

	list<int> matchesList;
	int idx = findPeaks(baseMass, 0, &matchesList, startPtr);
	matchesIdx.resize(matchesList.size());
	int idxUse = 0;
	for (list<int>::iterator mIt = matchesList.begin(); mIt != matchesList.end(); mIt++) {
		if (*mIt >= startIdx) {
			matchesIdx[idxUse] = *mIt;
			idxUse ++;
		}
	}
	matchesIdx.resize(idxUse);

	if (modTols) {
		revertTolerances();
	}

	return (short) idxUse;
	/*
    unsigned int baseIdx = (startIdx >= 0 and startIdx < peakList.size())
        ? (unsigned int)startIdx : 0;
    unsigned int numPeaks = peakList.size();
    unsigned int matchCount = 0;
    matchesIdx.resize(numPeaks);
    for (; baseIdx > 0 and peakList[baseIdx][0] >= baseMass - peakTol; baseIdx--)
      ;
    for (; baseIdx < numPeaks and peakList[baseIdx][0] < baseMass - peakTol; baseIdx++)
      ; // Get to first peak above the lower end of the tolerance window
    for (; baseIdx < numPeaks and peakList[baseIdx][0] <= baseMass + peakTol; baseIdx++)
      matchesIdx[matchCount++] = baseIdx;
    matchesIdx.resize(matchCount);
    return matchCount;
    */
  }
  // -------------------------------------------------------------------------
  int Spectrum::findClosest(float mass)
  {
    if (peakList.size() == 0)
      return -1;

    int lowerLim = 0, upperLim = (int)peakList.size() - 1, middle;
    while (upperLim > lowerLim + 1)
    { // Set lowerLim to the index of the largest monoisotopic mass < mass
      middle = (lowerLim + upperLim) / 2;
      if (peakList[middle][0] > mass)
        upperLim = middle;
      else
        lowerLim = middle;
    }
    if (fabs(peakList[lowerLim][0] - mass) < fabs(peakList[upperLim][0] - mass))
      return lowerLim;
    else
      return upperLim;
  }
  // -------------------------------------------------------------------------
  void Spectrum::mergePeakList(vector<TwoValues<float> > &newPeaks,
                               Spectrum *putHere)
  {
	  WARN_MSG("mergePeakList has not been updated for per-peak tolerances");
	  cerr << "mergePeakList has not been updated for per-peak tolerances\n";


    int idxUpdated, idxOld, idxNew;
    vector<TwoValues<float> > &updatedPeakList = (putHere == 0) ? peakList
        : putHere->peakList;
    vector<TwoValues<float> > *oldPeakList;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> > ;
      oldPeakList->resize(peakList.size());
      for (int i = 0; i < peakList.size(); i++)
        (*oldPeakList)[i] = peakList[i];
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      putHere->parentMass = parentMass;
      putHere->parentCharge = parentCharge;
    }

    // Sorted merge
    updatedPeakList.resize(oldPeakList->size() + newPeaks.size());
    idxUpdated = 0;
    idxOld = 0;
    idxNew = 0;
    while (idxOld < oldPeakList->size() or idxNew < newPeaks.size())
    {
      if (idxOld == oldPeakList->size())
      {
        updatedPeakList[idxUpdated++] = newPeaks[idxNew++];
        continue;
      }
      if (idxNew == newPeaks.size())
      {
        updatedPeakList[idxUpdated++] = (*oldPeakList)[idxOld++];
        continue;
      }

      if (fabs((*oldPeakList)[idxOld][0] - newPeaks[idxNew][0]) < 0.0001)
      {
        updatedPeakList[idxUpdated][0] = (*oldPeakList)[idxOld][0];
        updatedPeakList[idxUpdated++][1] = (*oldPeakList)[idxOld++][1]
            + newPeaks[idxNew++][1];
      }
      else if ((*oldPeakList)[idxOld][0] < newPeaks[idxNew][0])
        updatedPeakList[idxUpdated++] = (*oldPeakList)[idxOld++];
      else
        updatedPeakList[idxUpdated++] = newPeaks[idxNew++];
    }
    updatedPeakList.resize(idxUpdated); // Just in case some peaks get merged along the way

    if (putHere == 0) {
      delete oldPeakList;
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::mergeClosestPeaks(Spectrum& newSpec,
                                   short mergeType)
  {
    if (size() == 0) {
      peakList = newSpec.peakList;
      peakTols = newSpec.peakTols;
    }

    list<int> peaksAppend;
    map<int, list<MZRange> > peaksMerge;
    float mass1, intensity1, tolerance1, mass2, intensity2, tolerance2;
    int closestIdx;
    MZRange nextRange1, nextRange2;
    list<MZRange> mergePeaks;

    for (int i = 0; i < newSpec.size(); i++) {
      mass1 = newSpec[i][0];
      intensity1 = newSpec[i][1];
      tolerance1 = newSpec.peakTols[i];
      nextRange1.set(mass1, intensity1, tolerance1);

      closestIdx = findClosest(mass1);
      mass2 = peakList[closestIdx][0];
      intensity2 = peakList[closestIdx][1];
      tolerance2 = peakTols[closestIdx];
      nextRange2.set(mass2, intensity2, tolerance2);

      if (MZRange::EqualWithinRange(mass1, mass2, tolerance1 + tolerance2)) {
        if (peaksMerge.count(closestIdx) == 0) {
          mergePeaks.clear();
          mergePeaks.push_back(nextRange2);
          mergePeaks.push_back(nextRange1);
          peaksMerge[closestIdx] = mergePeaks;
        }
        else {
          peaksMerge[closestIdx].push_back(nextRange1);
        }
      }
      else {
        peaksAppend.push_back(i);
      }
    }

    for (map<int, list<MZRange> >::iterator mergeIt = peaksMerge.begin(); mergeIt
        != peaksMerge.end(); mergeIt++) {
      mass1 = peakList[mergeIt->first][0];
      intensity1 = peakList[mergeIt->first][1];
      tolerance1 = peakTols[mergeIt->first];
      nextRange1.set(mass1, intensity1, tolerance1);

      switch (mergeType) {
        case 0:
          nextRange1.MergeMZRanges(&mergeIt->second, 1);
        case 1:
          nextRange1.MergeMZRanges(&mergeIt->second, 2);
        case 2:
          nextRange1.MergeMZRanges(&mergeIt->second, 3);
        case 3:
          nextRange1.MergeMZRanges(&mergeIt->second, 0);
      }

      peakList[mergeIt->first][0] = nextRange1.getMass();
      peakList[mergeIt->first][1] = nextRange1.getIntensity();
      peakTols[mergeIt->first] = nextRange1.getTolerance();
    }

    int idx = peakList.size();
    peakList.resize(peakList.size() + peaksAppend.size());
    peakTols.resize(peakList.size());
    for (list<int>::iterator peakIt = peaksAppend.begin(); peakIt
        != peaksAppend.end(); peakIt++) {
      peakList[idx] = newSpec[*peakIt];
      peakTols[idx] = newSpec.peakTols[*peakIt];
      idx ++;
    }

    sortPeaks();
  }

  // -------------------------------------------------------------------------
  void Spectrum::massesToIndices(Spectrum &masses,
                                 vector<int> &indices,
                                 float peakTol)
  {
	vector<float> tempTols;
	bool modTols = false;

	if (peakTol >= 0) {
		modTols = true;
		rememberTolerances();
		setPeakTolerance(peakTol);
	}

    int idxClosest, idxPeaks = 0, idxMasses;
    float distClosest;// = peakTol + 1; // Supremum for distance to closest peak
    indices.resize(masses.size());
    for (idxMasses = 0; idxMasses < (int)masses.size(); idxMasses++)
    {
      while (idxPeaks > 0 and peakList[idxPeaks][0] > masses[idxMasses][0]
          - getTolerance(idxPeaks))
        idxPeaks--;
      while (idxPeaks < (int)peakList.size() and peakList[idxPeaks][0]
          < masses[idxMasses][0] - getTolerance(idxPeaks))
        idxPeaks++;
      idxClosest = -1;
      distClosest = getTolerance(idxPeaks) + 1.0;
      while (idxPeaks < (int)peakList.size() and abs(peakList[idxPeaks][0]
          - masses[idxMasses][0]) <= getTolerance(idxPeaks) + 0.0001)
      {
        if (abs(peakList[idxPeaks][0] - masses[idxMasses][0]) < distClosest)
        {
          idxClosest = idxPeaks;
          distClosest = abs(peakList[idxPeaks][0] - masses[idxMasses][0]);
        }
        idxPeaks++;
      }
      if (idxClosest >= 0)
        indices[idxMasses] = idxClosest;
      else
        indices[idxMasses] = -1;
      idxPeaks = min(idxPeaks, (int)peakList.size() - 1);
    }

    if (modTols) {
		revertTolerances();
	}
  }
  // -------------------------------------------------------------------------
  void Spectrum::selectIndices(vector<int> &idx)
  {
	  WARN_MSG("selectIndices has not been updated for per-peak tolerances");
	  cerr << "selectIndices has not been updated for per-peak tolerances\n";

    unsigned int i;
    for (i = 0; i < idx.size(); i++)
      if (idx[i] < peakList.size())
        peakList[i] = peakList[idx[i]];
      else
        break;
    resize(i);
  }
  // -------------------------------------------------------------------------
unsigned int Spectrum::filterLowIntensity(float minIntensity) {
	list<int> peaksToRemove;

	for (int i = 0; i < size(); i++) {
		if (peakList[i][1] < minIntensity) {
			peaksToRemove.push_back(i);
		}
	}

	removePeaks(peaksToRemove);

	return size();
}

  // -------------------------------------------------------------------------
/**
 * Removes peaks that have an intensity ranked worse than maxRank
 *   compared to all neighboring peaks +/- windowRadius
 * @param maxRank maximum allowable rank of each peak
 * @param windowRadius radius of peak comparison in the spectrum
 * @return
 */
void Spectrum::rankFilterPeaks(int maxRank, float windowRadius) {
	list<int> peaksToRemove;

	list<TwoValues<float> > rankedPeaks;
	int j, k, curRank, peakIdx;
	float mass, intensity, minRange, maxRange, prevInten;
	bool foundPeak;
	int numPeaks = size();
	for (int i = 0; i < numPeaks; i++) {
		rankedPeaks.clear();
		mass = peakList[i][0];
		intensity = peakList[i][1];
		minRange = mass - windowRadius;
		maxRange = mass + windowRadius;
		rankedPeaks.push_back(TwoValues<float> (intensity, (float) i));
		j = i - 1;
		while (j >= 0 && peakList[j][0] > minRange) {
			rankedPeaks.push_back(TwoValues<float> (peakList[j][1], (float) j));
			--j;
		}
		j = i + 1;
		while (j < numPeaks && peakList[j][0] < maxRange) {
			rankedPeaks.push_back(TwoValues<float> (peakList[j][1], (float) j));
			++j;
		}

		rankedPeaks.sort();
		rankedPeaks.reverse();
		curRank = 0;
		foundPeak = false;
		prevInten = -1.0;
		//cout << "ranked: ";
		for (list<TwoValues<float> >::iterator rankIt = rankedPeaks.begin(); rankIt
				!= rankedPeaks.end(); rankIt++) {
			//cout << "(" << (*rankIt)[0] << "," << (*rankIt)[1] << "), ";
			peakIdx = floatToInt((*rankIt)[1]);
			if ((*rankIt)[0] != prevInten) {
				curRank++;
				prevInten = (*rankIt)[0];
			}
			if (curRank > maxRank) {
				break;
			}
			if (peakIdx == i) {
				foundPeak = true;
				break;
			}
		}
		//cout << "\n";
		if (!foundPeak) {
			peaksToRemove.push_back(i);
		}
	}

	removePeaks(peaksToRemove);
}

  void Spectrum::output(ostream &output)
  {
    output << parentMass << " " << parentCharge << endl;
    for (unsigned int i = 0; i < peakList.size(); i++)
      output << peakList[i][0] << " " << peakList[i][1] << " " << peakTols[i] << endl;
  }
  // -------------------------------------------------------------------------
  void Spectrum::output_ms2(ostream &output)
  {
    //output << ":0.0.0\n";
    output << parentMass << " " << parentCharge << endl;
    for (unsigned int i = 0; i < peakList.size(); i++)
      output << peakList[i][0] << " " << peakList[i][1] << endl;
  }
  // -------------------------------------------------------------------------
  void Spectrum::reverse(float pmOffset, Spectrum *putHere)
  {
    vector<TwoValues<float> > &updatedPeakList = (putHere == 0) ? peakList
        : putHere->peakList;
    vector<float> &updatedPeakTols = (putHere == 0) ? peakTols
            : putHere->peakTols;

    vector<TwoValues<float> > *oldPeakList;
    vector<float> *oldPeakTols;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> > ;
      oldPeakTols = new vector<float>;
      oldPeakList->resize(peakList.size());
      oldPeakTols->resize(peakList.size());
      for (unsigned int i = 0; i < peakList.size(); i++) {
        (*oldPeakList)[i] = peakList[i];
        (*oldPeakTols)[i] = peakTols[i];
      }
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      oldPeakTols = &peakTols;
      (*putHere) = *this;
    }

    float totMass = parentMass - AAJumps::massHion + pmOffset;
    unsigned int numPeaks = oldPeakList->size();
    for (unsigned int i = 0; i < numPeaks; i++) {
      updatedPeakList[numPeaks - i - 1].set(totMass - (*oldPeakList)[i][0],
                                            (*oldPeakList)[i][1]);
      updatedPeakTols[numPeaks - i - 1] = (*oldPeakTols)[i];
    }

    if (putHere == 0) {
      delete oldPeakList;
      delete oldPeakTols;
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::rotate(float offset, float cterm)
  {
	  WARN_MSG("rotate has not been updated for per-peak tolerances");
	  cerr << "rotate has not been updated for per-peak tolerances\n";

    vector<TwoValues<float> > tmpPeaks;
    int peakIdx, v;

    if (offset < 0)
    {
      for (peakIdx = 0; peakIdx < (int)peakList.size(); peakIdx++)
        if (peakList[peakIdx][0] > -offset)
          break;
      int firstAbove = peakIdx;

      tmpPeaks.resize(firstAbove);
      for (peakIdx = 0; peakIdx < firstAbove; peakIdx++) // Compute peaks with masses lower than offset
        tmpPeaks[peakIdx].set(parentMass - cterm + peakList[peakIdx][0]
            + offset, peakList[peakIdx][1]);
      for (peakIdx = 0; peakIdx < (int)peakList.size() - firstAbove; peakIdx++) // Set peaks with masses higher than offset
      {
        v = peakIdx + firstAbove;
        peakList[peakIdx].set(peakList[v][0] + offset, peakList[v][1]);
      }
      v = peakList.size() - firstAbove;
      for (peakIdx = 0; peakIdx < firstAbove; peakIdx++) // Set peaks with masses lower than offset
        peakList[v + peakIdx] = tmpPeaks[peakIdx];
    }
    else
    {
      for (peakIdx = 0; peakIdx < (int)peakList.size(); peakIdx++)
        if (peakList[peakIdx][0] > parentMass - cterm - offset)
          break;
      int firstWrap = peakIdx, numWrap = peakList.size() - firstWrap;

      tmpPeaks.resize(numWrap);
      for (peakIdx = firstWrap; peakIdx < (int)peakList.size(); peakIdx++) // Compute peaks with masses higher than offset
        tmpPeaks[peakIdx - firstWrap].set(peakList[peakIdx][0] + offset
            - (parentMass - cterm), peakList[peakIdx][1]);
      for (peakIdx = (int)peakList.size() - 1; peakIdx >= numWrap; peakIdx--) // Set peaks with masses higher than offset
      {
        v = peakIdx - numWrap;
        peakList[peakIdx].set(peakList[v][0] + offset, peakList[v][1]);
      }
      for (peakIdx = 0; peakIdx < numWrap; peakIdx++) // Set peaks with masses lower than offset
        peakList[peakIdx] = tmpPeaks[peakIdx];
    }
  }
  // -------------------------------------------------------------------------
  void Spectrum::selectTopK(unsigned int topK, Spectrum *putHere)
  {
	  WARN_MSG("selectTopK has not been updated for per-peak tolerances and should be replaced by rankFilterPeaks");
	  cerr << "selectTopK has not been updated for per-peak tolerances and should be replaced by rankFilterPeaks\n";


    vector<TwoValues<float> > &updatedPeakList = (putHere == 0) ? peakList
        : putHere->peakList;
    vector<TwoValues<float> > *oldPeakList;
    unsigned int peakIdx;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> > ;
      oldPeakList->resize(peakList.size());
      for (unsigned int i = 0; i < peakList.size(); i++)
        (*oldPeakList)[i] = peakList[i];
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      (*putHere) = *this;
    }
    vector<TwoValues<float> > sortedInts(oldPeakList->size());

    if (topK < oldPeakList->size())
    {
      for (peakIdx = 0; peakIdx < oldPeakList->size(); peakIdx++)
        sortedInts[peakIdx].set((*oldPeakList)[peakIdx][1], peakIdx);
      sort(sortedInts.begin(), sortedInts.end());

      updatedPeakList.resize(topK);
      for (peakIdx = 1; peakIdx <= topK; peakIdx++)
        updatedPeakList[peakIdx - 1]
            = (*oldPeakList)[(unsigned int)sortedInts[sortedInts.size()
                - peakIdx][1]];

      sort(updatedPeakList.begin(), updatedPeakList.end());
    }

    if (putHere == 0)
      delete oldPeakList;
  }
  // -------------------------------------------------------------------------
  void Spectrum::selectTopK(unsigned int topK,
                            TwoValues<float> w,
                            Spectrum *putHere)
  {
	  WARN_MSG("selectTopK has not been updated for per-peak tolerances and should be replaced by rankFilterPeaks");
	  cerr << "selectTopK has not been updated for per-peak tolerances and should be replaced by rankFilterPeaks\n";

    vector<TwoValues<float> > &updatedPeakList = (putHere == 0) ? peakList
        : putHere->peakList;
    vector<TwoValues<float> > *oldPeakList;
    unsigned int peakIdx, idxLow = 0, idxHigh = 0, putIdx = 0;

    if (putHere == 0)
    {
      oldPeakList = new vector<TwoValues<float> > ;
      oldPeakList->resize(peakList.size());
      for (unsigned int i = 0; i < peakList.size(); i++)
        (*oldPeakList)[i] = peakList[i];
    }
    else
    {
      oldPeakList = &peakList; // No copying necessary if results goes into another object (putHere)
      (*putHere) = *this;
    }
    vector<TwoValues<float> > sortedInts;

    if (topK < oldPeakList->size())
    {
      for (peakIdx = 0; peakIdx < oldPeakList->size(); peakIdx++)
      {
        while (idxLow < peakIdx and (*oldPeakList)[idxLow][0]
            < (*oldPeakList)[peakIdx][0] + w[0])
          idxLow++;
        while (idxHigh < oldPeakList->size() and (*oldPeakList)[idxHigh][0]
            < (*oldPeakList)[peakIdx][0] + w[1])
          idxHigh++;

        if (topK < idxHigh - idxLow)
        {
          sortedInts.resize(idxHigh - idxLow);
          for (unsigned int neighIdx = idxLow; neighIdx < idxHigh; neighIdx++)
            sortedInts[neighIdx - idxLow].set((*oldPeakList)[neighIdx][1],
                                              neighIdx);
          sort(sortedInts.begin(), sortedInts.end());
          bool keep = false;
          for (unsigned int rank = 1; rank <= topK; rank++)
            if (sortedInts[sortedInts.size() - rank][1] == peakIdx)
            {
              keep = true;
              break;
            }
          if (not keep)
            continue;
        }

        updatedPeakList[putIdx++] = (*oldPeakList)[peakIdx];
      }
      updatedPeakList.resize(putIdx);
    }

    if (putHere == 0)
      delete oldPeakList;
  }

  // -------------------------------------------------------------------------
  void Spectrum::getSymmetricPeakPairs(float pmOffset, float tolerance, vector<
		vector<float> > &pairs, vector<vector<int> > &pairsIdx) {

	vector<float> tempTols;
	bool modTols = false;

	if (tolerance >= 0) {
		modTols = true;
		rememberTolerances();
		setPeakTolerance(tolerance);
	}

    float aaMass = parentMass + (pmOffset - idDist) * AAJumps::massHion;
    int curPair, // Index of the current pair being generated
        i, j, k; // iterators over the prefix/suffix peaks
    vector<float> peakScores; // Peak scores including scores of neigboring peaks within tolerance
    vector<bool> processed; // Keeps track of which peaks have already been used in some pair

    peakScores.resize(peakList.size());
    processed.resize(peakList.size());
    for (k = 0; k < (int)peakList.size(); k++)
    {
      processed[k] = false;
      peakScores[k] = peakList[k][1];
      // Removed 06/08/30 - It's reasonable but this is not the right place to do this
      //		for(i=k-1; i>=0 && peakList[i][0]>=peakList[k][0]-tolerance; i--) peakScores[k]+=peakList[i][1];
      //		for(i=k+1; i<peakList.size() && peakList[i][0]<=peakList[k][0]+tolerance; i++) peakScores[k]+=peakList[i][1];
    }

    //	for(k=0;k<peakList.size();k++) cerr<<k<<": "<<peakList[k][0]<<"\t"<<peakList[k][1]<<", "<<peakScores[k]<<"\n";

    unsigned int maxNumPairs = (unsigned int)round(((double)peakList.size()
        + 1.0) * (double)peakList.size() / 2);
    pairs.resize(maxNumPairs);
    pairsIdx.resize(maxNumPairs);
    //	pairs.resize(peakList.size()*(1+2*(int)ceil(tolerance/(0.1*idDist))));
    //	pairsIdx.resize(peakList.size()*(1+2*(int)ceil(tolerance/(0.1*idDist))));
    i = 0;
    j = peakList.size() - 1;
    curPair = 0;
    while (i <= j)
    {
      if (peakList[i][0] <= aaMass - peakList[j][0])
      {
        for (k = j; k > i && peakList[k][0] >= aaMass - peakList[i][0] - getTolerance(k) + getTolerance(i); k--)
        {
          pairs[curPair].resize(3);
          pairsIdx[curPair].resize(2);
          pairs[curPair][2] = peakScores[i] + peakScores[k];
          pairs[curPair][0] = peakList[i][0];
          pairs[curPair][1] = peakList[k][0];
          pairsIdx[curPair][0] = i;
          pairsIdx[curPair][1] = k;
          curPair++;
          processed[k] = true;
        }
        if (k == j && !processed[i])
        {
          pairs[curPair].resize(3);
          pairsIdx[curPair].resize(2);
          pairs[curPair][2] = peakScores[i];
          pairs[curPair][0] = peakList[i][0];
          pairs[curPair][1] = aaMass - peakList[i][0];
          pairsIdx[curPair][0] = i;
          pairsIdx[curPair][1] = -1;
          curPair++;
        }
        processed[i] = true;
        i++;
      }
      else
      {
        for (k = i; k < j && peakList[k][0] <= aaMass - peakList[j][0] + getTolerance(k) + getTolerance(j); k++)
        {
          pairs[curPair].resize(3);
          pairsIdx[curPair].resize(2);
          pairs[curPair][2] = peakScores[j] + peakScores[k];
          pairs[curPair][0] = peakList[k][0];
          pairs[curPair][1] = peakList[j][0];
          pairsIdx[curPair][0] = k;
          pairsIdx[curPair][1] = j;
          curPair++;
          processed[k] = true;
        }
        if (k == i && !processed[j])
        {
          pairs[curPair].resize(3);
          pairsIdx[curPair].resize(2);
          pairs[curPair][2] = peakScores[j];
          pairs[curPair][0] = aaMass - peakList[j][0];
          pairs[curPair][1] = peakList[j][0];
          pairsIdx[curPair][0] = -1;
          pairsIdx[curPair][1] = j;
          curPair++;
        }
        processed[j] = true;
        j--;
      }
    }
    pairs.resize(curPair);
    pairsIdx.resize(curPair);

    if (modTols) {
		revertTolerances();
	}
  }
  // -------------------------------------------------------------------------
  void Spectrum::makeSymmetric(float pmOffset, float tolerance,
		vector<int>* indicies) {

	bool modTols = false;

	if (tolerance >= 0) {
		modTols = true;
		rememberTolerances();
		setPeakTolerance(tolerance);
	}

	// Find all symmetric pairs in the spectrum
	vector<vector<float> > pairs;
	vector<vector<int> > pairsIdx;
	getSymmetricPeakPairs(pmOffset, -1.0, pairs, pairsIdx);

	map<TwoValues<float>*, int > indexMap;

	// Check which peaks already have a symmetric peak in the spectrum
	vector<bool> needsPair(size());
	unsigned int peakIdx;
	for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++) {
		needsPair[peakIdx] = true;
		if (indicies) {
			indexMap[&peakList[peakIdx]] = peakIdx;
		}
	}
	for (unsigned int pairIdx = 0; pairIdx < pairsIdx.size(); pairIdx++) {
		if (pairsIdx[pairIdx][0] >= 0 and pairsIdx[pairIdx][1] >= 0) {
			needsPair[pairsIdx[pairIdx][0]] = false;
			needsPair[pairsIdx[pairIdx][1]] = false;
		}
	}

	// Count how many new peaks are necessary
	unsigned int curCount = peakList.size(), newCount = 0;
	for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++) {
		newCount += (unsigned int) needsPair[peakIdx];
	}
	resize(curCount + newCount);
	if (modTols) {
		oldPeakTols.resize(curCount + newCount);
	}

	// Insert new symmetric peaks
	unsigned int newIdx = curCount;
	for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++)
		if (needsPair[peakIdx]) {
			peakList[newIdx].set(parentMass + (pmOffset - 1)
					* AAJumps::massHion - peakList[peakIdx][0],
					peakList[peakIdx][1]);
			peakTols[newIdx] = peakTols[peakIdx];
			if (indicies) {
				indexMap[&peakList[newIdx]] = -1;
			}
			if (modTols) {
				oldPeakTols[newIdx] = oldPeakTols[peakIdx];
			}
			if (peakList[newIdx][0] > 0) {
				newIdx++;
			}
		}
	resize(newIdx);

	if (modTols) {
		oldPeakTols.resize(newIdx);
		peakTols = oldPeakTols;
	}

	sortPeaks();

	if (indicies) {
		indicies->resize(size());
		for (int i = 0; i < size(); i++) {
			(*indicies)[i] = indexMap[&peakList[i]];
		}
	}
}
  // -------------------------------------------------------------------------
  bool Spectrum::compare(Spectrum &toSpec)
  {
    if (peakList.size() != toSpec.peakList.size() || parentMass
        != toSpec.parentMass)
      return false;
    for (int i = 0; i < peakList.size(); i++)
      if (peakList[i][0] != toSpec.peakList[i][0] || peakList[i][1]
          != toSpec.peakList[i][1] || peakTols[i] != toSpec.getTolerance(i))
        return false;
    return true;
  }
/*
  // -------------------------------------------------------------------------
  void Spectrum::mergeCommon(Spectrum &withSpectrum,
                             Spectrum *toSpec,
                             float shift,
                             float peakTol,
                             short mergeType)
  {
	  WARN_MSG("mergeCommon has not been updated for per-peak tolerances");
	  cerr << "mergeCommon has not been updated for per-peak tolerances\n";


    vector<int> idx1, idx2;
    unsigned int pivot;
    FindMatchPeaksAll2(*this, withSpectrum, shift, peakTol, idx1, idx2);
    Spectrum tmpToSpec;
    tmpToSpec.resize(idx1.size());
    tmpToSpec.copyNP(*this); // Temporary spectrum, in case toSpec==0

    float s, r1, r2;
    for (pivot = 0; pivot < idx1.size(); pivot++)
    {
      s = peakList[idx1[pivot]][1] + withSpectrum[idx2[pivot]][1];
      r1 = peakList[idx1[pivot]][1] / s;
      r2 = withSpectrum[idx2[pivot]][1] / s;
      tmpToSpec[pivot][0] = r1 * peakList[idx1[pivot]][0] + r2
          * (withSpectrum[idx2[pivot]][0] + shift);

      switch (mergeType)
      {
      case 0:
        tmpToSpec[pivot][1] = peakList[idx1[pivot]][1]
            + withSpectrum[idx2[pivot]][1];
        break;
      case 1:
        tmpToSpec[pivot][1] = max(peakList[idx1[pivot]][1],
                                  withSpectrum[idx2[pivot]][1]);
        break;
      case 2:
        tmpToSpec[pivot][1] = min(peakList[idx1[pivot]][1],
                                  withSpectrum[idx2[pivot]][1]);
        break;
      case 3:
        tmpToSpec[pivot][1] = peakList[idx1[pivot]][1];
        break;
      case 4:
        tmpToSpec[pivot][1] = withSpectrum[idx2[pivot]][1];
        break;
      }
    }
    sort(tmpToSpec.peakList.begin(), tmpToSpec.peakList.end());
    if (toSpec == (Spectrum *)0)
      (*this) = tmpToSpec;
    else
    {
      toSpec->copyNP(*this);
      (*toSpec) = tmpToSpec;
    }
  }
*/
  // -------------------------------------------------------------------------
  bool Spectrum::saveDTA(const char* outfile)
  {
    FILE* output = fopen(outfile, "w");
    fprintf(output, "%.3f", parentMass);
    fprintf(output, " ");
    fprintf(output, "%d", parentCharge);
    fprintf(output, "\n");
    for (int idx = 0; idx < peakList.size(); idx++)
    {
      fprintf(output, "%.3f", peakList[idx][0]);
      fprintf(output, " ");
      fprintf(output, "%.3f", peakList[idx][1]);
      fprintf(output, "\n");
    }
    fclose(output);
    return true;
  }
  // -------------------------------------------------------------------------
  float Spectrum::getTotalIonCurrent(void)
  {
    float total = 0;
    for (int i = 0; i < peakList.size(); i++)
    {
      total += peakList[i][1];
    }
    return total;
  }
  // -------------------------------------------------------------------------

  // MergeSpectra3 auxiliar variables:
  //static vector<vector<int> > gaux_peaksIdx(0); // Indices of selected peaks per spectrum
  //static vector<float> gaux_peaksScores(0); // Summed scores of the selected peaks per spectrum
  
}
