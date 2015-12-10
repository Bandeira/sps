#include "spectrum.h"
#include "alignment_scoring.h"
#include "label.h"
#include "batch.h"
#include "ms1.h"
#include "utils.h"

#include <errno.h>
#include <string>
#include <map>

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <list>

//const enum Spectrum::FragType { FragType_CID=0, FragType_ETD=1, FragType_HCD=2 };


//
// **************************************************************
//    Spectrum methods
// **************************************************************
//

Spectrum::Spectrum()
{
  //	std::cerr << "constructor" << endl;
  parentMass = 0;
  parentCharge = 0;
  scan = 0;
  msLevel = 0;
  msFragType = FragType_CID;
  peakList.resize(0);
  annotation.resize(0);
  ionTypes.resize(0);
  annotation_peptide.resize(0);
  resolution = 1;
  idDist = 1;
  info = (char *)0;
}

Spectrum::~Spectrum()
{
  //	std::cerr << "destructor" << endl;
  if (info != (char *)0)
  {
  }
  else
  {
    free(info);
    info = (char *)0;
  }
}

Spectrum::Spectrum(const Spectrum &other)
{
  //	std::cerr << "copy" << endl;
  parentMass = other.parentMass;
  parentCharge = other.parentCharge;
  scan = other.scan;
  resolution = other.resolution;
  idDist = other.idDist;
  peakList.resize(other.peakList.size());
  annotation.resize(other.annotation.size());
  annotation_peptide = other.annotation_peptide;
  ionTypes.resize(other.ionTypes.size());

  map<const ftIonFragment*, const ftIonFragment*> pointer_mapping; //mapping between old pointer (pointer_mapping.first) and new pointer (pointer_mapping.second

  for (unsigned int i = 0; i < other.peakList.size(); i++)
    peakList[i] = other.peakList[i];
  for (unsigned int i = 0; i < other.ionTypes.size(); i++)
  {
    ionTypes[i] = other.ionTypes[i];
    pointer_mapping[&(other.ionTypes[i])] = &(ionTypes[i]);
  }
  for (unsigned int i = 0; i < other.annotation.size(); i++)
  {
    //set old annotation fragment pointer to new annotation fragment pointer.
    const ftIonFragment* old_frag_ptr = other.annotation[i].first;
    annotation[i].first = pointer_mapping[old_frag_ptr];
    annotation[i].second = other.annotation[i].second;
  }
}

Spectrum &Spectrum::copyNP(const Spectrum &other)
{
  //	std::cerr << "copynp" << endl;
  parentMass = other.parentMass;
  parentCharge = other.parentCharge;
  scan = other.scan;
  msLevel = other.msLevel;
  msFragType = other.msFragType;
  resolution = other.resolution;
  idDist = other.idDist;
  return (*this);
}

Spectrum &Spectrum::operator=(const Spectrum &other)
{
  //	std::cerr << "assignment" << endl;
  parentMass = other.parentMass;
  parentCharge = other.parentCharge;
  scan = other.scan;
  resolution = other.resolution;
  idDist = other.idDist;
  peakList.resize(other.peakList.size());
  annotation.resize(other.annotation.size());
  annotation_peptide = other.annotation_peptide;
  ionTypes.resize(other.ionTypes.size());
  for (unsigned int i = 0; i < other.peakList.size(); i++)
    peakList[i] = other.peakList[i];

  map<const ftIonFragment*, const ftIonFragment*> pointer_mapping; //mapping between old pointer (pointer_mapping.first) and new pointer (pointer_mapping.second

  for (unsigned int i = 0; i < other.ionTypes.size(); i++)
  {
    ionTypes[i] = other.ionTypes[i];
    pointer_mapping[&(other.ionTypes[i])] = &(ionTypes[i]);
  }
  for (unsigned int i = 0; i < other.annotation.size(); i++)
  {
    //set old annotation fragment pointer to new annotation fragment pointer.
    const ftIonFragment* old_frag_ptr = other.annotation[i].first;
    annotation[i].first = pointer_mapping[old_frag_ptr];
    annotation[i].second = other.annotation[i].second;
  }
  return (*this);
}
//
//  Finds peaks in the spectrum inside a [-tolerance, tolerance] window around
//    targetMass. startIdx is a hint about where to start looking for targetMass
//    to avoid having to do binary search to find the peak closest to mass.
//    If the return set is empty then minIdx=-1 and maxIdx=-2;
//
inline void Spectrum::locateWithinTolerance(int startIdx,
                                            float targetMass,
                                            float tolerance,
                                            int &minIdx,
                                            int &maxIdx)
{
  minIdx = -1;
  maxIdx = -2;
  if (peakList[0][0] >= targetMass - tolerance)
    minIdx = 0;
  else
  {
    // Find first peak after targetMass
    for (minIdx = startIdx; minIdx < (int)peakList.size()
        && peakList[minIdx][0] < targetMass - tolerance; minIdx++)
      ;
    if (minIdx == (int)peakList.size())
    {
      minIdx = -1;
      return;
    }
    // Now find leftmost peak >= targetMass-tolerance
    for (; minIdx > 0 && peakList[minIdx - 1][0] >= targetMass - tolerance; minIdx--)
      ;
  }
  if (peakList[minIdx][0] > targetMass + tolerance)
  {
    minIdx = -1;
    return;
  }
  if (peakList[peakList.size() - 1][0] <= targetMass + tolerance)
    maxIdx = peakList.size() - 1;
  else
    for (maxIdx = minIdx; maxIdx < (int)peakList.size() - 1 && peakList[maxIdx
        + 1][0] <= targetMass + tolerance; maxIdx++)
      ;
}

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

void Spectrum::changeResolution(float newResolution, bool enforceRounding)
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
        }
      }
    }
    peakList.resize(j + 1);
  }
}

// addZPMpeaks - Ensures that the spectrum contains the ions b0/bk/y0/yk where k is the
//                 number of amino acids in the corresponding peptide. New peaks are added
//                 with a score of 0.1 if the b/y masses are not already present.
//
// tolerance - mass tolerance used when looking for b0/bk/y0/yk in the current spectrum
// ionOffest - used to specify b0/bk/y0/yk mass offset, depending on type of spectrum.
//               Use 0 for PRM spectra and AAJumps::massHion for MS/MS spectra
// includeY0k - If true then y0,yk are also added to the spectrum
// labels    - if !=NULL then the labels are shifted to match the correct peaks after any peak insertions
//
void Spectrum::addZPMpeaks(float tolerance,
                           float ionOffset,
                           bool includeY0k,
                           bool ctermH2O,
                           SpectrumPeakLabels *labels)
{
  float massB0 = ionOffset * idDist, massBk = parentMass - (ctermH2O
      ? AAJumps::massMH : AAJumps::massHion - ionOffset) * idDist, massY0 =
      (ctermH2O ? AAJumps::massH2O : 0 + ionOffset) * idDist, massYk =
      parentMass - (AAJumps::massHion - ionOffset) * idDist;

  list<TwoValues<float> > newPeaks;
  list<int> prevIndices;
  TwoValues<float> tmp(0, 0);
  unsigned int pivot;
  bool valIsPresent = false;
  for (pivot = 0; pivot < peakList.size() and peakList[pivot][0] <= massB0
      + tolerance; pivot++)
  {
    newPeaks.push_back(peakList[pivot]);
    prevIndices.push_back(pivot);
    valIsPresent = true;
  }
  if (!valIsPresent)
  {
    newPeaks.push_back(tmp.set(massB0, 0.1));
    prevIndices.push_back(-1);
  }
  valIsPresent = false;
  if (includeY0k)
  {
    for (; pivot < peakList.size() and peakList[pivot][0] < massY0 - tolerance; pivot++)
    {
      newPeaks.push_back(peakList[pivot]);
      prevIndices.push_back(pivot);
    }
    for (; pivot < peakList.size() and peakList[pivot][0] <= massY0 + tolerance; pivot++)
    {
      newPeaks.push_back(peakList[pivot]);
      prevIndices.push_back(pivot);
      valIsPresent = true;
    }
    if (!valIsPresent)
    {
      newPeaks.push_back(tmp.set(massY0, 0.1));
      prevIndices.push_back(-1);
    }
    valIsPresent = false;
  }
  else
    for (; pivot < peakList.size() and peakList[pivot][0] <= massY0 + tolerance; pivot++)
      ; // Skip Y0 masses
  for (; pivot < peakList.size() and peakList[pivot][0] < massBk - tolerance; pivot++)
  {
    newPeaks.push_back(peakList[pivot]);
    prevIndices.push_back(pivot);
  }
  for (; pivot < peakList.size() and peakList[pivot][0] <= massBk + tolerance; pivot++)
  {
    newPeaks.push_back(peakList[pivot]);
    prevIndices.push_back(pivot);
    valIsPresent = true;
  }
  if (!valIsPresent)
  {
    newPeaks.push_back(tmp.set(massBk, 0.1));
    prevIndices.push_back(-1);
  }
  valIsPresent = false;
  if (includeY0k)
  {
    for (; pivot < peakList.size() and peakList[pivot][0] < massYk - tolerance; pivot++)
    {
      newPeaks.push_back(peakList[pivot]);
      prevIndices.push_back(pivot);
    }
    for (; pivot < peakList.size() and peakList[pivot][0] <= massYk + tolerance; pivot++)
    {
      newPeaks.push_back(peakList[pivot]);
      prevIndices.push_back(pivot);
      valIsPresent = true;
    }
    if (!valIsPresent)
    {
      newPeaks.push_back(tmp.set(massYk, 0.1));
      prevIndices.push_back(-1);
    }
    valIsPresent = false;
  }

  if (labels != 0)
  {
    labels->resize(newPeaks.size());
    int labelsIdx = labels->size() - 1;
    for (list<int>::reverse_iterator iter = prevIndices.rbegin(); iter
        != prevIndices.rend(); iter++, labelsIdx--)
      if (*iter >= 0)
        (*labels)[labelsIdx] = (*labels)[*iter];
      else
        (*labels)[labelsIdx].set(0, 0, 0, 0, 1);
  }

  peakList.resize(newPeaks.size());
  list<TwoValues<float> >::iterator iter = newPeaks.begin();
  for (unsigned int i = 0; i < peakList.size(); i++, iter++)
    peakList[i] = *iter;
}

void Spectrum::maximizeZPMpeaks(float tolerance,
                                float ionOffset,
                                bool includeY0k)
{
  float massB0 = ionOffset * idDist, massBk = parentMass - (AAJumps::massMH
      - ionOffset) * idDist, massY0 = (AAJumps::massH2O + ionOffset) * idDist,
      massYk = parentMass - (AAJumps::massHion - ionOffset) * idDist;

  unsigned int pivot;
  float maxScore = 0;
  for (pivot = 0; pivot < peakList.size(); pivot++)
    if (peakList[pivot][1] > maxScore)
      maxScore = peakList[pivot][1];
  for (pivot = 0; pivot < peakList.size() and peakList[pivot][0] <= massB0
      + tolerance; pivot++)
    peakList[pivot][1] = maxScore;
  if (includeY0k)
  {
    for (; pivot < peakList.size() and peakList[pivot][0] < massY0 - tolerance; pivot++)
      ;
    for (; pivot < peakList.size() and peakList[pivot][0] <= massY0 + tolerance; pivot++)
      peakList[pivot][1] = maxScore;
  }
  for (; pivot < peakList.size() and peakList[pivot][0] < massBk; pivot++)
    ;
  for (; pivot < peakList.size() and peakList[pivot][0] <= massBk; pivot++)
    peakList[pivot][1] = maxScore;
  if (includeY0k)
  {
    for (; pivot < peakList.size() and peakList[pivot][0] < massYk - tolerance; pivot++)
      ;
    for (; pivot < peakList.size() and peakList[pivot][0] <= massYk + tolerance; pivot++)
      peakList[pivot][1] = maxScore;
  }
}

void Spectrum::normalize(float newTotalInt, bool removeNegatives)
{ // Normalizes total intensity to 100
  float totalIntensity = 0;
  for (unsigned int i = 0; i < peakList.size(); i++)
    totalIntensity += peakList[i][1];
  for (unsigned int i = 0; i < peakList.size(); i++)
    peakList[i][1] = newTotalInt * peakList[i][1] / totalIntensity;
}

void Spectrum::guessPrecursorZPM(Spectrum &parentMS1,
                                 float peakTol,
                                 short maxZ,
                                 IsoEnvelope &isoEnvs,
                                 bool strictMode)
{

  // Look for the best isotopic envelope match over all possible charges and monoisotopic masses
  float curMonoMass, chargeIncrement, bestPM = 0, curMatch, bestMatch = 1000000; // Match score of zero is optimal
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

//
// findMatches - Find all peaks within tolerance of baseMass
//
//  If startIdx is >=0 then the search starts at startIdx. Function returns the number of peaks within tolerance.
//
short Spectrum::findMatches(float baseMass,
                            float peakTol,
                            vector<int> &matchesIdx,
                            int startIdx)
{
  unsigned int baseIdx = (startIdx >= 0 and startIdx < peakList.size())
      ? (unsigned int)startIdx : 0, numPeaks = peakList.size(), matchCount = 0;
  matchesIdx.resize(numPeaks);
  for (; baseIdx > 0 and peakList[baseIdx][0] >= baseMass - peakTol; baseIdx--)
    ;
  for (; baseIdx < numPeaks and peakList[baseIdx][0] < baseMass - peakTol; baseIdx++)
    ; // Get to first peak above the lower end of the tolerance window
  for (; baseIdx < numPeaks and peakList[baseIdx][0] <= baseMass + peakTol; baseIdx++)
    matchesIdx[matchCount++] = baseIdx;
  matchesIdx.resize(matchCount);
  return matchCount;
}

//
//  findClosest - Finds the spectrum peak with peak mass closest to mass
//
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

void Spectrum::mergePeakList(vector<TwoValues<float> > &newPeaks,
                             Spectrum *putHere)
{
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

  if (putHere == 0)
    delete oldPeakList;
}

void Spectrum::mergePeakListRev(vector<TwoValues<float> > &newPeaks,
                                Spectrum *putHere)
{
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

  float aaMass = parentMass - 19;

  // Sorted merge of reversed peaks / endpoints
  updatedPeakList.resize(oldPeakList->size() + newPeaks.size());
  idxUpdated = 0;
  idxOld = oldPeakList->size() - 1;
  idxNew = newPeaks.size() - 1;
  while (idxOld >= 0 or idxNew >= 0)
  {
    if (idxOld < 0)
    {
      updatedPeakList[idxUpdated++].set(aaMass - newPeaks[idxNew][0],
                                        newPeaks[idxNew--][1]);
      continue;
    }
    if (idxNew < 0)
    {
      updatedPeakList[idxUpdated++].set(aaMass + 18 - (*oldPeakList)[idxOld][0],
                                        (*oldPeakList)[idxOld--][1]);
      continue;
    } // +18 because PRMs are assumed to be from y-ions
    if (fabs((aaMass + 18 - (*oldPeakList)[idxOld][0]) - (aaMass
        - newPeaks[idxNew][0])) < 0.0001)
    {
      updatedPeakList[idxUpdated++].set(max((float)0, aaMass + 18
          - (*oldPeakList)[idxOld][0]), (*oldPeakList)[idxOld--][1]
          + newPeaks[idxNew--][1]);
    }
    else if ((aaMass + 18 - (*oldPeakList)[idxOld][0]) < (aaMass
        - newPeaks[idxNew][0]))
    {
      updatedPeakList[idxUpdated++].set(max((float)0, aaMass + 18
          - (*oldPeakList)[idxOld][0]), (*oldPeakList)[idxOld--][1]);
    }
    else
    {
      updatedPeakList[idxUpdated++].set(max((float)0, aaMass
          - newPeaks[idxNew][0]), newPeaks[idxNew--][1]);
    }
  }
  updatedPeakList.resize(idxUpdated); // Just in case some peaks get merged along the way

  if (putHere == 0)
    delete oldPeakList;
}

// Converts a list of masses to a list of corresponding spectrum peak indices (closest mass)
//  masses is assumed to be sorted by increasing mass values.
//
void Spectrum::massesToIndices(vector<TwoValues<float> > &masses,
                               vector<int> &indices,
                               float peakTol)
{
  int idxClosest, idxPeaks = 0, idxMasses;
  float distClosest = peakTol + 1; // Supremum for distance to closest peak
  indices.resize(masses.size());
  for (idxMasses = 0; idxMasses < (int)masses.size(); idxMasses++)
  {
    while (idxPeaks > 0 and peakList[idxPeaks][0] > masses[idxMasses][0]
        - peakTol)
      idxPeaks--;
    while (idxPeaks < (int)peakList.size() and peakList[idxPeaks][0]
        < masses[idxMasses][0] - peakTol)
      idxPeaks++;
    idxClosest = -1;
    distClosest = peakTol + 1;
    while (idxPeaks < (int)peakList.size() and abs(peakList[idxPeaks][0]
        - masses[idxMasses][0]) <= peakTol + 0.0001)
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
}

void Spectrum::selectIndices(vector<int> &idx)
{
  unsigned int i;
  for (i = 0; i < idx.size(); i++)
    if (idx[i] < peakList.size())
      peakList[i] = peakList[idx[i]];
    else
      break;
  peakList.resize(i);
}

unsigned int Spectrum::filterLowInt(float minIntensity)
{
  unsigned int newIdx = 0, oldIdx;
  for (oldIdx = 0; oldIdx < peakList.size(); oldIdx++)
    if (peakList[oldIdx][1] >= minIntensity)
      peakList[newIdx++] = peakList[oldIdx];
  peakList.resize(newIdx);
  return newIdx;
}

void Spectrum::filterMasses(float minMass, float maxMass, float peakTol)
{
  unsigned int idxKeep, // iterator over retained peaks
      pivot; // iterator over unchecked peaks
  for (idxKeep = 0; idxKeep < peakList.size() and peakList[idxKeep][0]
      < minMass - peakTol; idxKeep++)
    ;
  if (maxMass < 0)
    maxMass = parentMass + peakTol;
  for (pivot = idxKeep; pivot < peakList.size(); pivot++)
    if (peakList[pivot][0] <= maxMass)
      peakList[idxKeep++] = peakList[pivot];
  peakList.resize(idxKeep);
}

void Spectrum::output(ostream &output)
{
  output << parentMass << " " << parentCharge << endl;
  for (unsigned int i = 0; i < peakList.size(); i++)
    output << peakList[i][0] << " " << peakList[i][1] << endl;
}

void Spectrum::output_ms2(ostream &output)
{
  //output << ":0.0.0\n";
  output << parentMass << " " << parentCharge << endl;
  for (unsigned int i = 0; i < peakList.size(); i++)
    output << peakList[i][0] << " " << peakList[i][1] << endl;
}

//
//  Computes the reversed spectrum and returns it in putHere (or reverses current if putHere==NULL)
//    The reverse of a peak with mass m is parentMass-1+pmOffset-m :
//    use pmOffset=0 for PRM spectra and pmOffset=2 for MS/MS spectra.
//
void Spectrum::reverse(float pmOffset, Spectrum *putHere)
{
  vector<TwoValues<float> > &updatedPeakList = (putHere == 0) ? peakList
      : putHere->peakList;
  vector<TwoValues<float> > *oldPeakList;

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

  float totMass = parentMass - AAJumps::massHion + pmOffset;
  unsigned int numPeaks = oldPeakList->size();
  for (unsigned int i = 0; i < numPeaks; i++)
    updatedPeakList[numPeaks - i - 1].set(totMass - (*oldPeakList)[i][0],
                                          (*oldPeakList)[i][1]);

  if (putHere == 0)
    delete oldPeakList;
}

// Adds offset to every peak mass and rotates negative masses to parentMass-cterm+mass
void Spectrum::rotate(float offset, float cterm)
{
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
      tmpPeaks[peakIdx].set(parentMass - cterm + peakList[peakIdx][0] + offset,
                            peakList[peakIdx][1]);
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

//
//  Selects the top k peaks in the spectrum.
//
void Spectrum::selectTopK(unsigned int topK, Spectrum *putHere)
{
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
          = (*oldPeakList)[(unsigned int)sortedInts[sortedInts.size() - peakIdx][1]];

    sort(updatedPeakList.begin(), updatedPeakList.end());
  }

  if (putHere == 0)
    delete oldPeakList;
}

//
//  A peak must be in the top k in a window of [w(0),w(1)] around its mass to be retained in the spectrum.
//
void Spectrum::selectTopK(unsigned int topK,
                          TwoValues<float> w,
                          Spectrum *putHere)
{
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

//  a) Rounds peak masses to the closest integer (after dividing by resolution)
//  b) Resulting peak intensities are added for peak with the same integer mass
//  c) Integer peak masses are converted back to float by multiplying by resolution
void Spectrum::roundMasses(float resolution)
{
  unsigned int peakIdx, pivot;

  // Step a)
  for (peakIdx = 0; peakIdx < peakList.size(); peakIdx++)
    peakList[peakIdx][0] = (float)round(peakList[peakIdx][0] / resolution);

  // Step b)
  for (peakIdx = 0, pivot = 1; pivot < peakList.size(); pivot++)
    if (fabs(peakList[peakIdx][0] - peakList[pivot][0]) < 0.00001)
      peakList[peakIdx][1] += peakList[pivot][1];
    else
      peakList[++peakIdx] = peakList[pivot];
  peakList.resize(peakIdx + 1);

  // Step c)
  for (peakIdx = 0; peakIdx < peakList.size(); peakIdx++)
    peakList[peakIdx][0] *= resolution;
}

//
//  Computes the sets of pairs in the spectrum and returns (in pairs) a list of
//    pairs sorted by minimum distance to the closest endpoint. Two peaks are
//    considered a pair if their masses add up to parentMass-idDist+pmOffset ->
//    use pmOffset=0 for PRM spectra and pmOffset=2 for MS/MS spectra.
//
//  pairs    - returns the masses of the peaks in the pairs (per pair)
//  pairsIdx - returns the indices of the peaks in the pairs (per pair)
//
void Spectrum::getPairs(float pmOffset,
                        float tolerance,
                        vector<vector<float> > &pairs,
                        vector<vector<int> > &pairsIdx)
{
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

  unsigned int maxNumPairs =
      (unsigned int)round(((double)peakList.size() + 1.0)
          * (double)peakList.size() / 2);
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
      for (k = j; k > i && peakList[k][0] >= aaMass - peakList[i][0] - 2
          * tolerance; k--)
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
      for (k = i; k < j && peakList[k][0] <= aaMass - peakList[j][0] + 2
          * tolerance; k++)
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
}

//
//  makeSymmetric - Forces a spectrum to be symmetric by adding symmetric peaks
//    whenever missing from the spectrum.
//
//  pmOffset - use 0 for PRM spectra, 2 for MS/MS spectra
//
void Spectrum::makeSymmetric(float pmOffset, float tolerance)
{
  // Find all symmetric pairs in the spectrum
  vector < vector<float> > pairs;
  vector < vector<int> > pairsIdx;
  getPairs(pmOffset, tolerance, pairs, pairsIdx);

  // Check which peaks already have a symmetric peak in the spectrum
  vector<bool> needsPair(peakList.size());
  unsigned int peakIdx;
  for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++)
    needsPair[peakIdx] = true;
  for (unsigned int pairIdx = 0; pairIdx < pairsIdx.size(); pairIdx++)
    if (pairsIdx[pairIdx][0] >= 0 and pairsIdx[pairIdx][1] >= 0)
    {
      needsPair[pairsIdx[pairIdx][0]] = false;
      needsPair[pairsIdx[pairIdx][1]] = false;
    }

  // Count how many new peaks are necessary
  unsigned int curCount = peakList.size(), newCount = 0;
  for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++)
    newCount += (unsigned int)needsPair[peakIdx];
  peakList.resize(curCount + newCount);

  // Insert new symmetric peaks
  unsigned int newIdx = curCount;
  for (peakIdx = 0; peakIdx < needsPair.size(); peakIdx++)
    if (needsPair[peakIdx])
    {
      peakList[newIdx].set(parentMass + (pmOffset - 1) * AAJumps::massHion
          - peakList[peakIdx][0], peakList[peakIdx][1]);
      if (peakList[newIdx][0] > 0)
        newIdx++;
    }
  peakList.resize(newIdx);
  sort(peakList.begin(), peakList.end());
}

bool Spectrum::compare(Spectrum &toSpec)
{
  if (peakList.size() != toSpec.peakList.size() || parentMass
      != toSpec.parentMass)
    return false;
  for (int i = 0; i < peakList.size(); i++)
    if (peakList[i][0] != toSpec.peakList[i][0] || peakList[i][1]
        != toSpec.peakList[i][1])
      return false;
  return true;
}

//
//  merge - merge the peaks in the current spectrum with those in withSpec and
//          output the result to toSpec.
//
//  scores1, scores2, scoresMerged are vectors with the same sizes as peakList,
//    withSpec.peakList and toSpec.peakList; scoresMerged is computed to add
//    the values in scores1,scores2 if the corresponding peaks with the same
//    indices (in peakList/withSpec.peakList) are merged into toSpec
//
void Spectrum::merge(vector<TwoValues<float> > &withPeaks,
                     Spectrum &toSpec,
                     vector<float> &scores1,
                     vector<float> &scores2,
                     vector<float> &scoresMerged)
{
  toSpec.copyNP(*this);
  toSpec.peakList.resize(peakList.size() + withPeaks.size());
  scoresMerged.resize(peakList.size() + withPeaks.size());
  unsigned int i = 0, j = 0, k = 0;
  while (k < toSpec.size() && (i < peakList.size() || j < withPeaks.size()))
  {
    if (i == peakList.size())
    {
      toSpec.peakList[k] = withPeaks[j];
      scoresMerged[k++] = scores2[j++];
      continue;
    }
    if (j == withPeaks.size())
    {
      toSpec.peakList[k] = peakList[i];
      scoresMerged[k++] = scores1[i++];
      continue;
    }
    if (peakList[i][0] == withPeaks[j][0])
    {
      toSpec.peakList[k].set(peakList[i][0], peakList[i][1] + withPeaks[j][1]);
      scoresMerged[k++] = scores1[i++] + scores2[j++];
    }
    else
    {
      if (peakList[i][0] < withPeaks[j][0])
      {
        toSpec.peakList[k] = peakList[i];
        scoresMerged[k++] = scores1[i++];
      }
      else
      {
        toSpec.peakList[k] = withPeaks[j];
        scoresMerged[k++] = scores2[j++];
      }
    }
  }
  toSpec.peakList.resize(k);
  scoresMerged.resize(k);
}

//
//  mergeCommon - merge the common peaks in the current spectrum with those in
//   withSpec and output the result to toSpec. Multiple peak matches (within tolerance)
//   are represented by different peaks in toSpec. Each peak in toSpec has its mass
//   given by the weighted average mass of the matched peaks. The intensity of
//   a merged peak is determined by the intensities of the paired peaks acording
//   to mergeType: 0 (sum), 1 (max), 2 (min), 3 (cur spec score only), 4 (withSpectrum score only)
//
void Spectrum::mergeCommon(Spectrum &withSpectrum,
                           Spectrum *toSpec,
                           float shift,
                           float peakTol,
                           short mergeType)
{
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

//
//  filterAdd - constructs a new spectrum (output) where each peak p in the current
//      object is added to output if it also exists in other (within tolerance).
//      The peak score in output is the score of p _plus_ the maximum score of any
//      matching peak in other _minus_ the absolute distance between them (in % of tolerance)
//
//  otherScores   - scores of the peaks selected from other
//  otherEPScores - summed scores of the peaks in other matched to the endpoints
//                   of the current spectrum (guaranteed matches)
//
void Spectrum::filterAdd(Spectrum &other,
                         float shift,
                         float tolerance,
                         Spectrum &output,
                         vector<float> &otherScores,
                         float &otherEPScores)
{
  int low = 0, high = 0; // Indices of the potential matching peaks in other
  int i, otherIdx, outputIdx = 0;

  output.copyNP(*this);

  // Check if m/z=0 in this spectrum matches any peaks in other
  /*	float maxScore=0;
   while(low<other.size() && other[low][0]+shift<-tolerance-0.00001) low++;
   while(high<other.size() && other[high][0]+shift<=tolerance+0.00001) high++; // high is index of first unreachable peak
   for(otherIdx=low+1; otherIdx<high; otherIdx++) maxScore=max(maxScore,other[otherIdx][1]);
   otherEPScores += maxScore;
   */
  // Check for matches of internal PRMs
  output.peakList.resize(peakList.size());
  otherScores.resize(peakList.size());
  for (i = 0; i < peakList.size(); i++)
  {
    while (low < other.size() && other[low][0] + shift < peakList[i][0]
        - tolerance - 0.00001)
      low++;
    while (high < other.size() && other[high][0] + shift <= peakList[i][0]
        + tolerance + 0.00001)
      high++; // high is index of first unreachable peak
    if (low == high)
      continue;
    output[outputIdx].set(peakList[i][0], peakList[i][1] + other[low][1]
        - fabs(peakList[i][0] - other[low][0] - shift) / tolerance);
    otherScores[outputIdx] = other[low][1];
    for (otherIdx = low + 1; otherIdx < high; otherIdx++)
      if (output[outputIdx][1] < peakList[i][1] + other[otherIdx][1]
          - fabs(peakList[i][0] - other[otherIdx][0] - shift) / tolerance)
      {
        output[outputIdx][1] = peakList[i][1] + other[otherIdx][1]
            - fabs(peakList[i][0] - other[otherIdx][0] - shift) / tolerance;
        otherScores[outputIdx] = other[otherIdx][1];
      }
    outputIdx++;
  }
  output.peakList.resize(outputIdx);
  otherScores.resize(outputIdx);

  // Check if m/z=parentMass-1 in this spectrum matches any peaks in other
  //
  //  Make idDist a member of Spectrum ??
  //
  /*	while(low<other.size() && other[low][0]+shift<parentMass-idDist-tolerance-0.00001) low++;
   while(high<other.size() && other[high][0]+shift<=parentMass-idDist+tolerance+0.00001) high++; // high is index of first unreachable peak
   for(maxScore=0, otherIdx=low+1; otherIdx<high; otherIdx++) maxScore=max(maxScore,other[otherIdx][1]);
   otherEPScores += maxScore;
   */
}

//
//  denovoLR - find a de-novo intepretation of the spectrum by finding the
//    highest scoring left-to-right path in the spectrum (does not exclude
//    forbidden pairs).
//
float Spectrum::denovoLR(AAJumps &jumps, float peakTol, vector<int> &matchedIdx)
{
  vector<TwoValues<float> > scores; // Pos[0] has current score, Pos[1] has best score so far
  vector<TwoValues<int> > indices; // Pos[0] has index of predecessor, Pos[1] has index of PRM with best path score so far
  int idxPeak, // Index of current peak
      idxJump, // Index of current jump
      idxPred; // Index of tentative predecessor

  scores.resize(peakList.size());
  indices.resize(peakList.size());
  for (idxPeak = 0; idxPeak < peakList.size(); idxPeak++)
  {
    scores[idxPeak].set(peakList[idxPeak][1], 0);
    indices[idxPeak].set(-1, idxPeak);
  }

  for (idxPeak = 0; idxPeak < peakList.size(); idxPeak++)
  {
    if (idxPeak == 0)
      scores[0].set(peakList[idxPeak][1], peakList[idxPeak][1]);
    else if (peakList[idxPeak][1] < scores[idxPeak - 1][1])
    {
      indices[idxPeak][1] = indices[idxPeak - 1][1];
      scores[idxPeak][1] = scores[idxPeak - 1][1];
    }

    for (idxJump = 0; idxJump < jumps.size() and jumps[idxJump] <= 0; idxJump++)
      ;
    for (idxPred = idxPeak - 1; idxPred >= 0 and peakList[idxPred][0]
        > peakList[idxPeak][0] - jumps[idxJump] + peakTol; idxPred--)
      ;
    //cerr<<"idxPred = "<<idxPred<<endl;
    while (idxPred >= 0 and idxJump < jumps.size())
    {
      for (; idxPred >= 0 and peakList[idxPred][0] > peakList[idxPeak][0]
          - jumps[idxJump] - peakTol; idxPred--)
      {
        //cerr<<"idxPeak = "<<idxPeak<<" considering pred "<<idxPred;
        if (scores[idxPred][0] + peakList[idxPeak][1] > scores[idxPeak][0])
        {
          scores[idxPeak][0] = scores[idxPred][0] + peakList[idxPeak][1];
          indices[idxPeak][0] = idxPred;
          //cerr << " and using its score.\n";
        } // else cerr << " and rejecting it!\n";
      }
      idxJump++;
      for (; idxJump < jumps.size() and idxPred >= 0 and peakList[idxPred][0]
          > peakList[idxPeak][0] - jumps[idxJump] + peakTol; idxPred--)
        ;
      //cerr<<"Looping for idxJump = "<<idxJump<<" out of max size "<<jumps.size()<<", idxPred = "<<idxPred<<endl;
    }
    //cerr<<"Out of the loop\n";

    // Check if path should link to highest scoring PRM with mass farther than the largest jump
    if (idxPred >= 0 and scores[idxPred][1] + peakList[idxPeak][1]
        > scores[idxPeak][0])
    {
      scores[idxPeak][0] = scores[idxPred][1] + peakList[idxPeak][1];
      indices[idxPeak][0] = indices[idxPred][1];
    }

    if (scores[idxPeak][0] > scores[idxPeak][1])
    {
      scores[idxPeak][1] = scores[idxPeak][0];
      indices[idxPeak][1] = idxPeak;
    }
    //cerr<<"idxPeak = "<<idxPeak<<", scores = ("<<scores[idxPeak][0]<<","<<scores[idxPeak][1]<<"), indices = ("<<indices[idxPeak][0]<<","<<indices[idxPeak][1]<<")\n";
  }

  // Fill matchedIdx
  int mIdx = 0;
  matchedIdx.resize(indices.size());
  idxPeak = indices[indices.size() - 1][1];
  while (idxPeak >= 0)
  {
    //cerr<<"matchedIdx["<<mIdx<<"] = "<<idxPeak<<endl;
    matchedIdx[mIdx++] = idxPeak;
    idxPeak = indices[idxPeak][0];
  }
  matchedIdx.resize(mIdx);

  return scores[scores.size() - 1][1];
}

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
/**
 * Helper function for annotate
 */
bool compare_probs(const ftIonFragment &i, const ftIonFragment &j)
{
  return (i.prob < j.prob);
}

/** helper function for annotate. Runs through match vector and
 *sets annotation vector
 *@param matches - vector of ion indices (from 0) to set matches for
 *@param annotation - vector of fragments to set
 *@param ionIdx - index of ion 0..N-1 where N is length of peptide. i.e. for b1, ionIdx = 1
 *@param currIonFrag - pointer to current ftIonFragment we're currently annotating
 */
void set_annotation_to_matches(vector<int> &matches,
                               vector<pair<const ftIonFragment*, short> > &annotation,
                               int ionIdx,
                               const ftIonFragment* currIonFrag)
{
  for (int i = 0; i < matches.size(); i++)
  {
    short annot_index = ionIdx + 1;
    const ftIonFragment* last_ion = annotation[matches[i]].first;
    if (last_ion && last_ion->prob >= currIonFrag->prob)
    {
      // if the last annotation had greater probabiliy, ignore current annotation
    }
    else
    {
      annotation[matches[i]].first = currIonFrag;
      annotation[matches[i]].second = annot_index;
    }
  }
}

/**
 * helper function for annotate
 * sets ionTypes based on inclusion list
 * @param ionNamesInclude -  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
 * then just include all ions in MS2Model. ex. "y,b,y++,b++"
 * @param inputIonTypes - definition of ion types: mass offsets, ion probabilities, prefix/suffix.
 * is copied into ionTypes
 * @param ionTypes - output vector for fragments
 */

bool copy_ion_types(string &ionNamesInclude,
                    vector<ftIonFragment> &inputIonTypes,
                    vector<ftIonFragment> &outputIonTypes)
{
  if (ionNamesInclude.compare("all") == 0)
  { //if we're including all ion types
    outputIonTypes.resize(inputIonTypes.size());
    copy(inputIonTypes.begin(), inputIonTypes.end(), outputIonTypes.begin());
    return inputIonTypes.size() == outputIonTypes.size(); // check to make sure our vector sizes match
  }
  else
  {
    vector < string > ionNames;

    splitText(ionNamesInclude.c_str(), ionNames, (const char*)",");

    for (int i = 0; i < inputIonTypes.size(); i++)
    {
      for (int j = 0; j < ionNames.size(); j++)
      {
        if (ionNames[j].compare(inputIonTypes[i].name) == 0)
        { //if current ion matches include list.
          ftIonFragment ion = inputIonTypes[i];
          outputIonTypes.push_back(ion); //copy ftIonFragment object.
          break;
        }
      }
    }
    return ionNames.size() == outputIonTypes.size();
  }
}

/**
 * add annotations for matched peaks from MS2ScoringModel
 * @param peptide amino acid sequence used to determine peak annotations
 * @param ionNamesInclude comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
 * then just include all ions in MS2Model. ex. "y,b,y++,b++"
 * @param inputIonTypes definition of ion types: mass offsets, ion probabilities, prefix/suffix.
 * is copied into ionTypes
 * @param prmOffset sets applied to theoretical ion masses to locate Prefix Residue Masses
 * @param srmOffset sets applied to theoretical ion masses to locate Suffix Residue Masses
 * @param peakTol - mass tolerance (in Da) when matching ion masses to peak masses
 */

void Spectrum::annotate(string &peptide,
                        string &ionNamesInclude,
                        MS2ScoringModel &inputIonTypes,
                        float prmOffset,
                        float srmOffset,
                        float peakTol)
{
  //clear existing annotations
  annotation.clear();
  ionTypes.clear();

  //copy ftIonFragments into Spectrum ionTypes
  if (!copy_ion_types(ionNamesInclude, inputIonTypes.probs, ionTypes))
  {
    cerr << "Unable to copy ionNames from MS2ScoringModel! ionNamesInclude "
        << ionNamesInclude << endl;
    return;
  }

  //resize annotation to contain annotations for each peak
  annotation.resize(peakList.size());

  //initialize annotate vector to be null values
  for (int i = 0; i < annotation.size(); i++)
  {
    annotation[i].first = (ftIonFragment*)NULL;
    annotation[i].second = 0;
  }

  //generate srm and prm masses, no offsets.
  vector<float> prm_masses;
  vector<float> srm_masses;

  //generate total peptide mass (summed value of aa masses)
  float peptide_mass;

  AAJumps jumps(1);

  jumps.getPRMandSRMMasses(peptide, prm_masses, srm_masses, peptide_mass);

  short peptide_length = prm_masses.size(); //length of peptide

  sort(ionTypes.begin(), ionTypes.end(), compare_probs); //sort fragment types in ascending order so that low probability
  //annotations will be overwritten.

  short ionIdx;
  int last_srm_index; //keep track of last srm index so other srm searches start from there
  int last_prm_index; //keep track of last prm index so other prm searches start from there

  for (ionIdx = 0; ionIdx < peptide_length - 1; ionIdx++)
  {

    vector<ftIonFragment>::const_iterator currIonFrag;

    for (currIonFrag = ionTypes.begin(); currIonFrag != ionTypes.end(); currIonFrag++)
    {
      vector<int> matches;
      if (currIonFrag->isIF) //we're looking at internal fragment
      {
        if (currIonFrag->isNTerm)
        {
          float curr_mass = (prm_masses[ionIdx] + currIonFrag->massOffset
              + prmOffset) / currIonFrag->charge;
#ifdef DEBUG
          std::cout << "mass " << curr_mass << endl;
          std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif

          findMatches(curr_mass, peakTol, matches, last_prm_index);
          if (matches.size() > 0)
            last_prm_index = matches[matches.size() - 1];
        }
        else
        {
          float curr_mass = (srm_masses[ionIdx] + currIonFrag->massOffset
              + srmOffset) / currIonFrag->charge;
#ifdef DEBUG
          std::cout << "mass " << curr_mass << endl;
          std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif
          findMatches(curr_mass, peakTol, matches, last_srm_index);
          if (matches.size() > 0)
            last_srm_index = matches[matches.size() - 1];
        }
      }
      set_annotation_to_matches(matches, annotation, ionIdx, &(*currIonFrag));
    }
  }

  int last_pm_index; //keep track of last pm index so other searches start from there.

  vector<ftIonFragment>::const_iterator currIonFrag;

  // looking at peptide mass shift (i.e., total peptide, total peptide - water, etc.)
  for (currIonFrag = ionTypes.begin(); currIonFrag != ionTypes.end(); currIonFrag++)
  {
    vector<int> matches;

    if (!currIonFrag->isIF)
    { // we're looking at peptpde mass shift
      float curr_mass = (peptide_mass + currIonFrag->massOffset + prmOffset
          + srmOffset) / currIonFrag->charge;
      findMatches(curr_mass, peakTol, matches, last_pm_index);

      if (matches.size() > 0)
        last_pm_index = matches[matches.size() - 1];
    }
    set_annotation_to_matches(matches, annotation, ionIdx, &(*currIonFrag));
  }

  //set annotation string
  annotation_peptide = peptide;
}

//
// **************************************************************
//    SpecSet methods
// **************************************************************
//

vector<list<int> > &SpecSet::getMassesHistogram(vector<list<int> > &results,
                                                float resolution)
{
  int i;
  double t, massIdx = 0; // Value is always integer but max is (double,double)

  for (i = 0; i < specs.size(); i++)
  {
    t = round(specs[i].parentMass / resolution);
    if (massIdx < t)
      massIdx = t;
  }
  results.resize((int)massIdx + 1);

  for (i = 0; i < specs.size(); i++)
    results[(int)(round(specs[i].parentMass / resolution))].push_front(i);
  return results;
}

template<class T> unsigned int SpecSet::extract(vector<T> &features,
                                                T featureValue,
                                                SpecSet &output)
{
  output.resize(max(features.size(), specs.size()));
  unsigned int keptSpectra = 0, pivot;

  for (pivot = 0; pivot < output.size(); pivot++)
    if (features[pivot] == featureValue)
      output[keptSpectra++] = specs[pivot];

  output.resize(keptSpectra);
  return keptSpectra;
}

void instantiate_SpecSet_extract()
{
  SpecSet specs, specsOut;

  vector<int> features;
  specs.extract(features, (int)0, specsOut);
}

bool SpecSet::LoadSpecSet_pkl_mic(const char* filename)
{
  BufferedLineReader blr;
  TwoValues<float> peak;
  Spectrum loadSpec;
  specs.resize(0);
  if (blr.Load(filename) <= 0 || blr.size() < 2)
  {
    cerr << "ERROR: Not enough lines loaded\n";
    return false;
  }
  bool inspec = false;
  for (unsigned int i = 0; i < blr.size(); i++)
  {

    vector<string> peak_vals;
    if (blr.getline(i) == NULL)
      break;
    const char* next_line = blr.getline(i);
    if (!splitText(next_line, peak_vals, "\t "))
    {
      cerr << "ATTENTION: stopped reading peaks before eof in " << filename
          << "\n";
      break;
    }
    if (peak_vals.size() == 0)
    {
      if (inspec)
      {
        sort(loadSpec.peakList.begin(), loadSpec.peakList.end());
        specs.push_back(loadSpec);
      }
      inspec = false;
      continue;
    }
    if (peak_vals.size() < 2)
    {
      cerr << "ATTENTION: found wierd line after header: <" << next_line
          << "> ---- skippping ...\n";
      continue;
    }

    if (!inspec && peak_vals.size() == 3)
    {
      loadSpec.peakList.resize(0);
      loadSpec.parentMass = (float)atof(peak_vals[0].c_str());
      loadSpec.parentMZ = (double)strtod(peak_vals[0].c_str(), NULL);
      loadSpec.parentCharge = (short)atoi(peak_vals[2].c_str());
      if (loadSpec.parentCharge > 0)
      {
        loadSpec.parentMass = loadSpec.parentMass
            * ((float)loadSpec.parentCharge) - ((float)loadSpec.parentCharge)
            + AAJumps::massHion;
      }
      inspec = true;
      continue;
    }

    if (inspec && peak_vals.size() < 2)
    {
      cerr << "ATTENTION: not enough peak info after header: <" << next_line
          << "> ---- skippping ...\n";
      continue;
    }

    if (inspec)
    {
      peak[0] = (float)atof(peak_vals[0].c_str());
      peak[1] = (float)atof(peak_vals[1].c_str());
      loadSpec.peakList.push_back(peak);
    }
  }
  if (inspec)
  {
    sort(loadSpec.peakList.begin(), loadSpec.peakList.end());
    specs.push_back(loadSpec);
  }
  return true;
}

/**
 * Get scan at scan number
 * @param scan_num - scan number (indexed from 1!)
 */
Spectrum* SpecSet::getScan(int scan_num)
{
  //cheating a bit, check to see whether scans are in order and we can just jump to the correct
  //scan
  Spectrum *curr_scan = &(specs[scan_num - 1]);

  if (curr_scan->scan == scan_num)
  {
    return curr_scan;
  }
  else
  {
    vector<Spectrum>::iterator spec_iter;
    for (spec_iter = specs.begin(); spec_iter != specs.end(); spec_iter++)
    {
      if (scan_num == spec_iter->scan)
      {
        return &(*spec_iter);
      }
    }
  }
  return (Spectrum*) NULL;
}

//
//  LoadSpecSet_mgf: Loads a set of spectra from a .mgf file. Recognizes and processes
//    these MGF tags: BEGIN IONS, END IONS, CHARGE, PEPMASS
//
//  Note: If a spectrum has more than 1 header then this function only uses the last header.
//
unsigned int SpecSet::LoadSpecSet_mgf(const char *filename)
{
  BufferedLineReader blr;
  resize(0);
  if (blr.Load(filename) <= 0 or blr.size() < 3)
    return 0; // A valid file must have at least 3 lines: BEGIN_IONS, END_IONS and one mass/intensity peak

  unsigned int lineIdx, specIdx, peakIdx;

  // Counts number of spectra and number of peaks per spectrum
  list<int> peaksPerSpec;
  int first;
  int numPeaks = 0; // Counts number of peaks in the spectrum
  for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
  {
    if (numPeaks >= 0 and strcmp("END IONS", blr.getline(lineIdx)) == 0)
    {
      peaksPerSpec.push_back(numPeaks);
      numPeaks = -1;
      continue;
    }
    if (strcmp("BEGIN IONS", blr.getline(lineIdx)) == 0)
    {
      numPeaks = 0;
      continue;
    }
    first = (int)blr.getline(lineIdx)[0];
    if (numPeaks >= 0 and first >= 48 and first <= 57)
      numPeaks++;
  }

  // Parse spectra
  resize(peaksPerSpec.size());
  lineIdx = 0;
  char *token, *line;
  for (specIdx = 0; specIdx < specs.size(); specIdx++)
  {
    // Skip empty lines
    while (lineIdx < blr.size() and (blr.getline(lineIdx)[0] == 0
        or strcmp("BEGIN IONS", blr.getline(lineIdx)) != 0))
      lineIdx++;
    if (lineIdx == blr.size())
    {
      cerr << "Error loading " << filename << " - " << specIdx
          << " spectra instead of " << specs.size() << "?\n";
      resize(0);
      return 0;
    }

    // Start of spectrum
    if (strcmp("BEGIN IONS", blr.getline(lineIdx)) != 0)
    {
      cerr << "ERROR: Expected BEGIN IONS, found '" << blr.getline(lineIdx)
          << "' (line " << lineIdx + 1 << ")\n";
      return 0;
    }
    else
      lineIdx++;

    // Read peaks/charge/parent mass
    specs[specIdx].resize(peaksPerSpec.front());
    peaksPerSpec.pop_front();
    peakIdx = 0;
    specs[specIdx].parentCharge = 0;
    while (lineIdx < blr.size() and strcmp("END IONS", blr.getline(lineIdx))
        != 0)
    {
      line = blr.getline(lineIdx++);
      if (line[0] >= 48 and line[0] <= 57)
      {
        token = strtok(line, " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse peak mass on line " << lineIdx << "!\n";
          resize(0);
          return 0;
        }
        specs[specIdx][peakIdx][0] = atof(token);
        token = strtok(NULL, " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse peak intensity on line " << lineIdx
              << "!\n";
          resize(0);
          return 0;
        }
        specs[specIdx][peakIdx][1] = atof(token);
        peakIdx++;
        continue;
      }

      if (strncmp("CHARGE=+", line, 8) == 0)
      {
        specs[specIdx].parentCharge = (short)atof(&line[8]);
        continue;
      }
      else if (strncmp("CHARGE=", line, 7) == 0)
      {
        specs[specIdx].parentCharge = (short)atof(&line[7]);
        continue;
      }

      if (strncmp("TITLE=Scan Number: ", line, 19) == 0)
        specs[specIdx].scan = (unsigned int)atof(&line[19]);
      if (strncmp("SCAN=", line, sizeof("SCAN=") - 1) == 0)
        specs[specIdx].scan = (unsigned int)atof(line + sizeof("SCAN=") - 1);
      if (strncmp("PEPMASS=", line, 8) == 0)
      {
        specs[specIdx].parentMass = atof(&line[8]);
        specs[specIdx].parentMZ = (double)strtod(&line[8], NULL);
      }
    }
    if (specs[specIdx].parentCharge > 0)
      specs[specIdx].parentMass = (specs[specIdx].parentMass
          * specs[specIdx].parentCharge) - (AAJumps::massHion
          * (specs[specIdx].parentCharge - 1));
    if (specs[specIdx].scan == 0)
    {
      specs[specIdx].scan = specIdx + 1;
    }
    lineIdx++;
  }

  return specs.size();
}

//
//  LoadSpecSet_ms2: Loads a set of spectra from a .ms2 file. Spectra must be separated
//    by at least one empty line.
//
//  Note: If a spectrum has more than 1 header then this function only uses the last header.
//
unsigned int SpecSet::LoadSpecSet_ms2(const char *filename)
{
  BufferedLineReader blr;
  resize(0);
  if (blr.Load(filename) <= 0 or blr.size() < 3)
    return 0; // A valid file must have at least 3 lines: 1 header (2 lines) + 1 (m/z,intensity) pair

  unsigned int lineIdx, specIdx, peakIdx;

  // Counts number of spectra and number of peaks per spectrum
  list<int> peaksPerSpec;
  int numPeaks = 1; // Counts number of peaks in the spectrum
  for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    if (blr.getline(lineIdx)[0] == ':')
    {
      if (numPeaks > 0)
        peaksPerSpec.push_back(numPeaks);
      numPeaks = -1;
    }
    else if (blr.getline(lineIdx)[0] != 0)
      numPeaks++;
  peaksPerSpec.push_back(numPeaks); // Number of peaks in the last spectrum
  peaksPerSpec.pop_front(); // First element is just the '1' used to initialize numPeaks

  // Parse spectra
  resize(peaksPerSpec.size());
  lineIdx = 0;
  char *token;
  for (specIdx = 0; specIdx < specs.size(); specIdx++)
  {
    // Skip empty lines
    while (blr.getline(lineIdx)[0] == 0 and lineIdx < blr.size())
      lineIdx++;
    if (lineIdx == blr.size())
    {
      cerr << "Error loading " << filename << " - " << specIdx
          << " spectra instead of " << specs.size() << "?\n";
      resize(0);
      return 0;
    }

    // Parse header(s)
    while (blr.getline(lineIdx)[0] == ':')
    {
      lineIdx++;
      token = strtok(blr.getline(lineIdx++), " \t");
      if (!token)
      {
        cerr << "Error loading " << filename
            << " - could not parse parent mass for spectrum " << specIdx + 1
            << "?\n";
        resize(0);
        return 0;
      }
      specs[specIdx].parentMass = (float)atof(token);
      specs[specIdx].parentMZ = (double)strtod(token, NULL);
      token = strtok(NULL, " \t");
      if (!token)
      {
        cerr << "Error loading " << filename
            << " - could not parse parent charge for spectrum " << specIdx + 1
            << "?\n";
        resize(0);
        return 0;
      }
      specs[specIdx].parentCharge = (short)atof(token);
      if (lineIdx == blr.size())
      {
        cerr << "Error loading " << filename << " - " << specIdx
            << " spectra instead of " << specs.size() << "?\n";
        resize(0);
        return 0;
      }
    }

    // Read spectrum peaks
    specs[specIdx].resize(peaksPerSpec.front());
    peaksPerSpec.pop_front();
    for (peakIdx = 0; peakIdx < specs[specIdx].size(); peakIdx++)
    {
      token = strtok(blr.getline(lineIdx++), " \t");
      if (!token)
      {
        cerr << "Error loading " << filename
            << " - could not parse peak mass for spectrum " << specIdx + 1
            << ", peak " << peakIdx + 1 << "?\n";
        resize(0);
        return 0;
      }
      specs[specIdx][peakIdx][0] = (float)atof(token);
      token = strtok(NULL, " \t");
      if (!token)
      {
        cerr << "Error loading " << filename
            << " - could not parse peak intensity for spectrum " << specIdx + 1
            << ", peak " << peakIdx + 1 << "?\n";
        resize(0);
        return 0;
      }
      specs[specIdx][peakIdx][1] = (float)atof(token);
      if (lineIdx == blr.size() and peakIdx < specs[specIdx].size() - 1)
      {
        cerr << "Error loading " << filename
            << " - end of file before end of spectrum " << specIdx + 1
            << " ended, got only " << peakIdx + 1 << " peaks instead of "
            << specs.size() << "?\n";
        resize(0);
        return 0;
      }
    }
  }

  return specs.size();
}

//
//  LoadSpecSet_prms: Loads a set of spectra from a .prms file, as output but the
//     current version of pepnovo_prms (2007/04/09). Text file with multiple spectra,
//     each delimited in the following format:
//
//  - ">> <original msms file index> <scan/index # in msms file> <single-spectrum file name>"
//      -- scan/index # in msms file is stored in the Spectrum.scan field.
//      -- The whole info line is stored in the Spectrum.info field.
//  - Set of peaks in "<peak mass> <peak intensity/score>" format
//  - Spectrum terminated by empty line
//  - Last pair of mass/intensity corresponds to the corrected sum-of-amino-acid-masses parent mass
//
//  Note: If a spectrum has more than 1 header then this function only uses the last header.
//
unsigned int SpecSet::LoadSpecSet_prms(const char *filename)
{
  unsigned int specIdx, lineIdx, peakIdx;
  list<unsigned int> specSizes;
  char *line;
  BufferedLineReader blr;

  // Load whole file into memory
  if (blr.Load(filename) <= 0)
    return 0;

  // Count # spectra in the file
  bool inSpectrum = false;
  peakIdx = 0;
  for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
  {
    line = blr.getline(lineIdx);
    if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
    {
      inSpectrum = true;
      peakIdx = 0;
      continue;
    }
    if (line[0] == 0)
    {
      if (inSpectrum)
        specSizes.push_back(peakIdx);
      inSpectrum = false;
    }
    else if (inSpectrum)
      peakIdx++;
  }
  if (inSpectrum)
    specSizes.push_back(peakIdx); // In case there is no empty line after the last spectrum
  specs.resize(specSizes.size());

  // Parse the text
  peakIdx = 0;
  specIdx = 0;
  list<unsigned int>::iterator sizesIter = specSizes.begin();
  char *token;
  for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
  {
    line = blr.getline(lineIdx);

    if (line[0] == 0)
    {
      if (inSpectrum)
      { // End of current spectrum
        if (peakIdx > 0 and peakIdx <= specs[specIdx].size())
          specs[specIdx].parentMass = specs[specIdx][peakIdx - 1][0]
              + AAJumps::massMH;
        specIdx++;
        inSpectrum = false;
      }
      continue;
    }

    // Check for start of a new spectrum
    if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
    {
      specs[specIdx].info = (char *)malloc(strlen(line) - 1);
      strcpy(specs[specIdx].info, &line[2]);
      char *tok = strtok(line, " \t");
      tok = strtok(NULL, " \t"); // Skip ">> " and "<file_index> "
      tok = strtok(NULL, " \t");
      specs[specIdx].scan = (unsigned int)strtoul(tok, NULL, 10);
      specs[specIdx].resize(*sizesIter);
      sizesIter++;
      peakIdx = 0;
      inSpectrum = true;
      continue;
    }

    // Peak <mass> <intensity> pair
    token = strtok(line, " \t");
    if (token[0] == 0)
    {
      cerr << "ERROR reading peak mass for peak " << peakIdx
          << " for the spectrum entitled (" << specs[specIdx].info
          << ") in file " << filename << "!\n";
      specs.resize(0);
      return 0;
    }
    specs[specIdx][peakIdx][0] = atof(token);
    token = strtok(NULL, " \t");
    if (token[0] == 0)
    {
      cerr << "ERROR reading peak intensity for peak " << peakIdx
          << " for the spectrum entitled (" << specs[specIdx].info
          << ") in file " << filename << "!\n";
      specs.resize(0);
      return 0;
    }
    specs[specIdx][peakIdx][1] = atof(token);
    peakIdx++;
  }

  return specs.size();
}

unsigned int SpecSet::LoadSpecSet_prmsv3(const char *filename)
{
  unsigned int specIdx, lineIdx, peakIdx;
  list<unsigned int> specSizes;
  char *line;
  BufferedLineReader blr;

  // Load whole file into memory
  if (blr.Load(filename) <= 0)
    return 0;

  // Count # spectra in the file
  bool inSpectrum = false;
  peakIdx = 0;
  for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
  {
    line = blr.getline(lineIdx);
    if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
    {
      inSpectrum = true;
      peakIdx = 0;
      continue;
    }
    if (line[0] != 0 and line[0] == '#' and inSpectrum)
      inSpectrum = false;
    if (line[0] == 0)
    {
      if (inSpectrum)
        specSizes.push_back(peakIdx);
      inSpectrum = false;
    }
    else if (inSpectrum and line[0] != 'C')
      peakIdx++;
  }
  if (inSpectrum)
    specSizes.push_back(peakIdx); // In case there is no empty line after the last spectrum
  specs.resize(specSizes.size());

  // Parse the text
  peakIdx = 0;
  specIdx = 0;
  list<unsigned int>::iterator sizesIter = specSizes.begin();
  char *token;
  for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
  { // Skip header lines
    line = blr.getline(lineIdx);
    if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      break;
  }

  for (; lineIdx < blr.size(); lineIdx++)
  {
    line = blr.getline(lineIdx);

    if (line[0] == 0)
    {
      if (inSpectrum)
      { // End of current spectrum
        //  -- Not valid for Pepnovo v3
        //				if(peakIdx>0 and peakIdx<=specs[specIdx].size()) specs[specIdx].parentMass = specs[specIdx][peakIdx-1][0]+AAJumps::massMH;
        specIdx++;
        inSpectrum = false;
      }
      continue;
    }

    // Check for start of a new spectrum
    if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
    {
      specs[specIdx].info = (char *)malloc(strlen(line) - 1);
      strcpy(specs[specIdx].info, &line[2]);
      char *tok = strtok(line, " \t");
      tok = strtok(NULL, " \t"); // Skip ">> " and "<file_index> "
      tok = strtok(NULL, " \t");
      specs[specIdx].scan = (unsigned int)strtoul(tok, NULL, 10);

      // Read charge/mass header
      line = blr.getline(++lineIdx);
      if (line[0] != 0 and line[0] == '#')
      {
        inSpectrum = false;
        continue;
      }
      specs[specIdx].parentCharge = (short)atoi(&line[8]);
      specs[specIdx].parentMass = (float)atof(&line[19]);
      peakIdx = 0;
      inSpectrum = true;

      specs[specIdx].resize(*sizesIter);
      sizesIter++;
      continue;
    }

    // Peak <mass> <intensity> pair
    token = strtok(line, " \t");
    if (token[0] == 0)
    {
      cerr << "ERROR reading peak mass for peak " << peakIdx
          << " for the spectrum entitled (" << specs[specIdx].info
          << ") in file " << filename << "!\n";
      specs.resize(0);
      return 0;
    }
    specs[specIdx][peakIdx][0] = atof(token);
    token = strtok(NULL, " \t");
    if (token[0] == 0)
    {
      cerr << "ERROR reading peak intensity for peak " << peakIdx
          << " for the spectrum entitled (" << specs[specIdx].info
          << ") in file " << filename << "!\n";
      specs.resize(0);
      return 0;
    }
    specs[specIdx][peakIdx][1] = atof(token);
    peakIdx++;
  }

  return specs.size();
}

unsigned int SpecSet::LoadSpecSet_pkl(const char *filename)
{
  int numSpecs = 1; // Assume file contains at least one spectrum and count one additional spectrum per blank line
  ifstream input(filename);
  if (!input.is_open() || !input.good())
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    specs.resize(0);
    return 0;
  }

  list<int> specSizes; // Register number of peaks per spectrum

  char *lineBuffer = (char *)malloc(1025);
  //    char *lineBuffer = new char [1024];
  int numPeaks = -1; // zero only after reading first tuple of (parent mass, intensity, charge)
  while (!input.eof() && !input.fail())
  {
    input.getline(lineBuffer, 1024, '\n');
    if (lineBuffer[0] == 0 or lineBuffer[0] == '\n' or lineBuffer[0] == '\r')
    {
      if (numPeaks >= 0)
      {
        numSpecs++;
        specSizes.push_back(numPeaks);
        numPeaks = -1;
      }
    }
    else
    {
      if (((int)lineBuffer[0] >= (int)'0' and (int)lineBuffer[0] <= (int)'9')
          or lineBuffer[0] == '-')
        numPeaks++;
    }
  }
  if (numPeaks >= 0)
  {
    specSizes.push_back(numPeaks);
  }
  else
  {
    numSpecs--; // In case there was no empty line at the end of the pkl file
  }

  specs.resize(numSpecs);
  //    input.seekg(0,ios_base::beg);
  //    input.close();  input.open(filename);
  //    Neither of the above worked, so ...
  input.close();
  ifstream input2(filename);

  float foo;
  for (int i = 0; i < numSpecs; i++)
  {
    input2 >> specs[i].parentMass >> foo >> specs[i].parentCharge; // Read precursor mass, total intensity, charge state
    if (specs[i].parentCharge > 0)
      specs[i].parentMass = specs[i].parentMass * specs[i].parentCharge
          - specs[i].parentCharge + AAJumps::massHion; // Convert to PM+19

    //        specs[i].numPeaks = specSizes.front();   specSizes.pop_front();
    //        specs[i].peakList.resize(specs[i].numPeaks);
    specs[i].peakList.resize(specSizes.front());
    specSizes.pop_front();

    for (int j = 0; j < specs[i].size(); j++)
      input2 >> specs[i][j][0] >> specs[i][j][1];
    //            input2 >> specs[i].peakList[j][0] >> specs[i].peakList[j][1];
    specs[i].scan = i + 1;
    input2.getline(lineBuffer, 1024, '\n'); // Read intermediate newline
  }

  //    delete lineBuffer;
  input2.close();
  free(lineBuffer);
  return 1;
}

short SpecSet::SaveSpecSet_info(const char *filename)
{
  ofstream output(filename);
  if (!output)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  for (int i = 0; i < specs.size(); i++)
    output << specs[i].info << "\n";

  output.close();
  return 1;
}

short SpecSet::SaveSpecSet_pkl(const char *filename)
{
  ofstream output(filename);
  if (!output)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  for (int i = 0; i < specs.size(); i++)
  {
    output << specs[i].parentMass << " -1 " << specs[i].parentCharge << "\n";
    for (int j = 0; j < specs[i].size(); j++)
    {
      output << specs[i][j][0] << " " << specs[i][j][1] << endl;
      //            output << specs[i].peakList[j][0] << " " << specs[i].peakList[j][1] << endl;
    }
    if (i < specs.size() - 1)
    {
      output << endl;
    }
  }

  output.close();
  return 1;
}

bool SpecSet::SaveAnnotations(const char* filename)
{
  FILE* out_buf = fopen(filename, "w");

  if (out_buf == NULL)
    return false;

  for (int i = 0; i < size(); i++)
  {
    if (specs[i].annotation_peptide.length() == 0)
      fprintf(out_buf, "\n");
    else
      fprintf(out_buf,
              "%s\n",
              specs[i].annotation_peptide.substr(2,
                                                 specs[i].annotation_peptide.length()
                                                     - 4).c_str());
  }
  fclose(out_buf);
  return true;
}

short SpecSet::SaveSpecSet_ms2(const char* filename)
{
  for (int i = 0; i < specs.size(); i++)
  {
    ostringstream outs;
    outs << filename << i << ".ms2";
    ofstream output(outs.str().c_str());
    if (!output)
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      return -1;
    }
    specs[i].output_ms2(output);
    output.close();
  }
  return 1;
}

short SpecSet::SaveSpecSet_mgf(const char* filename, const char* activation)
{
  FILE* output = fopen(filename, "w");
  if (!output)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  float mz, z;
  for (int i = 0; i < specs.size(); i++)
  {
    z = (float)specs[i].parentCharge;
    mz = specs[i].parentMass;
    if (z > 0)
      mz = (mz + (z - 1) * 1.0072763) / z;
    fprintf(output, "BEGIN IONS\nPEPMASS=%.5f\nCHARGE=+%.0f\n", mz, z);
    //output << "BEGIN IONS\nPEPMASS="<<mz<<"\nCHARGE=+"<<z<<"\n";
    if (specs[i].annotation_peptide.length() > 0)
    {
      fprintf(output, "PEPTIDE=%s\n", specs[i].annotation_peptide.c_str());
    }
    if (activation != 0)
    {
      fprintf(output, "ACTIVATION=%s\n", activation);
      //output<<"ACTIVATION="<<activation<<"\n";
    }
    if (specs[i].scan > 0)
    {
      fprintf(output, "TITLE=Scan Number: %d\n", specs[i].scan);
      //output<<"TITLE=Scan Number: "<<specs[i].scan<<"\n";
    }
    for (int j = 0; j < specs[i].size(); j++)
    {
      fprintf(output, "%.5f %.5f\n", specs[i][j][0], specs[i][j][1]);
      //output << specs[i][j][0] << " " << specs[i][j][1] << endl;
      //            output << specs[i].peakList[j][0] << " " << specs[i].peakList[j][1] << endl;
    }
    fprintf(output, "END IONS\n");
    //output << "END IONS\n";
    if (i < specs.size() - 1)
    {
      fprintf(output, "\n");
      //output << endl;
    }
  }
  fclose(output);
  //output.close();
  return 1;
}

//
//  SaveSpecSet_pklbin - saves a set of spectra in binary format. File format
//    1 int - number of spectra in the file
//    numSpecs shorts - number of peaks per spectrum in the file
//    arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
//
short SpecSet::SaveSpecSet_pklbin(const char *filename)
{
  FILE *fp;
  unsigned int numSpecs = specs.size();
  unsigned int i, p;

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numSpecs, sizeof(int), 1, fp); // Number of spectra in the file

  unsigned short *numPeaks = (unsigned short *)malloc(sizeof(short) * numSpecs);
  unsigned short maxNumPeaks = 0;
  for (i = 0; i < numSpecs; i++)
  {
    numPeaks[i] = specs[i].size();
    maxNumPeaks = max(maxNumPeaks, numPeaks[i]);
  }
  fwrite(numPeaks, sizeof(short), numSpecs, fp); // Number of peaks per spectrum in the file

  float *peaksBuffer = (float *)malloc(2 * (maxNumPeaks + 1) * sizeof(float));
  unsigned int pbIdx;
  for (i = 0; i < numSpecs; i++)
  {
    peaksBuffer[0] = specs[i].parentMass;
    peaksBuffer[1] = (float)specs[i].parentCharge;
    for (pbIdx = 2, p = 0; p < numPeaks[i]; p++)
    {
      peaksBuffer[pbIdx++] = specs[i][p][0];
      peaksBuffer[pbIdx++] = specs[i][p][1];
    }
    fwrite(peaksBuffer, sizeof(float), 2 * (numPeaks[i] + 1), fp); // [parentMass charge] followed by [masses intensities]
  }

  free(peaksBuffer);
  free(numPeaks);
  fclose(fp);
  return 1;
}

//  SaveSpecSet_pklbin - saves a set of spectra in binary format. File format
//    1 int  - Always set to zero (to differentiate from pklbin 1.0)
//    1 byte - sub-version of pklbin2
//    1 int  - number of spectra in the file
//    numSpecs shorts - number of peaks per spectrum in the file
//    arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
//
short SpecSet::SaveSpecSet_pklbin2(const char *filename)
{
  FILE *fp;
  unsigned int numSpecs = specs.size();
  unsigned int i, p;

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  i = 0;
  fwrite(&i, sizeof(int), 1, fp); // Empty 4 bytes
  char tmp = 0;
  fwrite(&tmp, 1, 1, fp); // Sub-version of pklbin2
  fwrite(&numSpecs, sizeof(int), 1, fp); // Number of spectra in the file

  unsigned int *numPeaks = (unsigned int *)malloc(sizeof(unsigned int)
      * numSpecs);
  unsigned int maxNumPeaks = 0;
  for (i = 0; i < numSpecs; i++)
  {
    numPeaks[i] = specs[i].size();
    maxNumPeaks = max(maxNumPeaks, numPeaks[i]);
  }
  fwrite(numPeaks, sizeof(unsigned int), numSpecs, fp); // Number of peaks per spectrum in the file

  float *peaksBuffer = (float *)malloc(2 * (maxNumPeaks + 1) * sizeof(float));
  unsigned int pbIdx;
  for (i = 0; i < numSpecs; i++)
  {
    peaksBuffer[0] = specs[i].parentMass;
    peaksBuffer[1] = (float)specs[i].parentCharge;
    for (pbIdx = 2, p = 0; p < numPeaks[i]; p++)
    {
      peaksBuffer[pbIdx++] = specs[i][p][0];
      peaksBuffer[pbIdx++] = specs[i][p][1];
    }
    fwrite(peaksBuffer, sizeof(float), 2 * (numPeaks[i] + 1), fp); // [parentMass charge] followed by [masses intensities]
  }

  free(peaksBuffer);
  free(numPeaks);
  fclose(fp);
  return 1;
}

//
//  LoadSpecSet_pklbin - loads a set of spectra in binary format. File format
//    1 int - number of spectra in the file
//    numSpecs shorts - number of peaks per spectrum in the file
//    arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
//
unsigned int SpecSet::LoadSpecSet_pklbin(const char *filename)
{
  FILE *fp;
  unsigned int numSpecs = 0;
  unsigned short *numPeaks; // Pointer to array containing 1+number of peaks per spectrum
  float *data; // Pointer to array containing spectrum data
  unsigned int dataSize = 0; // Size of the data buffer
  unsigned int i, p, dataIdx, count, numValues;

  fp = fopen(filename, "r");
  if (fp == 0)
  {
    cerr << "ERROR: " << strerror(errno) << " '" << filename << "'\n";
    return 0;
  }

  count = fread(&numSpecs, sizeof(unsigned int), 1, fp); // Number of spectra in the file
  if (count != 1)
  {
    cerr << "ERROR reading " << filename << "!\n";
    return 0;
  }

  numPeaks = (unsigned short *)malloc(sizeof(unsigned short) * numSpecs); // Number of peaks per spectrum
  if (numPeaks == (unsigned short *)0)
  {
    cerr << "ERROR: Not enough memory for " << numSpecs << " spectra.\n";
    fclose(fp);
    return 0;
  }
  count = fread(numPeaks, sizeof(unsigned short), numSpecs, fp);
  if (count != numSpecs)
  {
    cerr << "ERROR reading " << filename << "!\n";
    free(numPeaks);
    return 0;
  }
  for (i = 0; i < numSpecs; i++)
  {
    if (numPeaks[i] > dataSize)
      dataSize = numPeaks[i];
  }
  dataSize = 2 * dataSize + 2; // Each spectrum 'line' has 2 values and there's one additional 'line' with parent mass/charge
  data = (float *)malloc(sizeof(float) * dataSize);

  specs.resize(numSpecs);

  for (i = 0; i < numSpecs; i++)
  {
    numValues = 2 * numPeaks[i] + 2;
    count = fread(data, sizeof(float), numValues, fp);
    if (count != numValues)
    {
      cerr << "ERROR reading " << filename << "!\n";
      free(numPeaks);
      free(data);
      return 0;
    }

    specs[i].parentMass = data[0];
    specs[i].parentCharge = (int)data[1];

    if ( specs[i].parentCharge > 0 ) {
      specs[i].parentMZ = (specs[i].parentMass
          + ((specs[i].parentCharge - 1.0) * AAJumps::massHion))
          / specs[i].parentCharge;
    }
    else
    {
      specs[i].parentMZ = specs[i].parentMass;
    }

    specs[i].scan = i + 1;
    specs[i].resize(numPeaks[i]);
    for (p = 0, dataIdx = 2; dataIdx < numValues; dataIdx += 2, p++)
      specs[i][p].set(data[dataIdx], data[dataIdx + 1]);
  }

  free(numPeaks);
  free(data);
  fclose(fp);
  return 1;
}

TwoValues<int> countConsecutiveMP(Spectrum& contig1,
                                  Spectrum& contig2,
                                  float shiftFor,
                                  float shiftRev,
                                  float pmTol,
                                  float minOverlapArea,
                                  bool debug)
{
  vector<int> idx1_f, idx2_f, idx1_r, idx2_r;
  Spectrum rev2;
  contig2.reverse(0.0 - AAJumps::massH2O, &rev2);
  FindMatchPeaksAll2(contig1, contig2, shiftFor, pmTol, idx1_f, idx2_f);
  FindMatchPeaksAll2(contig1, rev2, shiftRev, pmTol, idx1_r, idx2_r);

  vector<bool> c1_f(contig1.size()), c1_r(contig1.size()),
      c2_f(contig2.size()), c2_r(contig2.size());

  for (int i = 0; i < idx1_f.size(); i++)
  {
    c1_f[idx1_f[i]] = true;
  }
  for (int i = 0; i < idx1_r.size(); i++)
  {
    c1_r[idx1_r[i]] = true;
  }
  for (int i = 0; i < idx2_f.size(); i++)
  {
    c2_f[idx2_f[i]] = true;
  }
  for (int i = 0; i < idx2_r.size(); i++)
  {
    c2_r[idx2_r[i]] = true;
  }

  int cons_peaks = 0, max_cons_f1 = 0, max_cons_f2 = 0, max_cons_r1 = 0,
      max_cons_r2 = 0;
  for (int i = 0; i < c1_f.size(); i++)
  {
    if (c1_f[i])
    {
      cons_peaks++;
      if (cons_peaks > max_cons_f1)
      {
        max_cons_f1 = cons_peaks;
      }
    }
    else
    {
      cons_peaks = 0;
    }
  }
  cons_peaks = 0;
  for (int i = 0; i < c1_r.size(); i++)
  {
    if (c1_r[i])
    {
      cons_peaks++;
      if (cons_peaks > max_cons_r1)
      {
        max_cons_r1 = cons_peaks;
      }
    }
    else
    {
      cons_peaks = 0;
    }
  }
  cons_peaks = 0;
  for (int i = 0; i < c2_f.size(); i++)
  {
    if (c2_f[i])
    {
      cons_peaks++;
      if (cons_peaks > max_cons_f2)
      {
        max_cons_f2 = cons_peaks;
      }
    }
    else
    {
      cons_peaks = 0;
    }
  }
  cons_peaks = 0;
  for (int i = 0; i < c2_r.size(); i++)
  {
    if (c2_r[i])
    {
      cons_peaks++;
      if (cons_peaks > max_cons_r2)
      {
        max_cons_r2 = cons_peaks;
      }
    }
    else
    {
      cons_peaks = 0;
    }
  }
  max_cons_f1 = min(max_cons_f1, max_cons_f2);
  max_cons_r1 = min(max_cons_r1, max_cons_r2);
  return TwoValues<int> (max_cons_f1, max_cons_r2);
}

TwoValues<float> getContigOverlapScores(Spectrum& contig1,
                                        Spectrum& contig2,
                                        float shiftFor,
                                        float shiftRev,
                                        float pmTol,
                                        float minOverlapArea,
                                        bool debug)
{
  Spectrum contig = contig1;
  Spectrum f_contig_o = contig2;

  Spectrum r_contig_o;
  contig2.reverse(0.0 - AAJumps::massH2O, &r_contig_o);

  Spectrum f_contig = f_contig_o;
  Spectrum r_contig = r_contig_o;
  unsigned int idxUse = 0;
  for (unsigned int i = 0; i < f_contig_o.size(); i++)
  {
    float mass = f_contig_o[i][0] + shiftFor;
    if (mass > 0.0 - pmTol && mass < contig.parentMass + pmTol)
    {
      f_contig.peakList[idxUse] = f_contig_o.peakList[i];
      idxUse++;
    }
  }
  f_contig.resize(idxUse);
  idxUse = 0;
  for (unsigned int i = 0; i < r_contig_o.size(); i++)
  {
    float mass = r_contig_o[i][0] + shiftRev;
    if (mass > 0.0 - pmTol && mass < contig.parentMass + pmTol)
    {
      r_contig.peakList[idxUse] = r_contig_o.peakList[i];
      idxUse++;
    }
  }
  r_contig.resize(idxUse);
  idxUse = 0;

  if (debug)
  {
    cout << "Forward Shift: " << shiftFor << "\nReverse Shift: " << shiftRev
        << "\nSpectrum:\n";
    contig.output(cout);
    cout << "\nForward Contig2:\n";
    for (unsigned int i = 0; i < f_contig.size(); i++)
    {
      cout << f_contig[i][0] << " + " << shiftFor << " = " << f_contig[i][0]
          + shiftFor << "\n";
    }
    cout << "\nReverse Contig2:\n";
    for (unsigned int i = 0; i < r_contig.size(); i++)
    {
      cout << r_contig[i][0] << " + " << shiftRev << " = " << r_contig[i][0]
          + shiftRev << "\n";
    }
  }

  Spectrum f_merged;
  Spectrum f_merged2;
  Spectrum r_merged;
  Spectrum r_merged2;

  if (f_contig.size() > 0)
  {
    contig.mergeCommon(f_contig, &f_merged, shiftFor, pmTol, 3);
    contig.mergeCommon(f_contig, &f_merged2, shiftFor, pmTol, 4);
  }
  if (r_contig.size() > 0)
  {
    contig.mergeCommon(r_contig, &r_merged, shiftRev, pmTol, 3);
    contig.mergeCommon(r_contig, &r_merged2, shiftRev, pmTol, 4);
  }

  float f_2_min = 0.0 - shiftFor - pmTol;
  float f_2_max = contig1.peakList.back()[0] - shiftFor + pmTol;
  float r_2_min = 0.0 - shiftRev - pmTol;
  float r_2_max = contig1.peakList.back()[0] - shiftRev + pmTol;
  float f_1_min = shiftFor - pmTol;
  float f_1_max = shiftFor + contig2.peakList.back()[0] + pmTol;
  float r_1_min = shiftRev - pmTol;
  float r_1_max = shiftRev + contig2.peakList.back()[0] + pmTol;

  if (debug)
  {
    cout << "pmTol = " << pmTol << endl;
    cout << "f_1_min = " << f_1_min << ", f_1_max = " << f_1_max << endl;
    cout << "f_2_min = " << f_2_min << ", f_2_max = " << f_2_max << endl;
    cout << "r_1_min = " << r_1_min << ", r_1_max = " << r_1_max << endl;
    cout << "r_2_min = " << r_2_min << ", r_2_max = " << r_2_max << endl;

    cout << "\nMerged Forward 1: \n";
    f_merged.output(cout);
    cout << "\nMerged Forward 2: \n";
    f_merged2.output(cout);
    cout << "\nMerged Reverse 1: \n";
    r_merged.output(cout);
    cout << "\nMerged Reverse 2: \n";
    r_merged2.output(cout);
    /*
     cout << "\nContig1  parent mass: " << contig.parentMass << "\nContig2 parent mass: " << f_contig.parentMass << " = " << r_contig.parentMass << "\n";
     cout << "Forward 1 min - max: " << f_1_min << " " << f_1_max << "\n";
     cout << "Reverse 1 min - max: " << r_1_min << " " << r_1_max << "\n";
     cout << "Forward 2 min - max: " << f_2_min << " " << f_2_max << "\n";
     cout << "Reverse 2 min - max: " << r_2_min << " " << r_2_max << "\n\n\n";
     */
  }

  float f_intensity = 0, f_intensity2 = 0;
  for (unsigned int i = 0; i < contig.size(); i++)
  {
    if (contig[i][0] > f_1_min && contig[i][0] < f_1_max)
    {
      f_intensity += contig[i][1];
    }
  }
  for (unsigned int i = 0; i < f_contig.size(); i++)
  {
    if (f_contig[i][0] > f_2_min && f_contig[i][0] < f_2_max)
    {
      f_intensity2 += f_contig[i][1];
    }
  }

  float r_intensity = 0, r_intensity2 = 0;
  for (unsigned int i = 0; i < contig.size(); i++)
  {
    if (contig[i][0] > r_1_min && contig[i][0] < r_1_max)
    {
      r_intensity += contig[i][1];
    }
  }
  for (unsigned int i = 0; i < r_contig.size(); i++)
  {
    if (r_contig[i][0] > r_2_min && r_contig[i][0] < r_2_max)
    {
      r_intensity2 += r_contig[i][1];
    }
  }

  float m_f_intensity = 0;
  for (unsigned int i = 0; i < f_merged.size(); i++)
    m_f_intensity += f_merged[i][1];
  float m_r_intensity = 0;
  for (unsigned int i = 0; i < r_merged.size(); i++)
    m_r_intensity += r_merged[i][1];
  float m_f_intensity2 = 0;
  for (unsigned int i = 0; i < f_merged2.size(); i++)
    m_f_intensity2 += f_merged2[i][1];
  float m_r_intensity2 = 0;
  for (unsigned int i = 0; i < r_merged2.size(); i++)
    m_r_intensity2 += r_merged2[i][1];

  float minF = min((m_f_intensity / f_intensity), (m_f_intensity2
      / f_intensity2));
  float minR = min((m_r_intensity / r_intensity), (m_r_intensity2
      / r_intensity2));

  if (debug)
  {
    cout << "m_f_intensity = " << m_f_intensity << ", m_f_intensity2 = "
        << m_f_intensity2 << endl;
    cout << "f_intensity = " << f_intensity << ", f_intensity2 = "
        << f_intensity2 << endl;
    cout << "m_r_intensity = " << m_r_intensity << ", m_r_intensity2 = "
        << m_r_intensity2 << endl;
    cout << "r_intensity = " << r_intensity << ", r_intensity2 = "
        << r_intensity2 << endl;
    cout << "minF = " << minF << ", minR = " << minR << endl;
  }

  if (min(f_1_max, contig.parentMass) - max((float)0, f_1_min) < minOverlapArea)
  {
    minF = -1.0;
  }
  else if (f_merged.size() == 0)
  {
    minF = 0;
  }

  if (min(r_1_max, contig.parentMass) - max((float)0, r_1_min) < minOverlapArea)
  {
    minR = -1.0;
  }
  else if (r_merged.size() == 0)
  {
    minR = 0;
  }

  return TwoValues<float> (minF, minR);
}

TwoValues<int> getContigOverlapPeaks(Spectrum& contig1,
                                     Spectrum& contig2,
                                     float shiftFor,
                                     float shiftRev,
                                     float pmTol,
                                     float minOverlapArea,
                                     bool debug)
{
  vector<int> idx1_f, idx2_f, idx1_r, idx2_r;
  Spectrum rev2;
  contig2.reverse(0.0 - AAJumps::massH2O, &rev2);
  FindMatchPeaksAll2(contig1, contig2, shiftFor, pmTol, idx1_f, idx2_f);
  FindMatchPeaksAll2(contig1, rev2, shiftRev, pmTol, idx1_r, idx2_r);

  vector<bool> c1_f(contig1.size()), c1_r(contig1.size()),
      c2_f(contig2.size()), c2_r(contig2.size());

  for (int i = 0; i < idx1_f.size(); i++)
  {
    c1_f[idx1_f[i]] = true;
  }
  for (int i = 0; i < idx1_r.size(); i++)
  {
    c1_r[idx1_r[i]] = true;
  }
  for (int i = 0; i < idx2_f.size(); i++)
  {
    c2_f[idx2_f[i]] = true;
  }
  for (int i = 0; i < idx2_r.size(); i++)
  {
    c2_r[idx2_r[i]] = true;
  }

  int max_f1 = 0, max_f2 = 0, max_r1 = 0, max_r2 = 0;
  for (int i = 0; i < c1_f.size(); i++)
  {
    if (c1_f[i])
    {
      max_f1++;
    }
  }
  for (int i = 0; i < c1_r.size(); i++)
  {
    if (c1_r[i])
    {
      max_r1++;
    }
  }
  for (int i = 0; i < c2_f.size(); i++)
  {
    if (c2_f[i])
    {
      max_f2++;
    }
  }
  for (int i = 0; i < c2_r.size(); i++)
  {
    if (c2_r[i])
    {
      max_r2++;
    }
  }
  max_f1 = min(max_f1, max_f2);
  max_r1 = min(max_r1, max_r2);
  return TwoValues<int> (max_f1, max_r1);
}

//
//  MakeSymmetric - convert a list of masses/scores to a symmetric list of masses/scores,
//    i.e. outSpec[i]+outSpec[n-i-1]=parentMass+18, adding entries with score zero if necessary.
//  Notes:
//    1 - if outSpec==NULL or outSpec==inSpec then the resulting contents replace those of inSpec
//    2 - if indices!=NULL then each position contains the inSpec index of the entry in outSpec (or -1 if not in inSpec)
//    3 - Current implementation is tailored to PRM spectra. Set parameter parentMass to sum of aa masses+18+3
//          to apply this function to MS/MS spectra.
//
void MakeSymmetric(vector<TwoValues<float> > *inSpec,
                   float parentMass,
                   float peakTol,
                   vector<TwoValues<float> > *outSpec,
                   vector<int> *indices)
{
  bool replaceContents = (inSpec == outSpec || outSpec == 0);
  if (replaceContents)
    outSpec = new vector<TwoValues<float> > ;

  outSpec->resize(2 * inSpec->size());
  if (indices)
    indices->resize(2 * inSpec->size());

  int inBottom = 0, inTop = inSpec->size() - 1, // Used to traverse inSpec
      outBottom = 0, outTop = outSpec->size() - 1; // Used to fill in the symmetric outSpec
  float symMass;
  //cerr << "[inBottom,mass], [inTop,mass], [outBottom,mass], [outTop,mass]\n";
  while (inBottom <= inTop)
  {
    //cerr << "["<<inBottom<<","<<(*inSpec)[inBottom][0]<<"], ["<<inTop<<","<<(*inSpec)[inTop][0]<<"], "; cerr.flush();
    if (abs((*inSpec)[inBottom][0] - (parentMass - AAJumps::massHion
        - (*inSpec)[inTop][0])) <= peakTol)
    { // symmetric peaks
      (*outSpec)[outBottom++] = (*inSpec)[inBottom++];
      if (indices)
        (*indices)[outBottom - 1] = inBottom - 1;
      if (inBottom <= inTop) // if condition avoids the case where a peak is its own symmetric
      {
        (*outSpec)[outTop--] = (*inSpec)[inTop--];
        if (indices)
          (*indices)[outTop + 1] = inTop + 1;
      }
    }
    else if ((*inSpec)[inBottom][0] < (parentMass - AAJumps::massHion
        - (*inSpec)[inTop][0]))
    {
      symMass = parentMass - AAJumps::massHion - (*inSpec)[inBottom][0];
      (*outSpec)[outBottom++] = (*inSpec)[inBottom++];
      (*outSpec)[outTop--].set(symMass, 0);
      if (indices)
      {
        (*indices)[outBottom - 1] = inBottom - 1;
        (*indices)[outTop + 1] = -1;
      }
    }
    else
    {
      symMass = parentMass - AAJumps::massHion - (*inSpec)[inTop][0];
      (*outSpec)[outBottom++].set(symMass, 0);
      (*outSpec)[outTop--] = (*inSpec)[inTop--];
      if (indices)
      {
        (*indices)[outBottom - 1] = -1;
        (*indices)[outTop + 1] = inTop + 1;
      }
    }
    //cerr << "["<<outBottom-1<<","<<(*outSpec)[outBottom-1][0]<<"], ["<<outTop+1<<","<<(*outSpec)[outTop+1][0]<<"]\n";
  }

  // Remove the possible empty positions in the middle of the array (by shifting the top elements down)
  int difference = outTop + 1 - outBottom;
  if (difference > 0)
    for (int pivot = outTop + 1; pivot < outSpec->size(); pivot++)
    {
      (*outSpec)[pivot - difference] = (*outSpec)[pivot];
      if (indices)
        (*indices)[pivot - difference] = (*indices)[pivot];
    }
  outSpec->resize(outSpec->size() - difference);
  if (indices)
    indices->resize(outSpec->size() - difference);

  if (replaceContents)
  {
    inSpec->resize(outSpec->size());
    for (int i = 0; i < outSpec->size(); i++)
      (*inSpec)[i] = (*outSpec)[i];
    delete outSpec;
  }
}

// MergeSpectra3 auxiliar variables:
static vector<vector<int> > gaux_peaksIdx(0); // Indices of selected peaks per spectrum
static vector<float> gaux_peaksScores(0); // Summed scores of the selected peaks per spectrum

//
//  MergeSpectra3 - merges the 3 lists of peaks (after applying the shifts) and
//      adds the scores of all peaks with the same peak masses.
//
//  minPeakCount - A peak is reported in the merged spectrum if it occurs in at least
//                   minPeakCount spectra (1,2 or 3).
//
unsigned int MergeSpectra3(Spectrum &spec1,
                           Spectrum &spec2,
                           Spectrum &spec3,
                           float shift12,
                           float shift13,
                           float peakTol,
                           int minPeakCount,
                           Spectrum &mergedSpec,
                           vector<vector<int> > &peaksIndices)
{
  vector<float> offsets(3); // Spectrum shift offsets from the start of the merged spectrum (spectra 1,2,3 in pos 0,1,2)
  unsigned int intTol = (unsigned int)round(10 * peakTol);
  float tolOffset; // Keeps track of how far a mass is from the generating spectrum peak (within tolerance)
  unsigned int specIdx, peakIdx, tolIdx, specMassPos;
  Spectrum *specs[3]; // Pointers to spec1,spec2 and spec3 - to simplify the merging iterations
  specs[0] = &spec1;
  specs[1] = &spec2;
  specs[2] = &spec3;

  offsets[0] = min(shift12, shift13);
  if (offsets[0] >= 0)
    offsets[0] = 0;
  else
    offsets[0] = -offsets[0];
  mergedSpec.parentMass = spec1.parentMass + offsets[0];
  offsets[1] = offsets[0] + shift12;
  mergedSpec.parentMass = max(spec2.parentMass + offsets[1],
                              mergedSpec.parentMass);
  offsets[2] = offsets[0] + shift13;
  mergedSpec.parentMass = max(spec3.parentMass + offsets[2],
                              mergedSpec.parentMass);
  unsigned int newSize = (unsigned int)round(10 * mergedSpec.parentMass)
      + intTol + 1;
  if (newSize > gaux_peaksScores.size())
  {
    gaux_peaksScores.resize(newSize);
    gaux_peaksIdx.resize(newSize);
  }

  // Fill in the peaks masses and scores
  for (peakIdx = 0; peakIdx < gaux_peaksIdx.size(); peakIdx++)
  {
    gaux_peaksScores[peakIdx] = 0;
    gaux_peaksIdx[peakIdx].resize(3);
    for (specIdx = 0; specIdx < 3; specIdx++)
      gaux_peaksIdx[peakIdx][specIdx] = -1;
  }
  for (specIdx = 0; specIdx < 3; specIdx++)
    for (peakIdx = 0; peakIdx < specs[specIdx]->size(); peakIdx++)
    {
      specMassPos = (unsigned int)round(10 * (offsets[specIdx]
          + (*specs[specIdx])[peakIdx][0]));
      if (offsets[specIdx] + (*specs[specIdx])[peakIdx][0] >= peakTol)
        tolOffset = -peakTol;
      else
        tolOffset = -(offsets[specIdx] + (*specs[specIdx])[peakIdx][0]);
      for (tolIdx = max(0, (int)specMassPos - (int)intTol); tolIdx
          <= specMassPos + intTol; tolIdx++, tolOffset += .1)
      {
        if (gaux_peaksIdx[tolIdx][specIdx] >= 0)
        {
          if ((*specs[specIdx])[gaux_peaksIdx[tolIdx][specIdx]][1]
              >= (*specs[specIdx])[peakIdx][1])
            continue;
          else
            gaux_peaksScores[tolIdx]
                -= ((*specs[specIdx])[gaux_peaksIdx[tolIdx][specIdx]][1]
                    - fabs((*specs[specIdx])[gaux_peaksIdx[tolIdx][specIdx]][0]
                        + offsets[specIdx] - ((float)tolIdx) / 10.0));
        }
        gaux_peaksIdx[tolIdx][specIdx] = peakIdx;
        gaux_peaksScores[tolIdx] += (*specs[specIdx])[peakIdx][1]
            - fabs(tolOffset); // Rightmost term is a small linear penalty for tolerance offsets
      }
    }

  // Find all non-empty positions and define the merged spectrum
  unsigned int count = 0;
  for (peakIdx = 0; peakIdx < gaux_peaksScores.size(); peakIdx++)
    count += (gaux_peaksScores[peakIdx] > 0 and (((gaux_peaksIdx[peakIdx][0]
        >= 0) + (gaux_peaksIdx[peakIdx][1] >= 0) + (gaux_peaksIdx[peakIdx][2]
        >= 0)) >= minPeakCount)) ? 1 : 0;
  mergedSpec.resize(count);
  peaksIndices.resize(count);
  unsigned int mergedIdx = 0;
  for (peakIdx = 0; peakIdx < gaux_peaksScores.size(); peakIdx++)
    if (gaux_peaksScores[peakIdx] > 0)
    {
      if (((gaux_peaksIdx[peakIdx][0] >= 0) + (gaux_peaksIdx[peakIdx][1] >= 0)
          + (gaux_peaksIdx[peakIdx][2] >= 0)) >= minPeakCount)
      {
        mergedSpec[mergedIdx].set(peakIdx / 10.0, gaux_peaksScores[peakIdx]);
        peaksIndices[mergedIdx].resize(3);
        for (int i = 0; i < 3; i++)
          peaksIndices[mergedIdx][i] = gaux_peaksIdx[peakIdx][i];
        mergedIdx++;
      }
    }

  return count;
}

void SpecSet::SaveScanNums(const char *filename)
{
  vector<unsigned int> scanNums(specs.size());
  for (unsigned int i = 0; i < specs.size(); i++)
    scanNums[i] = specs[i].scan;
  Save_binArray(filename, scanNums);
}

//
//  SaveSpecSet_pklbin - saves a set of spectra in binary format. File format
//    1 int - number of spectra in the file
//    numSpecs shorts - number of peaks per spectrum in the file
//    arrays of numPeaks+1 pairs of floats - [parentMass charge] + array of [mass intensity]
//
short SaveSpecSet_pklbin(const char *filename, vector<Spectrum *> specs)
{
  FILE *fp;
  unsigned int numSpecs = specs.size();
  unsigned int i, p;

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numSpecs, sizeof(int), 1, fp); // Number of spectra in the file

  unsigned short *numPeaks = (unsigned short *)malloc(sizeof(short) * numSpecs);
  unsigned short maxNumPeaks = 0;
  for (i = 0; i < numSpecs; i++)
  {
    numPeaks[i] = specs[i]->size();
    maxNumPeaks = max(maxNumPeaks, numPeaks[i]);
  }
  fwrite(numPeaks, sizeof(short), numSpecs, fp); // Number of peaks per spectrum in the file

  float *peaksBuffer = (float *)malloc(2 * (maxNumPeaks + 1) * sizeof(float));
  unsigned int pbIdx;
  for (i = 0; i < numSpecs; i++)
  {
    peaksBuffer[0] = specs[i]->parentMass;
    peaksBuffer[1] = (float)specs[i]->parentCharge;
    for (pbIdx = 2, p = 0; p < numPeaks[i]; p++)
    {
      peaksBuffer[pbIdx++] = (*specs[i])[p][0];
      peaksBuffer[pbIdx++] = (*specs[i])[p][1];
    }
    fwrite(peaksBuffer, sizeof(float), 2 * (numPeaks[i] + 1), fp); // [parentMass charge] followed by [masses intensities]
  }

  free(peaksBuffer);
  free(numPeaks);
  fclose(fp);
  return 1;
}

