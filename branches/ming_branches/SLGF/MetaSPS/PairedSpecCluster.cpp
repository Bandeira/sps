/*
 * PairedSpecCluster.cpp
 *
 *  Created on: Nov 30, 2011
 *      Author: aguthals
 */

#include "PairedSpecCluster.h"
#include "Logger.h"

namespace specnets
{
  void PairedSpecCluster::initialize()
  {
    numSpecs = -1;
    int numCID = -1;
    int numHCD = -1;
    int numETD = -1;
    SpecSet* templateSpecs;

    if (spectraSetCID != 0) {
      numSpecs = spectraSetCID->size();
      numCID = spectraSetCID->size();
      templateSpecs = spectraSetCID;
    }
    if (spectraSetHCD != 0) {
      numSpecs = spectraSetHCD->size();
      numHCD = spectraSetHCD->size();
      templateSpecs = spectraSetHCD;
    }
    if (spectraSetETD != 0) {
      numSpecs = spectraSetETD->size();
      numETD = spectraSetETD->size();
      templateSpecs = spectraSetETD;
    }
    if (numSpecs < 0) {
      ERROR_MSG("No spectra have been loaded, initialization of PairedSpecCluster failed!");
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
      mergedSet = new SpecSet(1);
    }
    else {
      mergedSet->resize(1);
    }

    if (mergedLabels == 0) {
      mergedLabels = new vector<vector<short> > (1);
    }
    else {
      mergedLabels->resize(1);
    }

    if (spectraSetCID != 0) {
      if (usedPeaksCID == 0) {
        usedPeaksCID = new vector<vector<short> > (numSpecs);
      }
      initializeUsedPeaks(usedPeaksCID, spectraSetCID);
    }
    else {
      if (usedPeaksCID != 0) {
        delete usedPeaksCID;
      }
    }

    if (spectraSetHCD != 0) {
      if (usedPeaksHCD == 0) {
        usedPeaksHCD = new vector<vector<short> > (numSpecs);
      }
      initializeUsedPeaks(usedPeaksHCD, spectraSetHCD);
    }
    else {
      if (usedPeaksHCD != 0) {
        delete usedPeaksHCD;
      }
    }

    if (spectraSetETD != 0) {
      if (usedPeaksETD == 0) {
        usedPeaksETD = new vector<vector<short> > (numSpecs);
      }
      initializeUsedPeaks(usedPeaksETD, spectraSetETD);
    }
    else {
      if (usedPeaksETD != 0) {
        delete usedPeaksETD;
      }
    }

    float totalScore = 0, numScores = 0, averageScore = 0, standardDev = 0,
        endPtScore = 0;
    list<float> scores;
    for (int i = 0; i < numSpecs; i++) {
      float massB0 = 0;
      float massBk = (*templateSpecs)[i].parentMass - AAJumps::massMH;
      float massY0 = AAJumps::massH2O;
      float massYk = (*templateSpecs)[i].parentMass - AAJumps::massHion;

      if (spectraSetCID != 0) {
        for (int j = 0; j < (*spectraSetCID)[i].size(); j++) {
          if (MZRange::EqualWithinRange(massB0,
                                        (*spectraSetCID)[i][j][0],
                                        (*spectraSetCID)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massBk,
                                           (*spectraSetCID)[i][j][0],
                                           (*spectraSetCID)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massY0,
                                           (*spectraSetCID)[i][j][0],
                                           (*spectraSetCID)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massYk,
                                           (*spectraSetCID)[i][j][0],
                                           (*spectraSetCID)[i].getTolerance(j))) {
            (*usedPeaksCID)[i][j] = 3;
            continue;
          }
          totalScore += (*spectraSetCID)[i][j][1];
          scores.push_back((*spectraSetCID)[i][j][1]);
          numScores += 1.0;
        }
      }

      if (spectraSetHCD != 0) {
        for (int j = 0; j < (*spectraSetHCD)[i].size(); j++) {
          if (MZRange::EqualWithinRange(massB0,
                                        (*spectraSetHCD)[i][j][0],
                                        (*spectraSetHCD)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massBk,
                                           (*spectraSetHCD)[i][j][0],
                                           (*spectraSetHCD)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massY0,
                                           (*spectraSetHCD)[i][j][0],
                                           (*spectraSetHCD)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massYk,
                                           (*spectraSetHCD)[i][j][0],
                                           (*spectraSetHCD)[i].getTolerance(j))) {
            (*usedPeaksHCD)[i][j] = 3;
            continue;
          }
          totalScore += (*spectraSetHCD)[i][j][1];
          scores.push_back((*spectraSetHCD)[i][j][1]);
          numScores += 1.0;
        }
      }

      float massZk = (*templateSpecs)[i].parentMass - AAJumps::massMH
          - AAJumps::massNH;

      if (spectraSetETD != 0) {
        for (int j = 0; j < (*spectraSetETD)[i].size(); j++) {
          if (MZRange::EqualWithinRange(massB0,
                                        (*spectraSetETD)[i][j][0],
                                        (*spectraSetETD)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massBk,
                                           (*spectraSetETD)[i][j][0],
                                           (*spectraSetETD)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massY0,
                                           (*spectraSetETD)[i][j][0],
                                           (*spectraSetETD)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massYk,
                                           (*spectraSetETD)[i][j][0],
                                           (*spectraSetETD)[i].getTolerance(j))
              || MZRange::EqualWithinRange(massZk,
                                           (*spectraSetETD)[i][j][0],
                                           (*spectraSetETD)[i].getTolerance(j))) {
            (*usedPeaksETD)[i][j] = 3;
            continue;
          }
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

    endPtScore = averageScore + (2.0 * standardDev);

    int insertIdx;
    (*mergedSet)[0].copyNP((*templateSpecs)[0]);
    (*mergedSet)[0].resize(0);
    float massB0 = 0;
    float massBk = (*mergedSet)[0].parentMass - AAJumps::massMH;
    float massY0 = AAJumps::massH2O;
    float massYk = (*mergedSet)[0].parentMass - AAJumps::massHion;
    insertIdx = (*mergedSet)[0].insertPeak(massB0, endPtScore, 0);
    insertIdx = (*mergedSet)[0].insertPeak(massBk,
                                           endPtScore,
                                           (*mergedSet)[0].parentMassTol);
    insertIdx = (*mergedSet)[0].insertPeak(massY0, endPtScore, 0);
    insertIdx = (*mergedSet)[0].insertPeak(massYk,
                                           endPtScore,
                                           (*mergedSet)[0].parentMassTol);
    (*mergedLabels)[0].resize(0);
    (*mergedLabels)[0].resize(4, (short) 3);

  }

  void PairedSpecCluster::moveSamePRMs(Spectrum::FragType fragCIDHCD)
  {
    for (int i = 0; i < numSpecs; i++) {
      moveSamePRMsPair(0, i, i, fragCIDHCD);
    }
  }

  void PairedSpecCluster::moveSameSRMs(bool checkPRM,
                                       Spectrum::FragType fragCIDHCD)
  {
    list<int> idxCheck;
    for (int i = 0; i < numSpecs; i++) {
      idxCheck.clear();
      idxCheck.push_back(i);
      moveSameSRMsPair(0,
                       i,
                       i,
                       checkPRM,
                       fragCIDHCD,
                       &idxCheck,
                       &idxCheck,
                       &idxCheck);
    }
  }

  void PairedSpecCluster::movePRMsSRMs(Spectrum::FragType fragCIDHCD)
  {
    for (int i = 0; i < numSpecs; i++) {
      movePRMsSRMsPair(0, i, i, fragCIDHCD);
    }
  }

  void PairedSpecCluster::moveLeftovers(Spectrum::FragType fragCIDHCDETD)
  {

    for (int i = 0; i < numSpecs; i++) {
      moveLeftoversPair(0, i, fragCIDHCDETD);
    }
  }

  void PairedSpecCluster::moveSamePRMsAllPairs(Spectrum::FragType fragCIDHCD)
  {
    for (int i = 0; i < numSpecs; i++) {
      for (int j = i + 1; j < numSpecs; j++) {
        moveSamePRMsPair(0, i, j, fragCIDHCD);
      }
    }
  }

  void PairedSpecCluster::moveSameSRMsAllPairs(bool checkPRM,
                                               Spectrum::FragType fragCIDHCD)
  {
    list<int> idxCheck;
    for (int i = 0; i < numSpecs; i++) {
      idxCheck.push_back(i);
    }
    for (int i = 0; i < numSpecs; i++) {
      for (int j = i + 1; j < numSpecs; j++) {
        moveSameSRMsPair(0,
                         i,
                         j,
                         checkPRM,
                         fragCIDHCD,
                         &idxCheck,
                         &idxCheck,
                         &idxCheck);
      }
    }
  }

  void PairedSpecCluster::movePRMsSRMsAllPairs(Spectrum::FragType fragCIDHCD)
  {
    for (int i = 0; i < numSpecs; i++) {
      for (int j = i + 1; j < numSpecs; j++) {
        movePRMsSRMsPair(0, i, j, fragCIDHCD);
      }
    }
  }
}
