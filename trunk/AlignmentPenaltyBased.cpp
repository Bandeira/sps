#include "AlignmentPenaltyBased.h"
#include "Logger.h"

#include <algorithm>
//#include <cmath>
#include <iostream>
#include <string>
#include <limits.h>
#include <stdio.h>


#define DEBUG_RANGE 0
#define DEBUG_RANGE2 0
#define DEBUG_SPECS 0
#define DEBUG_ALIGN 0
#define DEBUG_ALIGN2 0
#define DEBUG_ALIGN3 0
#define DEBUG_GAP 0
#define DEBUG_GAP2 0
#define DEBUG_GAP3 0
#define DEBUG_GAP4 0
#define DEBUG_GAP5 0


namespace specnets
{
  //--------------------------------------------------------------------------
  string massToString(float mass)
  {
	char strMod[100];
       sprintf(strMod, "%.3f", mass);
	return string(strMod);
  }


  //--------------------------------------------------------------------------
  void fillVectors(size_t           startIndex,
                   size_t           index,
                   size_t           aaCounter, 
                   string           strAA,
                   vector<float>  & scores,
                   vector<int>    & aaNum,
                   vector<int>    & backPtr,
                   vector<string> & matchType,
                   PenaltyMatrix *  penaltyMatrix,
                   float            avgPeakIntensity,
                   float            mult,
                   float            maxMod,
                   float            minMod,
                   bool             lastAA,
                   bool             addPrevious)
  {
    // No penalty for hitting DB mass exactly
    if (index >= 0 && 
        index < scores.size() &&
        (lastAA || index != scores.size() - 1)) {
      float newPenalty = 0;
      if (addPrevious) {
        newPenalty += scores[startIndex];
      }
      if (newPenalty > scores[index]) {
        scores[index] = newPenalty;
        if (DEBUG_GAP4) DEBUG_MSG(index << "  " << scores[index]);
        aaNum[index] = aaCounter;
        matchType[index] = strAA;
        backPtr[index] = startIndex;
        if (DEBUG_GAP2) DEBUG_MSG(matchType[index] << " at " << index);
      }
    }

    if (DEBUG_GAP2) DEBUG_VAR(aaCounter);
    if (DEBUG_GAP2) DEBUG_VAR(lastAA);
    if (DEBUG_GAP2) DEBUG_VAR(scores.size());

    std::map<float, float> penalties;
    penaltyMatrix->getPenalties(strAA, penalties, avgPeakIntensity);
    std::map<float, float>::iterator itrm = penalties.begin();
    std::map<float, float>::iterator itrm_end = penalties.end();
    // Loop over all the possible mods (and their penalties)
    for (;itrm != itrm_end; itrm++) {

      float massDiff = itrm->first;
      float penalty = itrm->second;
      if (DEBUG_GAP2) DEBUG_MSG(massDiff << ", " << penalty);

      if (penalty == 0.0) {
        continue;
      }

      size_t indexPenalty = (size_t)(index + massDiff);
      //if (DEBUG_GAP2) DEBUG_VAR(indexPenalty);

      // Penalties are negative so a "higher" penalty is better
      float newPenalty = mult * penalty;
      if (DEBUG_GAP2) DEBUG_MSG(mult << "  " << penalty);
      if (addPrevious) {
        if (DEBUG_GAP2) DEBUG_MSG(startIndex << "  " << scores[startIndex]);
        newPenalty += scores[startIndex];
      }
      if (indexPenalty >= 0 && 
          indexPenalty < scores.size() && 
          newPenalty > scores[indexPenalty] &&
          (lastAA || indexPenalty != scores.size() - 1)) {
        scores[indexPenalty] = newPenalty;
        if (DEBUG_GAP4) DEBUG_MSG(indexPenalty << "  " << scores[indexPenalty]);
        aaNum[indexPenalty] = aaCounter;
        backPtr[indexPenalty] = startIndex;
        matchType[indexPenalty] = "(" + strAA + "," + massToString(massDiff) + ")";
        if (DEBUG_GAP4) DEBUG_MSG(matchType[indexPenalty] << " at " << indexPenalty 
                                 << " diff[" << massDiff << "] pen[" << newPenalty << "]");
      }
    }

    // Fill in all possible unknown modifications
    int minRealDelta = (int)(penaltyMatrix->getMass(strAA) - 57.0);
    if (DEBUG_GAP2) DEBUG_VAR(index);
    if (DEBUG_GAP2) DEBUG_VAR(minRealDelta);
    if (DEBUG_GAP2) DEBUG_VAR(minMod);
    int startUnk = max((int)(index - minRealDelta), (int)(index + minMod)); // min mod is negative
    startUnk = max(0, startUnk);  // Don't go beyond beginning of array
    int endUnk = min((int)(index + maxMod), (int)scores.size());
    float unkPenalty = penaltyMatrix->getUnknownPenalty(avgPeakIntensity);
    if (DEBUG_GAP2) DEBUG_VAR(strAA);
    if (DEBUG_GAP2) DEBUG_VAR(lastAA);
    if (DEBUG_GAP2) DEBUG_VAR(startUnk);
    if (DEBUG_GAP2) DEBUG_VAR(endUnk);
    if (DEBUG_GAP2) DEBUG_VAR(unkPenalty);
    if (DEBUG_GAP2) DEBUG_VAR(addPrevious);
    for (int iUnk = startUnk; iUnk < endUnk; iUnk++) {
      float newPenalty = unkPenalty;
      if (addPrevious) {
        newPenalty += scores[startIndex];
      }
      if (newPenalty > scores[iUnk] &&
          (lastAA || iUnk != scores.size() - 1)) {
        scores[iUnk] = newPenalty;
        if (DEBUG_GAP4) DEBUG_MSG(iUnk << "  " << scores[iUnk]);
        aaNum[iUnk] = aaCounter;
        backPtr[iUnk] = startIndex;
        float massDiff = (float)iUnk - (float)index;
        matchType[iUnk] = "(" + strAA + "," + massToString(massDiff) + ")";
        if (DEBUG_GAP3) DEBUG_MSG(matchType[iUnk] << " at " << iUnk << " diff[" << massDiff << "] pen[" << newPenalty << "]");
      }
    }

    return;
  }

  //--------------------------------------------------------------------------
  float computeGapPenalty(Spectrum & 	spec,
                          Spectrum & 	dbSpec,
                          char * 		dbSeq,
                          int 		dbSpecIdx,
                          int 		specIdx,
                          int 		predDbSpecIdx,
                          int 		predSpecIdx,
                          PenaltyMatrix * modPenaltyMatrix,
                          PenaltyMatrix * blossumPenaltyMatrix,
                          float           avgPeakIntensity,
                          float	       alpha,
                          float           beta,
                          float           maxMod,
                          float           minMod,
                          string & matchString)
  {
    if (DEBUG_GAP) DEBUG_VAR(predDbSpecIdx);
    if (DEBUG_GAP) DEBUG_VAR(dbSpecIdx);
    if (DEBUG_GAP) DEBUG_VAR(predSpecIdx);
    if (DEBUG_GAP) DEBUG_VAR(specIdx);

    float bestScore = 0.0;
    float massLength = (int)(spec[specIdx][0] - spec[predSpecIdx][0] + 0.5) + 1;
    if (DEBUG_GAP) DEBUG_VAR(massLength);
    vector<float> scores((size_t)massLength);
    vector<int>   aaNum((size_t)massLength);
    vector<int>   backPtr((size_t)massLength);
    vector<string> matchType((size_t)massLength);

    // Initialize the vectors
    for (int i = 0; i < scores.size(); i++) {
      scores[i] = -(float)INT_MAX;
      aaNum[i] = 0;
      backPtr[i] = -1;
    }

    int dbGapSize = dbSpecIdx - predDbSpecIdx;
    if (DEBUG_GAP) DEBUG_VAR(dbGapSize);
    char aa = dbSeq[predDbSpecIdx];
    if (DEBUG_GAP) DEBUG_VAR(aa);
    string strAA("X");  // Create a dummy string
    strAA[0] = aa;      // Set the string to the AA char

    float aaMass = modPenaltyMatrix->getMass(strAA);
    if (DEBUG_GAP) DEBUG_VAR(aaMass);

    // Sanity check for 0 masses
    if (aaMass == 0.0) {
      WARN_MSG("0 Mass for AA [" <<  strAA << "]");

      return -(float)INT_MAX;
    }

    fillVectors(-1, (size_t)aaMass, 1, strAA, scores, aaNum, backPtr, matchType, 
                modPenaltyMatrix, avgPeakIntensity, alpha, maxMod, minMod, dbGapSize == 1, false);

    fillVectors(-1, (size_t)aaMass, 1, strAA, scores, aaNum, backPtr, matchType, 
                blossumPenaltyMatrix, avgPeakIntensity, beta, maxMod, minMod, dbGapSize == 1, false);

    if (DEBUG_GAP) {
      for (int j = 0; j < scores.size(); j++) {
	 cout << scores[j] << "\t";
      }
      cout << endl;
      for (int j = 0; j < scores.size(); j++) {
	 cout << aaNum[j] << "\t";
      }
      cout << endl;
      for (int j = 0; j < scores.size(); j++) {
	 cout << backPtr[j] << "\t";
      }
      cout << endl;
      for (int j = 0; j < scores.size(); j++) {
	 cout << matchType[j] << "\t";
      }
      cout << endl;
      cout << endl;
    }

    int aaCounter = 1;
    for (int i = predDbSpecIdx + 1; i < dbSpecIdx; i++) {
      char aa = dbSeq[i];
      if (DEBUG_GAP) DEBUG_VAR(aa);
      string strAA("X");  // Create a dummy string
      strAA[0] = aa;      // Set the string to the AA char
      float aaMass = modPenaltyMatrix->getMass(strAA);
      if (DEBUG_GAP) DEBUG_VAR(aaMass);

      // Sanity check for 0 masses
      if (aaMass == 0.0) {
        return -(float)INT_MAX;
      }

      aaCounter++;
      if (DEBUG_GAP) DEBUG_VAR(aaCounter);

      for (int startIndex = 1; startIndex < scores.size(); startIndex++) {
      
      	 // We can only start at places that were set by the previous AA
      	 // This way we don't skip any AA in the database sequence
      	 if (aaNum[startIndex] != aaCounter - 1) {
          continue;
      	 }
        size_t index = startIndex + (size_t)aaMass;

        fillVectors(startIndex, index, aaCounter, strAA, scores, aaNum, backPtr, matchType, 
                    modPenaltyMatrix, avgPeakIntensity, alpha, maxMod, minMod, dbGapSize == aaCounter, true);
        fillVectors(startIndex, index, aaCounter, strAA, scores, aaNum, backPtr, matchType, 
                    blossumPenaltyMatrix, avgPeakIntensity, beta, maxMod, minMod, dbGapSize == aaCounter, true);

      } // for (int startIndex = 1; startIndex < scores.size(); startIndex++) {

    } // for (int i = predDbSpecIdx + 1; i < dbSpecIdx; i++) {

    if (DEBUG_GAP) {
      for (int j = 0; j < scores.size(); j++) {
        cout << scores[j] << "\t";
      }
      cout << endl;
      for (int j = 0; j < scores.size(); j++) {
        cout << aaNum[j] << "\t";
      }
      cout << endl;
      for (int j = 0; j < scores.size(); j++) {
        cout << backPtr[j] << "\t";
      }
      cout << endl;
      for (int j = 0; j < scores.size(); j++) {
	 cout << matchType[j] << "\t";
      }
      cout << endl;
      cout << endl;
    }

    size_t backIndex = (size_t)(massLength - 1);
    while (backIndex != -1) {
      if (DEBUG_GAP) DEBUG_VAR(backIndex);
      matchString = matchType[backIndex] + matchString; 
      if (DEBUG_GAP) DEBUG_VAR(matchString);
      backIndex = backPtr[backIndex];
    }
    //matchString = "[" + matchString +"]";

    if (DEBUG_GAP) DEBUG_VAR(scores[(size_t)(massLength - 1)]);
    return scores[(size_t)(massLength - 1)];
  }

  
  //--------------------------------------------------------------------------
  void alignmentPenaltyBased(Spectrum & spec,
                             Spectrum & dbSpec,
                             char * dbSeq,
                             int dbIndex,
                             int matchOrientation,
                             set<float> & startRange,
                             PenaltyMatrix * modPenaltyMatrix,
                             PenaltyMatrix * blossumPenaltyMatrix,
                             float penaltyAlpha,
                             float penaltyBeta,
                             short maxNumMods,
                             int minMatchedPeaks,
                             int maxGapSize,
                             float pmTolerance,
                             float tolerance,
                             float minSpecDist,
                             float maxDbSpecMod,
                             float minDbSpecMod,
                             bool enforceEndpeaks)
  {
    if (spec.size() == 0 or spec.size() < minMatchedPeaks) {
      WARN_MSG("Spectrum size [" <<  spec.size() << 
               "] is less than minimum matched peaks [" << minMatchedPeaks << "]");
      return;
    }

#if DEBUG_ALIGN
    DEBUG_VAR(enforceEndpeaks);
    DEBUG_VAR(maxGapSize);
    DEBUG_VAR(penaltyAlpha);
    DEBUG_VAR(penaltyBeta);
    DEBUG_VAR(startRange.size());
    DEBUG_VAR(pmTolerance);
    DEBUG_VAR(tolerance);
    DEBUG_VAR(maxDbSpecMod);
    DEBUG_VAR(minDbSpecMod);
#endif

    // Compute average peak intensity
    float avgPeakIntensity = 0;
    for (int i = 0; i < spec.size(); i++) {
      avgPeakIntensity += spec[i][1];
    }
    avgPeakIntensity /= spec.size();
    if (DEBUG_ALIGN) DEBUG_VAR(avgPeakIntensity);

    vector<vector< float > > matchMatrix; 
    vector<vector< pair<int,int> > > matchPtr; 
    vector<vector< string > > matchType;

    // Initialize starting flag array
    vector<vector<char> > startFlags(spec.size());
    for (int i = 0; i < spec.size(); i++) {
      vector<char> newArray(dbSpec.size());
      startFlags[i] = newArray;
      for (int j = 0; j < dbSpec.size(); j++) {
        // If the startRange set is empty then allow any start position
        if (startRange.size() == 0) { 
          startFlags[i][j] = 1;
        } else {
          startFlags[i][j] = 0;
        }
      }
    }

    // Find all the valid starting points in the matrix    
    set<float>::iterator itr = startRange.begin();
    set<float>::iterator itrEnd = startRange.end();
    for (; itr != itrEnd; itr++) {
      int minIdx1, maxIdx1, minIdx2, maxIdx2;
      float minStartMass = *itr + minDbSpecMod;
      float maxStartMass = *itr + maxDbSpecMod;
      if (DEBUG_RANGE) DEBUG_VAR(minStartMass);
      if (DEBUG_RANGE) DEBUG_VAR(maxStartMass);
      float lastMass = 0;
      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
        float massDiff = spec[specIdx][0] - lastMass;
        if (DEBUG_RANGE) DEBUG_VAR(massDiff);
        minStartMass += massDiff;
        maxStartMass += massDiff;
        if (DEBUG_RANGE) DEBUG_VAR(minStartMass);
        if (DEBUG_RANGE) DEBUG_VAR(maxStartMass);
        lastMass = spec[specIdx][0];
        // Find peaks close to desired bounds
        float T = 200.0;

        list<int> matches1;
        list<int> matches2;
        dbSpec.setPeakTolerance(0);
        dbSpec.findPeaks(minStartMass, T, &matches1);
        dbSpec.findPeaks(maxStartMass, T, &matches2);
        if (DEBUG_RANGE) DEBUG_VAR(matches1.size());
        if (DEBUG_RANGE) DEBUG_VAR(matches1.front());
        if (DEBUG_RANGE) DEBUG_VAR(matches1.back());
        if (DEBUG_RANGE) DEBUG_VAR(matches2.size());
        if (DEBUG_RANGE) DEBUG_VAR(matches2.front());
        if (DEBUG_RANGE) DEBUG_VAR(matches2.back());
        minIdx1 = (matches1.size() == 0) ? -1 : matches1.front();
        maxIdx1 = (matches1.size() == 0) ? -1 : matches1.back();
        minIdx2 = (matches2.size() == 0) ? -1 : matches2.front();
        maxIdx2 = (matches2.size() == 0) ? -1 : matches2.back();
/*
        dbSpec.locateWithinTolerance(0,
                                     minStartMass,
                                     T,
                                     minIdx1,
                                     maxIdx1);
        dbSpec.locateWithinTolerance(0,
                                     maxStartMass,
                                     T,
                                     minIdx2,
                                     maxIdx2);
*/
        // If peak is not found, index could be -1.. so make sure at least 0
        minIdx1 = max<int>(minIdx1,0);
        minIdx2 = max<int>(minIdx2,0);
        // Make sure min peak is INSIDE the desired range
        while (minIdx1 < dbSpec.size() && 
               dbSpec[minIdx1][0] < minStartMass - tolerance - AAJumps::massH2O) minIdx1++;

        // Make sure max peak is OUTSIDE (above) the desired range
        while (minIdx2 < dbSpec.size() && 
               dbSpec[minIdx2][0] < maxStartMass + tolerance + AAJumps::massH2O) minIdx2++;

        // Mark all peaks in range as valid starting points in matrix        
        if (DEBUG_RANGE) DEBUG_VAR(minIdx1);
        if (DEBUG_RANGE) DEBUG_VAR(minIdx2);
        for (int i = minIdx1; i < minIdx2 && i < dbSpec.size(); i++) {
          startFlags[specIdx][i] = 1;
        }
      } // for (int specIdx = 0; specIdx < spec.size(); specIdx++)
    } // for (; itr != itrEnd; itr++) {

    matchMatrix.resize(dbSpec.size());
    matchPtr.resize(dbSpec.size());
    matchType.resize(dbSpec.size());

#if DEBUG_SPECS
    DEBUG_MSG(0 << ", " << dbSpec[0][0] << ", " << dbSpec[0][0]);
    for (int dbSpecIdx = 1; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      DEBUG_MSG(dbSpecIdx << ", " << dbSpec[dbSpecIdx][0] << ", " << dbSpec[dbSpecIdx][0] - dbSpec[dbSpecIdx-1][0]);
    }
    DEBUG_MSG(0 << ", " << spec[0][0] << ", " << spec[0][1]);
    for (int specIdx = 1; specIdx < spec.size(); specIdx++) {
      DEBUG_MSG(specIdx << ", " << spec[specIdx][0] - spec[specIdx-1][0] << ", " << spec[specIdx][1]);
    }
#endif

    // Initialize the matrices    
    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      matchMatrix[dbSpecIdx].resize(spec.size());
      matchPtr[dbSpecIdx].resize(spec.size());
      matchType[dbSpecIdx].resize(spec.size());
      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
        matchMatrix[dbSpecIdx][specIdx] = -(float)INT_MAX;
        matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(-1,-1);
        matchType[dbSpecIdx][specIdx] = "?";
      }
    }

    // Compute the match between the two spectra
    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {

      for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
    
        if (DEBUG_RANGE2) DEBUG_MSG(dbSpecIdx << ", " << specIdx << ", " << (int)startFlags[specIdx][dbSpecIdx]);

        if (!enforceEndpeaks or specIdx == 0) {
          // Matching can start on any mass, assume no predecessor 
          matchMatrix[dbSpecIdx][specIdx] = spec[specIdx][1];
        }

        if (DEBUG_RANGE2) DEBUG_MSG(matchMatrix[dbSpecIdx][specIdx]);

        for (int predDbSpecIdx = dbSpecIdx - 1;
             (predDbSpecIdx >= 0) && (dbSpecIdx - predDbSpecIdx < maxGapSize);
              predDbSpecIdx--) {

          for (int predSpecIdx = specIdx - 1; predSpecIdx >= 0; predSpecIdx--) {

            if (!startFlags[predSpecIdx][predDbSpecIdx]) {
              // if location is not a valid start, skip scoring
              if (DEBUG_RANGE2) DEBUG_MSG("Skipping " << predDbSpecIdx << ", " << predSpecIdx);
              continue;
            }
            if (matchMatrix[predDbSpecIdx][predSpecIdx] == -(float)INT_MAX) {
              // Must build on a real score (can make anything out of -infinity)
              if (DEBUG_RANGE2) DEBUG_MSG("Skipping " << predDbSpecIdx << ", " << predSpecIdx);
              continue;
            }
        
#if DEBUG_ALIGN
            DEBUG_TRACE;
            DEBUG_VAR(predDbSpecIdx);
            DEBUG_VAR(dbSpecIdx);
            DEBUG_VAR(predSpecIdx);
            DEBUG_VAR(specIdx);
#endif
            float dbGapMass = dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0];
            float dbAAMinMass = (dbSpecIdx - predDbSpecIdx) * 57.0;
            float specGapMass = spec[specIdx][0] - spec[predSpecIdx][0];
            if (DEBUG_ALIGN) DEBUG_VAR(dbAAMinMass);
            if (DEBUG_ALIGN) DEBUG_VAR(specGapMass);

            // Is this a larger gap (and it can be filled with the right number of AA's min 57 size)
            if (dbSpecIdx - predDbSpecIdx != 1 && dbAAMinMass < specGapMass && specGapMass < dbGapMass - minDbSpecMod) {
              if (DEBUG_ALIGN2) DEBUG_MSG("Gap match");
              string matchString;
              float penalty = computeGapPenalty(spec,
                                                  dbSpec,
                                                  dbSeq,
                                                  dbSpecIdx,
                                                  specIdx,
                                                  predDbSpecIdx,
                                                  predSpecIdx,
                                                  modPenaltyMatrix,
                                                  blossumPenaltyMatrix,
                                                  avgPeakIntensity,
                                                  penaltyAlpha,
                                                  penaltyBeta,
                                                  maxDbSpecMod,
                                                  minDbSpecMod,
                                                  matchString);
              if (DEBUG_ALIGN2) DEBUG_VAR(penalty);
              float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                spec[specIdx][1] + penalty;
              if (DEBUG_ALIGN2) DEBUG_VAR(score);
              if (score > matchMatrix[dbSpecIdx][specIdx]) {
                matchMatrix[dbSpecIdx][specIdx] = score;
                if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                matchType[dbSpecIdx][specIdx] = matchString;
                if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
              }
            } else {
            
              float deltaDbMass = dbSpec[dbSpecIdx][0] - dbSpec[predDbSpecIdx][0];
              if (deltaDbMass > 400.0) continue;
              float deltaSpecMass = spec[specIdx][0] - spec[predSpecIdx][0];
              if (deltaSpecMass > 400.0) continue;
              float deltaMasses = deltaSpecMass - deltaDbMass;
            
              char aa = dbSeq[predDbSpecIdx];
              string strAA("X");  // Create a dummy string
              strAA[0] = aa;      // Set the string to the AA char

#if DEBUG_ALIGN
              DEBUG_VAR(deltaDbMass);
              DEBUG_VAR(deltaSpecMass);
              DEBUG_VAR(deltaMasses);
              DEBUG_VAR(strAA);
              DEBUG_TRACE;
#endif

              if (deltaMasses < minDbSpecMod || deltaMasses > maxDbSpecMod) {
                continue;
              }

              if (abs(deltaMasses) < tolerance) {
              
                if (DEBUG_ALIGN2) DEBUG_MSG("Exact match");
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                      spec[specIdx][1];
                if (DEBUG_ALIGN2) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = strAA;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }

              } else if (modPenaltyMatrix->isKnown(strAA, deltaMasses)) {
              
                if (DEBUG_ALIGN2) DEBUG_MSG("Modification match");
                float penaltyMod = (*modPenaltyMatrix)(strAA, deltaMasses, avgPeakIntensity);
                if (DEBUG_ALIGN2) DEBUG_VAR(penaltyMod);
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                        spec[specIdx][1] +
                                        penaltyAlpha * penaltyMod;
                if (DEBUG_ALIGN2) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = "(" + strAA + "," + massToString(deltaMasses) + ")";
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }
#if 0
              } else if (blossumPenaltyMatrix->isKnown(strAA, deltaMasses)) {
              
                if (DEBUG_ALIGN2) DEBUG_MSG("Blossum match");
                float penaltyBlossum = (*blossumPenaltyMatrix)(strAA, deltaMasses, avgPeakIntensity);
                if (DEBUG_ALIGN2) DEBUG_VAR(penaltyBlossum);
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                        spec[specIdx][1] +
                                        penaltyBeta * penaltyBlossum;
                if (DEBUG_ALIGN2) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = "(" + strAA + "," + massToString(deltaMasses) + ")";
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }
#endif
              } else if (!modPenaltyMatrix->isKnown(strAA, deltaMasses) &&
                         modPenaltyMatrix->getMass(strAA) + deltaMasses >= 57.0) {

                if (DEBUG_ALIGN2) DEBUG_MSG("Unknown match");
                float penaltyMod = modPenaltyMatrix->getUnknownPenalty(avgPeakIntensity);
                if (DEBUG_ALIGN2) DEBUG_VAR(penaltyMod);
                float score = matchMatrix[predDbSpecIdx][predSpecIdx] +
                                  spec[specIdx][1] +
                                  penaltyMod;
                if (DEBUG_ALIGN2) DEBUG_VAR(score);
                if (score > matchMatrix[dbSpecIdx][specIdx]) {
                  matchMatrix[dbSpecIdx][specIdx] = score;
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[predDbSpecIdx][predSpecIdx]);
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchMatrix[dbSpecIdx][specIdx]);
                  matchPtr[dbSpecIdx][specIdx] = make_pair<int,int>(predDbSpecIdx,predSpecIdx);
                  matchType[dbSpecIdx][specIdx] = "(" + strAA + "," + massToString(deltaMasses) + ")";
                  if (DEBUG_ALIGN2) DEBUG_VAR(matchType[dbSpecIdx][specIdx]);
                }

              }  // else if (modPenaltyMatrix->isKnown(strAA, deltaMasses))

            } // if (dbSpecIdx - predDbSpecIdx != 1 && specIdx - predSpecIdx != 1) {

          } // predSpecIdx

        } // predDbSpecIdx

      } // specIdxS

    } // dbSpecIdx

#if DEBUG_ALIGN3
    for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
      bool outputRow = false;
      // Only output "interesting" rows
      for (int specIdx = 1; specIdx < spec.size(); specIdx++) {
        if (matchMatrix[dbSpecIdx][specIdx] != -(float)INT_MAX) {
          outputRow = true;
          break;
        }
      }
      if (1) {
        cout << dbSpecIdx << "\t";
        for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
          cout << matchMatrix[dbSpecIdx][specIdx] << "\t";
        }
        cout << endl;
      }
    }
#endif
          
    float bestScore = -(float)INT_MAX;
    pair<int,int> bestMatchPtr;

    if (enforceEndpeaks) {
      for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
        if (matchMatrix[dbSpecIdx][spec.size()-1] > bestScore) {
          bestScore = matchMatrix[dbSpecIdx][spec.size()-1];
          bestMatchPtr = make_pair<int,int>(dbSpecIdx,spec.size()-1);
#if DEBUG_ALIGN2
          DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);
#endif
        }
      }
    } else {
      for (int dbSpecIdx = 0; dbSpecIdx < dbSpec.size(); dbSpecIdx++) {
        for (int specIdx = 0; specIdx < spec.size(); specIdx++) {
          if (matchMatrix[dbSpecIdx][specIdx] > bestScore) {
            bestScore = matchMatrix[dbSpecIdx][specIdx];
            bestMatchPtr = make_pair<int,int>(dbSpecIdx,specIdx);
#if DEBUG_ALIGN2
            DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);
#endif
          }
        }
      }
    }
    
#if DEBUG_ALIGN2
    DEBUG_MSG(bestMatchPtr.first << ", " << bestMatchPtr.second);
#endif

    if (bestScore != -(float)INT_MAX) {

      // Create a PSM for this match      
      psmPtr p(new PeptideSpectrumMatch);
      p->m_score = bestScore;
      p->m_dbIndex = dbIndex,
      p->m_matchOrientation = matchOrientation;
      p->m_spectrum = &spec;

      int lastDbIndex = bestMatchPtr.first;
      
#if DEBUG_ALIGN2
      // Copy the peak matches to the PSM
      DEBUG_VAR(bestMatchPtr.first);
      DEBUG_VAR(bestMatchPtr.second);
      DEBUG_VAR(matchType[bestMatchPtr.first][bestMatchPtr.second]);
#endif
      p->m_annotation.clear();
      // Check for gap at end
      if (bestMatchPtr.second != spec.size() - 1) {
        float endMass = spec[spec.size() - 1][0] - spec[bestMatchPtr.second][0];
        p->m_annotation = "[" + massToString(endMass) + "]";
      }

      pair<int,int> prevMatchPtr = bestMatchPtr;
      while (bestMatchPtr.second != -1) {
        TwoValues<int> tmpIndices;
        tmpIndices[0] = bestMatchPtr.second;
        tmpIndices[1] = bestMatchPtr.first;
        p->m_matchedPeaks.push_back(tmpIndices);
        if (matchType[bestMatchPtr.first][bestMatchPtr.second] != "?") {
          p->m_annotation = matchType[bestMatchPtr.first][bestMatchPtr.second] + p->m_annotation;
	 }
        prevMatchPtr = bestMatchPtr;
        bestMatchPtr = matchPtr[bestMatchPtr.first][bestMatchPtr.second];

#if DEBUG_ALIGN2
        DEBUG_VAR(bestMatchPtr.first);
        DEBUG_VAR(bestMatchPtr.second);
        if (bestMatchPtr.first != -1 || bestMatchPtr.second != -1) {
          DEBUG_VAR(matchType[bestMatchPtr.first][bestMatchPtr.second]);
        }
#endif

      }
      int firstDbIndex = prevMatchPtr.first;
      p->m_startMass = dbSpec[firstDbIndex][0];

      // Check for gap at beginning
      if (prevMatchPtr.second != 0) {
        p->m_annotation = "[" + massToString(spec[prevMatchPtr.second][0]) + "]" + p->m_annotation;
      }

#if DEBUG_ALIGN2
      DEBUG_VAR(firstDbIndex);
      DEBUG_VAR(lastDbIndex);
#endif

      if (lastDbIndex - firstDbIndex > 2) {

#if DEBUG_ALIGN2
        DEBUG_VAR(p->m_annotation);
        DEBUG_VAR(p->m_matchedPeaks.size());
#endif

        reverse(p->m_matchedPeaks.begin(),p->m_matchedPeaks.end());

        spec.psmList.push_back(p);
      } 

    }  // if (bestScore != -(float)INT_MAX)

    return;
  } // alignmentPenaltyBased()

} // namespace specnets

