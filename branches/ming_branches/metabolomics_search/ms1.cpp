#include "ms1.h"
#include "alignment_scoring.h"
#include "batch.h"
#include <cmath>

namespace specnets
{

  void BinnedMatch(Spectrum &spec1,
                   Spectrum &spec2,
                   TwoValues<float> &matchedIC)
  {
    unsigned int idxPeaks1 = 0, idxPeaks2 = 0;

    matchedIC.set(0, 0);
    while (idxPeaks1 < spec1.size() and idxPeaks2 < spec2.size()) {
      if (fabs(spec1[idxPeaks1][0] - spec2[idxPeaks2][0]) <= 0.00001) {
        matchedIC[0] += spec1[idxPeaks1++][1];
        matchedIC[1] = spec2[idxPeaks2++][1];
      }
      else {
        if (spec1[idxPeaks1][0] > spec2[idxPeaks2][0])
          idxPeaks2++;
        else
          idxPeaks1++;
      }
    }
  }

  float BinnedMatchProduct(Spectrum &spec1, Spectrum &spec2)
  {
    unsigned int idxPeaks1 = 0, idxPeaks2 = 0;
    float product = 0;

    while (idxPeaks1 < spec1.size() and idxPeaks2 < spec2.size()) {
      if (fabs(spec1[idxPeaks1][0] - spec2[idxPeaks2][0]) <= 0.00001)
        product += spec1[idxPeaks1++][1] * spec2[idxPeaks2++][1];
      else {
        if (spec1[idxPeaks1][0] > spec2[idxPeaks2][0])
          idxPeaks2++;
        else
          idxPeaks1++;
      }
    }
    return product;
  }

  void AlignChromatography(SpecSet &run1, SpecSet &run2, float peakTol, vector<
      TwoValues<unsigned int> > &matchedScans, bool binSpectra)
  {

    vector<float> data[2]; // One vector is used to buffer the previous line in the 2D alignment matrix while
    //   the other contains the values being calculated for the current line.
    //   Which is which is determined by the variables idxCurr, idxPrev
    unsigned int idxCurr, idxPrev, // see data
        idxRun1, idxRun2; // spectrum indices in each run
    vector<float> tic1, tic2, // total ion current for each spectrum in each run. Used to compute % matched intensity
        tic1sq, tic2sq; // sum of squared peak intensities for each spectrum in each run. Used to compute normalized dot product
    vector<int> idxMatched1, idxMatched2; // Variables used in FindMatchPeaksAll / ScoreOverlap6
    vector<vector<bool> > predIsDiagonal, // Pos [i,j]==true if predecessor of [i,j] is [i-1,j-1]
        predIsLeft; // Pos [i,j]==true if predecessor of [i,j] is [i,j-1]
    // Pos predIsLeft[i,j]==false and predIsDiagonal[i,j]=false if predecessor of [i,j] is [i-1,j]

    // Initialize variables
    data[0].resize(run2.size());
    data[1].resize(run2.size());
    for (idxRun2 = 0; idxRun2 < run2.size(); idxRun2++) {
      data[0][idxRun2] = 0;
      data[1][idxRun2] = 0;
    }
    predIsDiagonal.resize(run1.size());
    predIsLeft.resize(run1.size());
    for (idxRun1 = 0; idxRun1 < run1.size(); idxRun1++) {
      predIsDiagonal[idxRun1].resize(run2.size());
      predIsLeft[idxRun1].resize(run2.size());
      for (idxRun2 = 0; idxRun2 < run2.size(); idxRun2++) {
        predIsDiagonal[idxRun1][idxRun2] = false;
        predIsLeft[idxRun1][idxRun2] = false;
      }
    }

    tic1.resize(run1.size());
    tic2.resize(run2.size());
    tic1sq.resize(run1.size());
    tic2sq.resize(run2.size());
    for (idxRun1 = 0; idxRun1 < run1.size(); idxRun1++) {
      if (binSpectra)
        run1[idxRun1].setResolution(1, true);
      tic1[idxRun1] = 0;
      tic1sq[idxRun1] = 0;
      for (unsigned int pivot = 0; pivot < run1[idxRun1].size(); pivot++) {
        tic1[idxRun1] += run1[idxRun1][pivot][1];
        tic1sq[idxRun1] += run1[idxRun1][pivot][1] * run1[idxRun1][pivot][1];
      }
    }
    for (idxRun2 = 0; idxRun2 < run2.size(); idxRun2++) {
      if (binSpectra)
        run2[idxRun2].setResolution(1, true);
      tic2[idxRun2] = 0;
      tic2sq[idxRun2] = 0;
      for (unsigned int pivot = 0; pivot < run2[idxRun2].size(); pivot++) {
        tic2[idxRun2] += run2[idxRun2][pivot][1];
        tic2sq[idxRun2] += run2[idxRun2][pivot][1] * run2[idxRun2][pivot][1];
      }
    }

    // Align the runs
    idxPrev = 0;
    idxCurr = 1;
    unsigned int bandWidth = 5000;
    float tmpScore; // temporary variable used to compute the match score between two spectra
    idxMatched1.reserve(20000);
    idxMatched2.reserve(20000);
    unsigned int limRun2left, limRun2rigth; // Limits for banded alignment
    for (idxRun1 = 0; idxRun1 < run1.size(); idxRun1++) {
      if (idxRun1 > bandWidth)
        limRun2left = idxRun1 - bandWidth;
      else
        limRun2left = 0;
      if (idxRun1 + bandWidth < run2.size())
        limRun2rigth = idxRun1 + bandWidth + 1;
      else
        limRun2rigth = run2.size();
      for (idxRun2 = limRun2left; idxRun2 < limRun2rigth; idxRun2++) {
        // Find match score
        if (tic1[idxRun1] > 0 and tic2[idxRun2] > 0) {
          /*  Maximum non-intersecting peak matches. Need to fix ScoreOverlap6 to guarantee non-intersection when minIPdist < peakTol
           * 				FindMatchPeaksAll(run1[idxRun1], run2[idxRun2], 0, peakTol, idx1all, idx2all);
           // NOTE: minimum inter-peak distance set to 0.1 Da
           ScoreOverlap6(run1[idxRun1], run2[idxRun2], 0, peakTol, idxMatched1, idxMatched2, 0.1);
           tmpScore=0; for(unsigned int pivot=0; pivot<idxMatched1.size(); pivot++) tmpScore+=run1[idxRun1][idxMatched1[pivot]][1];
           data[idxCurr][idxRun2] = tmpScore/tic1[idxRun1];
           tmpScore=0; for(unsigned int pivot=0; pivot<idxMatched2.size(); pivot++) tmpScore+=run2[idxRun2][idxMatched2[pivot]][1];
           data[idxCurr][idxRun2] += tmpScore/tic2[idxRun2]; */

          // Summed explained intensities for binned spectra
          /*				TwoValues<float> matchedIC;
           BinnedMatch(run1[idxRun1], run2[idxRun2], matchedIC);
           data[idxCurr][idxRun2] = matchedIC[0]/tic1[idxRun1] + matchedIC[1]/tic2[idxRun2];
           */
          // Normalized dot-product between binned spectra
          tmpScore = BinnedMatchProduct(run1[idxRun1], run2[idxRun2]);
          data[idxCurr][idxRun2] = tmpScore / sqrt(tic1sq[idxRun1]
              * tic2sq[idxRun2]);
        }
        else
          data[idxCurr][idxRun2] = 0;

        // Find best predecessor
        if (idxRun1 == 0) {
          if (idxRun2 == 0 or data[idxCurr][idxRun2] > data[idxCurr][idxRun2
              - 1])
            predIsDiagonal[idxRun1][idxRun2] = true;
          else {
            predIsLeft[idxRun1][idxRun2] = true;
            data[idxCurr][idxRun2] == data[idxCurr][idxRun2 - 1];
          }
          continue;
        }
        // Spectrum i is matched to spectrum j iff predIsDiagonal[i][j]==true
        if ((data[idxCurr][idxRun2] + data[idxPrev][idxRun2 - 1] + 0.00001)
            >= data[idxCurr][idxRun2 - 1] and // +0.00001 is a puny bias towards the diagonal to avoid zig-zagging around empty scans
            (data[idxCurr][idxRun2] + data[idxPrev][idxRun2 - 1] + 0.00001)
                >= data[idxPrev][idxRun2]) {
          predIsDiagonal[idxRun1][idxRun2] = true;
          data[idxCurr][idxRun2] += data[idxPrev][idxRun2 - 1];
          continue;
        }
        if (data[idxCurr][idxRun2 - 1] > data[idxPrev][idxRun2]) {
          predIsLeft[idxRun1][idxRun2] = true;
          data[idxCurr][idxRun2] = data[idxCurr][idxRun2 - 1];
        }
        else
          data[idxCurr][idxRun2] = data[idxPrev][idxRun2];
      }
      if ((idxRun1 % 100) == 0)
        cerr << "Done with scan " << run1[idxRun1].scan << " (run1)\n";
      idxPrev = (idxPrev + 1) % 2;
      idxCurr = (idxCurr + 1) % 2;
    }

    // Get list of matched scans
    unsigned int matchIdx = 0;
    matchedScans.resize(min(run1.size(), run2.size()));
    idxRun1 = run1.size() - 1;
    idxRun2 = run2.size() - 1;
    while (idxRun1 > 0 and idxRun2 > 0) {
      if (predIsDiagonal[idxRun1][idxRun2]) {
        //			cerr<<"D"; cerr.flush();
        matchedScans[matchIdx++].set(idxRun1--, idxRun2--);
        continue;
      }
      if (predIsLeft[idxRun1][idxRun2]) { /*cerr<<"L("<<run1[idxRun1].scan<<","<<run2[idxRun2].scan<<")"; cerr.flush(); */
        idxRun2--;
      }
      else { /*cerr<<"U("<<run1[idxRun1].scan<<","<<run2[idxRun2].scan<<")"; cerr.flush(); */
        idxRun1--;
      }
    }
    matchedScans.resize(matchIdx);
  }

  //
  //  Isotopic envelope classes & functions
  //

  // findEnvelope - finds the index of the envelope in envelopes with the
  //   monoisotopic mass closest to mass
  inline unsigned int IsoEnvelope::findEnvelope(float mass)
  {
    unsigned int lowerLim, upperLim, middle;

    if (envelopes.size() == 0) {
      cerr << "ERROR in findEnvelope (ms1.cpp): No envelope masses loaded!\n";
      exit(-1);
    }

    lowerLim = 0;
    upperLim = envelopes.size() - 1;
    while (upperLim > lowerLim + 1) { // Set lowerLim to the index of the largest monoisotopic mass < mass
      middle = (lowerLim + upperLim) / 2;
      if (envelopes[middle][0] > mass)
        upperLim = middle;
      else
        lowerLim = middle;
    }
    if (fabs(envelopes[lowerLim][0] - mass) < fabs(envelopes[upperLim][0]
        - mass))
      return lowerLim;
    else
      return upperLim;
  }

  //
  //  makeStrict - converts a vector [mass p1 p2 p3] to [mass 0 p1 0 p2 0 p3]
  //
  inline void IsoEnvelope::makeStrict(unsigned int idxEnv, vector<float> &out)
  {
    out.resize(1 + envelopeSize + envelopeSize);
    out[0] = envelopes[idxEnv][0];
    for (unsigned int pivot = 1, counter = 1; counter <= envelopeSize; pivot
        += 2, counter++) {
      out[pivot] = 0;
      out[pivot + 1] = envelopes[idxEnv][counter];
    }
    normalize(out, true); // Need non-zero values in every bin to avoid division errors in the hypothesis test
  }

  float IsoEnvelope::GetMonoisotopicMass(float in_mass, unsigned short charge)
  {
    return (in_mass * ((float) ((short) charge))) - ((float) ((short) charge))
        + AAJumps::massHion;
  }

  //
  //  LoadModel - loads a set of monoisotopic masses and associated average isotopic envelopes.
  //    File is generated in Matlab (getIsotopeTable.m/saveIsotopeTable.m) and the format is:
  //     number of lines, cols in data (2x int32)
  //     isotopic mass, propensity of isotopic peak at isotopic mass + [0 1 2 3 4] (6x float)
  //
  bool IsoEnvelope::LoadModel(const char *filename)
  {
    if (Load_binArray(filename, envelopes) > 0)
      if (envelopes.size() > 0) {
        envelopeSize = envelopes[0].size() - 1;
        for (unsigned int i = 0; i < envelopes.size(); i++)
          normalize(envelopes[i], true, 1);
        return true;
      }
    return false;
  }

  float IsoEnvelope::ScoreEnvelope(float monoisotopicMass,
                                   vector<float> &massEnvelope,
                                   bool strictMode)
  {
    vector<float> *curEnvelope;
    unsigned int idxEnv = findEnvelope(monoisotopicMass), curEnvelopeSize;
    float matchScore = 0;

    if (strictMode) {
      curEnvelope = new vector<float> ;
      makeStrict(idxEnv, *curEnvelope);
      curEnvelopeSize = curEnvelope->size() - 1;
    }
    else {
      curEnvelope = &envelopes[idxEnv];
      curEnvelopeSize = envelopeSize;
    }

    //cerr<<"curEnvelope: "; for(unsigned int i=0; i<curEnvelope->size(); i++) cerr<<(*curEnvelope)[i]<<" "; cerr<<"\n"; cerr.flush();

    if (massEnvelope.size() != curEnvelopeSize) {
      cerr
          << "ERROR in ScoreEnvelope (ms1.cpp): Envelope dimensions don't match!\n";
      exit(-1);
    }
    for (unsigned int pivot = 0; pivot < curEnvelopeSize; pivot++)
      matchScore += massEnvelope[pivot] * log(massEnvelope[pivot]
          / (*curEnvelope)[pivot + 1]);

    if (strictMode)
      delete curEnvelope;

    return matchScore;
  }

  void IsoEnvelope::ExtractEnvelope(float monoisotopicMass,
                                    unsigned short charge,
                                    Spectrum &spec,
                                    float peakWidth,
                                    vector<float> &intensities,
                                    bool strictMode)
  {
    if (strictMode)
      intensities.resize(envelopeSize + envelopeSize);
    else
      intensities.resize(envelopeSize);
    for (unsigned int pivot = 0; pivot < intensities.size(); pivot++)
      intensities[pivot] = 0;

    // Find monoisotopic mass and get its summed intensity
    float increment = 1.0 / ((float) charge);
    vector<int> matches;
    int curSpecIdx, massEnvIdx = 0;
    if (strictMode) {
      // Collect intensity between monoisotopicMass-increment-peakWidth and monoisotopicMass-peakWidth (should total roughly zero for the correct Z/PM)
      curSpecIdx = spec.findClosest(monoisotopicMass - increment - peakWidth);
      if (curSpecIdx < 0)
        return;
      while (curSpecIdx > 0 and spec[curSpecIdx][0] >= monoisotopicMass
          - increment - peakWidth)
        curSpecIdx--;
      curSpecIdx++;
      for (; curSpecIdx < (int) spec.size() and spec[curSpecIdx][0]
          < monoisotopicMass - peakWidth; curSpecIdx++)
        intensities[massEnvIdx] += spec[curSpecIdx][1];
      massEnvIdx++;
    }
    else {
      curSpecIdx = spec.findClosest(monoisotopicMass);
      if (curSpecIdx < 0)
        return;
    }
    spec.findMatches(monoisotopicMass, peakWidth, matches, curSpecIdx);
    for (unsigned int pivot = 0; pivot < matches.size(); pivot++)
      intensities[massEnvIdx] += spec[matches[pivot]][1];
    massEnvIdx++;
    if (matches.size() > 0)
      curSpecIdx = matches[matches.size() - 1] + 1; // Set curSpecIdx to first peak with
    else if (spec[curSpecIdx][0] < monoisotopicMass)
      curSpecIdx++; //  mass > monoisotopicMass+peakWidth

    // Get intensities for the remaining peaks in the envelope
    float nextMass = monoisotopicMass + increment;
    unsigned int envIdx;
    for (envIdx = 1; envIdx < envelopeSize; envIdx++) {
      // Advance to the next peak in the envelope
      for (; curSpecIdx < (int) spec.size() and spec[curSpecIdx][0] < nextMass
          - peakWidth; curSpecIdx++)
        if (strictMode)
          intensities[massEnvIdx] += spec[curSpecIdx][1]; // Add the intensities of the peaks between the current and next isotopic peaks
      if (strictMode)
        massEnvIdx++;

      // Get the intensity for the next peak in the envelope
      for (; curSpecIdx < (int) spec.size() and spec[curSpecIdx][0] <= nextMass
          + peakWidth; curSpecIdx++)
        intensities[massEnvIdx] += spec[curSpecIdx][1];
      massEnvIdx++;
      nextMass += increment;
    }
  }

  void IsoEnvelope::ExtractEnvelope(float monoisotopicMass,
                                    unsigned short charge,
                                    Spectrum &spec,
                                    float peakWidth,
                                    set<int>& ignorePeaks,
                                    vector<float> &intensities,
                                    vector<list<unsigned int> >& indicesUsed,
                                    bool usePPM,
                                    bool strictMode)
  {
    float peakTol = (usePPM) ? (monoisotopicMass * peakWidth * PPM_FACTOR) : peakWidth;
    indicesUsed.clear();
    set<unsigned int> used;
    if (strictMode) {
      intensities.resize(envelopeSize + envelopeSize);
      indicesUsed.resize(envelopeSize + envelopeSize);
    }
    else {
      intensities.resize(envelopeSize);
      indicesUsed.resize(envelopeSize);
    }
    for (unsigned int pivot = 0; pivot < intensities.size(); pivot++)
      intensities[pivot] = 0;

    // Find monoisotopic mass and get its summed intensity
    float increment = 1.0 / ((float) charge);
    vector<int> matches;
    int curSpecIdx, massEnvIdx = 0;
    if (strictMode) {
      // Collect intensity between monoisotopicMass-increment-peakWidth and monoisotopicMass-peakWidth (should total roughly zero for the correct Z/PM)
      curSpecIdx = spec.findClosest(monoisotopicMass - increment - peakTol);
      if (curSpecIdx < 0)
        return;
      while (curSpecIdx > 0 and spec[curSpecIdx][0] >= monoisotopicMass
          - increment - peakTol)
        curSpecIdx--;
      curSpecIdx++;
      for (; curSpecIdx < (int) spec.size() and spec[curSpecIdx][0]
          < monoisotopicMass - peakTol; curSpecIdx++) {
        intensities[massEnvIdx] += spec[curSpecIdx][1];
        if (used.count(curSpecIdx) == 0)
          indicesUsed[massEnvIdx].push_back(curSpecIdx);
        else
          used.insert(curSpecIdx);
      }
      massEnvIdx++;
    }
    else {
      curSpecIdx = spec.findClosest(monoisotopicMass);
      if (curSpecIdx < 0)
        return;
    }

    spec.findMatches(monoisotopicMass, peakTol, matches, curSpecIdx);
    for (unsigned int pivot = 0; pivot < matches.size(); pivot++) {
      intensities[massEnvIdx] += spec[matches[pivot]][1];
      if (used.count(matches[pivot]) == 0)
        indicesUsed[massEnvIdx].push_back(matches[pivot]);
      else
        used.insert(matches[pivot]);
    }
    massEnvIdx++;
    if (matches.size() > 0)
      curSpecIdx = matches[matches.size() - 1] + 1; // Set curSpecIdx to first peak with
    else if (spec[curSpecIdx][0] < monoisotopicMass)
      curSpecIdx++; //  mass > monoisotopicMass+peakWidth

    // Get intensities for the remaining peaks in the envelope
    float nextMass = monoisotopicMass + increment;
    unsigned int envIdx;
    for (envIdx = 1; envIdx < envelopeSize; envIdx++) {
      // Advance to the next peak in the envelope
      peakTol = (usePPM) ? (nextMass * peakWidth * PPM_FACTOR) : peakWidth;
      for (; curSpecIdx < (int) spec.size() and spec[curSpecIdx][0] < nextMass
          - peakTol; curSpecIdx++)
        if (ignorePeaks.count(curSpecIdx) > 0) {
          continue;
        }
        if (strictMode) {
          intensities[massEnvIdx] += spec[curSpecIdx][1]; // Add the intensities of the peaks between the current and next isotopic peaks
          if (used.count(curSpecIdx) == 0)
            indicesUsed[massEnvIdx].push_back(curSpecIdx);
          else
            used.insert(curSpecIdx);
        }
      if (strictMode)
        massEnvIdx++;

      // Get the intensity for the next peak in the envelope
      for (; curSpecIdx < (int) spec.size() and spec[curSpecIdx][0] <= nextMass
          + peakTol; curSpecIdx++) {
        if (ignorePeaks.count(curSpecIdx) > 0) {
          continue;
        }
        intensities[massEnvIdx] += spec[curSpecIdx][1];
        if (used.count(curSpecIdx) == 0)
          indicesUsed[massEnvIdx].push_back(curSpecIdx);
        else
          used.insert(curSpecIdx);
      }
      massEnvIdx++;
      nextMass += increment;
    }
  }

  void IsoEnvelope::normalize(vector<float> &values,
                              bool avoidZeros,
                              unsigned int startIdx)
  {
    float totalIntensity = 0;
    for (unsigned int pivot = startIdx; pivot < values.size(); pivot++) {
      if (avoidZeros and values[pivot] <= 0.0001)
        values[pivot] = 0.0001;
      totalIntensity += values[pivot];
    }
    for (unsigned int pivot = startIdx; pivot < values.size(); pivot++)
      values[pivot] /= totalIntensity;
  }

}
