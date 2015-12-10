#include "tags.h"
#include "aminoacid.h"
#include<cmath>

namespace specnets
{
  bool cmp_tag(Tag t1, Tag t2)
  {
    return t1.score > t2.score;
  } // Used to sort tags by descending tag score

  // FindMassNeighbors - indices[i] contains the indices of the peaks in spec
  //    such that abs( mass - jumps[i] - spec[any peak in indices[i]][0] ) <= peakTol
  //    idxSupremum - index of a peak in spec with a mass higher than mass-min(jumps) (optional)
  void FindMassNeighbors(float mass, Spectrum &spec, AAJumps &jumps, vector<
      list<short> > &indices, float peakTol, int idxSupremum = -1)
  {
    if (idxSupremum < 0 or idxSupremum >= (int) spec.size())
      idxSupremum = (int) spec.size() - 1;
    if (idxSupremum < 0)
      return;

    unsigned int jumpIdx;
    int peakIdx = idxSupremum;
    indices.resize(jumps.size());
    for (jumpIdx = 0; jumpIdx < jumps.size(); jumpIdx++)
      indices[jumpIdx].clear();
    if (spec.size() == 0)
      return;
    float curMass;

    for (jumpIdx = 0; jumpIdx < jumps.size(); jumpIdx++) {
      curMass = mass - jumps[jumpIdx];
      if (curMass < -peakTol)
        break;
      if (peakIdx < 0)
        peakIdx = 0;
      while (peakIdx < (int) spec.size() and spec[peakIdx][0] <= curMass
          + peakTol + 0.0001)
        peakIdx++;
      if (peakIdx >= (int) spec.size())
        peakIdx = (int) spec.size() - 1;
      while (peakIdx >= 0 and spec[peakIdx][0] > curMass + peakTol + 0.0001)
        peakIdx--;
      while (peakIdx >= 0 and spec[peakIdx][0] > curMass - peakTol - 0.0001)
        indices[jumpIdx].push_back(peakIdx--);
    }
  }

  unsigned int ExtractTags(Spectrum &spec,
                           list<Tag> &tags,
                           float peakTol,
                           unsigned int tagLen,
                           unsigned int maxNumJumps,
                           unsigned int maxNumTags)
  {
    vector < vector<vector<list<Tag> > > > partialTags; // Table containing all partial tags
    // Position [i,j,k] contains all the tags of length k
    //   with j di-peptide jumps and ending on the i-th spectrum peak
    // dimension 1 - Index of the spectrum peak where the tags end
    // dimension 2 - Number of used di-peptide jumps (usually 0-2)
    // dimension 3 - Tag length
    AAJumps jumps(1);
    unsigned int diIdx, peakIdx, lenIdx, aaIdx, numPeaks = spec.size(),
        numJumps = jumps.size();
    vector<list<short> > aaNeighs(numJumps), // Peaks with mass one aa to the left of current peak
        diNeighs(numJumps); // Peaks with mass two aa to the left of current peak
    Tag curTag(tagLen);

    // Initializations
    partialTags.resize(numPeaks);
    tags.clear();
    for (peakIdx = 0; peakIdx < numPeaks; peakIdx++) {
      partialTags[peakIdx].resize(maxNumJumps + 1);
      for (diIdx = 0; diIdx <= maxNumJumps; diIdx++) {
        partialTags[peakIdx][diIdx].resize(tagLen);
        for (lenIdx = 0; lenIdx < tagLen; lenIdx++)
          partialTags[peakIdx][diIdx][lenIdx].clear();
      }
    }

    // Generate tags
    for (peakIdx = 0; peakIdx < numPeaks; peakIdx++) {
      FindMassNeighbors(spec[peakIdx][0],
                        spec,
                        jumps,
                        aaNeighs,
                        peakTol,
                        (int) peakIdx);
      unsigned int lastAAneigh = peakIdx; // Index of the lowest-mass detected aa neighbor of peakIdx (used for diNeighs searches)

      for (aaIdx = 0; aaIdx < numJumps; aaIdx++) {

        // No neighbors 1-aa away, look for neighbors 2 aa masses away
        if (aaNeighs[aaIdx].size() == 0 and maxNumJumps > 0) {
          FindMassNeighbors(spec[peakIdx][0] - jumps[aaIdx],
                            spec,
                            jumps,
                            diNeighs,
                            peakTol,
                            (int) lastAAneigh);

          for (unsigned int aaDiIdx = 0; aaDiIdx < numJumps; aaDiIdx++) {
            for (list<short>::iterator neigh = diNeighs[aaDiIdx].begin(); neigh
                != diNeighs[aaDiIdx].end(); neigh++) {

              // Initialize tags of length 2 with no middle peak
              curTag.sequence[0] = (char) aaDiIdx;
              curTag.sequence[1] = (char) aaIdx;
              curTag.score = spec[*neigh][1] + spec[peakIdx][1];
              curTag.flankingPrefix = spec[*neigh][0];
              curTag.flankingSuffix = spec.parentMass - AAJumps::massMH
                  - spec[peakIdx][0];
              partialTags[peakIdx][1][1].push_back(curTag);
              if (tagLen == 2)
                tags.push_back(curTag);

              // Extend to tags of length 3+
              for (diIdx = 1; diIdx <= maxNumJumps; diIdx++) {
                for (lenIdx = 2; lenIdx < tagLen; lenIdx++) {
                  if (partialTags[*neigh][diIdx - 1][lenIdx - 2].size() == 0)
                    continue;
                  for (list<Tag>::iterator tagIter = partialTags[*neigh][diIdx
                      - 1][lenIdx - 2].begin(); tagIter
                      != partialTags[*neigh][diIdx - 1][lenIdx - 2].end(); tagIter++) {
                    curTag = *tagIter;
                    curTag.sequence[lenIdx - 1] = (char) aaDiIdx;
                    curTag.sequence[lenIdx] = (char) aaIdx;
                    curTag.score += spec[peakIdx][1];
                    curTag.flankingSuffix = spec.parentMass - AAJumps::massMH
                        - spec[peakIdx][0];
                    partialTags[peakIdx][diIdx][lenIdx].push_back(curTag);
                    if (lenIdx == tagLen - 1)
                      tags.push_back(curTag);
                  }
                }
              }
            }
          }
          continue;
        }

        // Handle neighbors 1-aa away
        for (list<short>::iterator neigh = aaNeighs[aaIdx].begin(); neigh
            != aaNeighs[aaIdx].end(); neigh++) {
          lastAAneigh = *neigh;

          // Initialize tags of length 1
          curTag.sequence[0] = (char) aaIdx;
          curTag.score = spec[*neigh][1] + spec[peakIdx][1];
          curTag.flankingPrefix = spec[*neigh][0];
          curTag.flankingSuffix = spec.parentMass - AAJumps::massMH
              - spec[peakIdx][0];
          partialTags[peakIdx][0][0].push_back(curTag);
          if (tagLen == 1)
            tags.push_back(curTag);

          // Extend to tags of length 2+
          for (diIdx = 0; diIdx <= maxNumJumps; diIdx++) {
            for (lenIdx = 1; lenIdx < tagLen; lenIdx++) {
              if (partialTags[*neigh][diIdx][lenIdx - 1].size() == 0)
                continue;
              for (list<Tag>::iterator tagIter =
                  partialTags[*neigh][diIdx][lenIdx - 1].begin(); tagIter
                  != partialTags[*neigh][diIdx][lenIdx - 1].end(); tagIter++) {
                curTag = *tagIter;
                curTag.sequence[lenIdx] = (char) aaIdx;
                curTag.score += spec[peakIdx][1];
                curTag.flankingSuffix = spec.parentMass - AAJumps::massMH
                    - spec[peakIdx][0];
                partialTags[peakIdx][diIdx][lenIdx].push_back(curTag);
                if (lenIdx == tagLen - 1)
                  tags.push_back(curTag);
              }
            }
          }
        }
      }
    }

    tags.sort(cmp_tag);
    list<Tag>::iterator iter = tags.begin();
    if (maxNumTags > 0 and tags.size() > maxNumTags) {
      for (unsigned int count = 0; count < maxNumTags; count++)
        iter++;
      tags.erase(iter, tags.end());
    }

    return (tags.size());
  }

  unsigned int ExtractTags(char *sequence,
                           vector<Tag> &tags,
                           unsigned int tagLen)
  {
    tags.resize(0);
    if (!sequence or strlen(sequence) < tagLen)
      return 0;
    unsigned int maxStart = strlen(sequence) - tagLen;
    tags.resize(maxStart + 1);
    for (unsigned int tagStart = 0; tagStart <= maxStart; tagStart++) {
      tags[tagStart].sequence.resize(tagLen);
      for (unsigned int pivot = 0; pivot < tagLen; pivot++)
        tags[tagStart].sequence[pivot] = sequence[tagStart + pivot];
    }
    return maxStart + 1;
  }

  /*
   * MatchTagsToSpectra - Matches a set of amino acid sequence tags to each of a set of spectra
   *   and returns the indices and match-score of the matched tags.
   *
   * specs     - set of spectra
   * tags      - set of tags
   * peakTol   - peak mass tolerance when matching tags to the spectrum
   * maxCharge - maximum fragment charge to consider when matching tags to the spectrum
   * missingPeakPenalty - score penalty for missing a peak in the tag
   * noisePeakPenalty   - score penalty factor (times score of noise peaks) applied to unmatched
   *                       spectrum peaks between first and last tag masses
   * maxNumMissedPeaks  - maximum number of missed tag peaks
   * matchedTagsIdx     - indices of matched tags per spectrum
   * matchedTagsScore   - scores of matched tags per spectrum
   * matchedTagsPos     - index of the leftmost spectrum peak where the tag was matched
   */
  void MatchTagsToSpectra(SpecSet &specs,
                          list<Tag> &tags,
                          float peakTol,
                          float maxCharge,
                          float missingPeakPenalty,
                          float noisePeakPenalty,
                          unsigned int maxNumMissedPeaks,
                          vector<list<unsigned int> > &matchedTagsIdx,
                          vector<list<float> > &matchedTagsScores,
                          vector<list<unsigned int> > &matchedTagsPos)
  {
    matchedTagsIdx.resize(specs.size());
    matchedTagsScores.resize(specs.size());
    matchedTagsPos.resize(specs.size());
    for (unsigned int i = 0; i < specs.size(); i++) {
      matchedTagsIdx[i].clear();
      matchedTagsScores[i].clear();
      matchedTagsPos[i].clear();
    }
    if (specs.size() == 0 or tags.size() == 0)
      return;
    vector < vector<float> > tagMasses(tags.size()); // Mass-vector representation of all tags
    vector<float> curTagMasses;
    unsigned int specIdx, peakIdx, tagsIdx, tagMassIdx;
    float tagMatchScore, cumTagMass;

    //	for(tagsIdx=0; tagsIdx<tags.size(); tagsIdx++) getMasses(tags[tagsIdx].sequence, tagMasses[tagsIdx]);
    list<Tag>::iterator tagsIter = tags.begin();
    for (tagsIdx = 0; tagsIter != tags.end(); tagsIter++, tagsIdx++)
      getMasses(tagsIter->sequence, tagMasses[tagsIdx]);

    for (specIdx = 0; specIdx < specs.size(); specIdx++) {
      for (float charge = 1.0; charge < maxCharge + 0.0001; charge
          = (float) round(charge + 1.0)) {
        for (tagsIdx = 0; tagsIdx < tagMasses.size(); tagsIdx++) {
          //			for(tagsIdx=6; tagsIdx<7; tagsIdx++) {
          curTagMasses.resize(tagMasses[tagsIdx].size() + 1);
          curTagMasses[0] = 0;
          for (tagMassIdx = 0; tagMassIdx < tagMasses[tagsIdx].size(); tagMassIdx++)
            curTagMasses[tagMassIdx + 1] = tagMasses[tagsIdx][tagMassIdx]
                / charge;

          for (peakIdx = 0; peakIdx < specs[specIdx].size(); peakIdx++) {
            vector<int> matchesIdx;
            unsigned int pivot, lastMatchedPeak = peakIdx, missedPeaks = 0;
            cumTagMass = specs[specIdx][peakIdx][0];
            tagMatchScore = 0;
            for (tagMassIdx = 0; tagMassIdx < curTagMasses.size()
                and missedPeaks <= maxNumMissedPeaks; tagMassIdx++) {
              cumTagMass += curTagMasses[tagMassIdx];

              // Look for a matching tag peak
              if (specs[specIdx].findMatches(cumTagMass,
                                             peakTol,
                                             matchesIdx,
                                             lastMatchedPeak)) {
                for (pivot = 0; pivot < matchesIdx.size(); pivot++)
                  tagMatchScore += specs[specIdx][matchesIdx[pivot]][1];
                if (tagMassIdx == 0) {
                  lastMatchedPeak = matchesIdx[matchesIdx.size() - 1];
                  continue;
                }
              }
              else {
                tagMatchScore += missingPeakPenalty;
                missedPeaks++;
                //cerr<<"["<<specIdx<<","<<peakIdx<<","<<charge<<"]: missed cumTagMass = "<<cumTagMass<<", tagMassIdx = "<<tagMassIdx<<endl; cerr.flush();
              }

              // Penalize for noise peaks
              for (pivot = lastMatchedPeak + 1; pivot < specs[specIdx].size()
                  and specs[specIdx][pivot][0] < cumTagMass - peakTol - 0.0001; pivot++)
                tagMatchScore += (noisePeakPenalty * specs[specIdx][pivot][1]);
              if (matchesIdx.size())
                lastMatchedPeak = matchesIdx[matchesIdx.size() - 1];
              else {
                for (lastMatchedPeak = pivot - 1; lastMatchedPeak
                    < specs[specIdx].size()
                    and specs[specIdx][lastMatchedPeak][0] < cumTagMass
                        + peakTol + 0.0001; lastMatchedPeak++)
                  ;
                lastMatchedPeak--;
              }
            }

            if (missedPeaks <= maxNumMissedPeaks) {
              //cerr<<"["<<specIdx<<","<<peakIdx<<","<<charge<<"]: missedPeaks = "<<missedPeaks<<", tagsIdx = "<<tagsIdx<<", tagMatchScore = "<<tagMatchScore<<"\n"; cerr.flush();
              matchedTagsIdx[specIdx].push_back(tagsIdx);
              matchedTagsScores[specIdx].push_back(tagMatchScore);
              matchedTagsPos[specIdx].push_back(peakIdx);
            }
            
          } // for (peakIdx = 0;
          
        } // for (tagsIdx = 0;
        
      } // for (float charge = 1.0; 
      
    } // for (specIdx = 0;
    
  } // MatchTagsToSpectra()
  
  
} // namespace specnets

