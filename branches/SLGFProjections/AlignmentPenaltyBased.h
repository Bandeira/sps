#ifndef ALIGNMENT_PENALTY_BASED
#define ALIGNMENT_PENALTY_BASED

#include "spectrum.h"
#include "PenaltyMatrix.h"
#include <set>

namespace specnets
{
  using namespace std;

  void alignmentPenaltyBased(Spectrum &spec,
                             Spectrum &dbSpec,
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
                             bool enforceEndpeaks = false);
}

#endif

