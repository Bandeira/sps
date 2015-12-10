#include <vector>
#include <algorithm>
#include <cmath>

#include "spectrum.h"

float spec_align(Spectrum *spec1, Spectrum *spec2, float peakTol, 
                  Spectrum *matched1, Spectrum *matched2, int maxAAJump,
                  float sameVertexPenalty=-5, float ptmPenalty=-5, bool forceSymmetry=true,
                  bool addZPMmatches=false);

