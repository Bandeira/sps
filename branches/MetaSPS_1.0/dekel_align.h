#ifndef DEKEL_ALIGN_H
#define DEKEL_ALIGN_H

#include <vector>
#include <algorithm>
#include <cmath>

#include "spectrum.h"

//#define DEBUG 1

#ifdef DEBUG

/**
 * TODO: add description
 *
 *@param spec1
 *@param spec2
 *@param peakTol
 *@param matched1
 *@param matched2
 *@param maxAAJump
 *@param sameVertexPenalty
 *@param ptmPenalty
 *@param forceSymmetry
 *@param adZPMmatches
 *@param debug
 *@return
 */
float dekel_align(Spectrum *spec1, Spectrum *spec2, float peakTol,
		Spectrum *matched1, Spectrum *matched2, int maxAAJump,
		float sameVertexPenalty=-5, float ptmPenalty=-5, bool forceSymmetry=true,
		bool addZPMmatches=false, ostream &debug=cerr);
#else

/**
 * TODO: add description
 *
 *@param spec1
 *@param spec2
 *@param peakTol
 *@param matched1
 *@param matched2
 *@param maxAAJump
 *@param sameVertexPenalty
 *@param ptmPenalty
 *@param forceSymmetry
 *@param addZPMmatches
 *@return
 */
float dekel_align(Spectrum *spec1, Spectrum *spec2, float peakTol,
		Spectrum *matched1, Spectrum *matched2, int maxAAJump,
		float sameVertexPenalty = -5, float ptmPenalty = -5,
		bool forceSymmetry = true, bool addZPMmatches = false);
#endif

#endif
