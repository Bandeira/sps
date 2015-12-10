#include "alignment_scoring.h"
#include "inputParams.h"

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <cmath>


//
//  FindMatchPeaksAll - like findMatchPeaksAll.m. Indices of the matched PRMs are returned in idx1 and idx2.
//
int FindMatchPeaksAll(Spectrum &spec1, Spectrum &spec2, float shift, float tolerance, vector<int> &idx1, vector<int> &idx2){
    int i,j;            // Iterators over the peaks indices
    int low=0,high=0;   // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1

    idx1.resize(0);   idx2.resize(0);
    for(i=0; i<(int)spec1.size(); i++) {
    	while ((low<(int)spec2.size()) && (spec1[i][0]-tolerance-0.000001) > (spec2[low][0]+shift)) low++;
    	while ((high<(int)spec2.size()) && (spec1[i][0]+tolerance+0.000001) >= (spec2[high][0]+shift)) high++;  // high is index of first unreachable peak
    	for (j=low; j<high; j++) { idx1.push_back(i); idx2.push_back(j); }
    }
    return idx1.size();
}

//
//  FindMatchPeaksAll2 - like findMatchPeaksAll.m. Indices of the matched PRMs are returned in idx1 and idx2.
//
int FindMatchPeaksAll2(Spectrum &spec1, Spectrum &spec2, float shift, float tolerance, vector<int> &idx1, vector<int> &idx2){
    int i,j;            // Iterators over the peaks indices
    int low=0,high=0;   // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1

    idx1.resize(0);   idx2.resize(0);
    for(i=0; i<(int)spec1.size(); i++) {
        while (low>0 && (spec1[i][0]-tolerance-0.000001) < (spec2[low][0]+shift)) low--;
        while ((low<(int)spec2.size()) && (spec1[i][0]-tolerance-0.000001) > (spec2[low][0]+shift)) low++;
        while ((high<(int)spec2.size()) && (spec1[i][0]+tolerance+0.000001) >= (spec2[high][0]+shift)) high++;  // high is index of first unreachable peak
        for (j=low; j<high; j++) { idx1.push_back(i); idx2.push_back(j); }
    }
    return idx1.size();
}

//
//  FindMatchPeaks - Finds all _peaks_ that match at least one peak on the
//        other spectrum (i.e. peak masses within tolerance of one another).
//        Note that the output of this function does not have a one-to-one
//        correspondence like FindMatchPeaksAll.
//        Indices of the matched PRMs are returned in idx1 and idx2.
//
//  NOTE: idx1 and idx2 should be pre-allocated to have enough storage and avoid resizing operations in the middle
//
void FindMatchPeaks(Spectrum &spec1, Spectrum &spec2, float shift, float tolerance, vector<int> &idx1, vector<int> &idx2){
    unsigned int i,j;            // Iterators over the peaks indices
    unsigned int low=0,high=0;   // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1
	vector<bool> match1(spec1.size()), match2(spec2.size());
	for(i=0;i<spec1.size();i++) match1[i]=false;   for(i=0;i<spec2.size();i++) match2[i]=false;

    for(i=0; i<spec1.size(); i++) {
        while ((low<spec2.size()) && (spec1[i][0]-tolerance-0.000001) > (spec2[low][0]+shift)) low++;
        while ((high<spec2.size()) && (spec1[i][0]+tolerance+0.000001) >= (spec2[high][0]+shift)) high++;  // high is index of first unreachable peak
        for (j=low; j<high; j++) { match1[i]=true; match2[j]=true; }
    }

    idx1.resize(0);   idx2.resize(0);
	for(i=0;i<spec1.size();i++) if(match1[i]) idx1.push_back(i);
	for(i=0;i<spec2.size();i++) if(match2[i]) idx2.push_back(i);
}

void CyclicAlign(Spectrum &spec1, Spectrum &spec2, float minShift, float minAAmass, float ionOffset,
                    float peakTol, float resolution, float ctermMass, vector<TwoValues<float> > &scoredShifts) {
	unsigned int pivot, peakIdx;
	Spectrum specDouble;    specDouble = spec1;
	unsigned int numPeaks = spec1.size(), numPeaks2=2*numPeaks;
	float massOffset = spec1.parentMass-ctermMass-AAJumps::massHion+ionOffset; // Offset that guarantees that bn overlaps with b0

	specDouble.resize(numPeaks2);   specDouble.parentMass=specDouble.parentMass*2-ctermMass-AAJumps::massHion;
	for(peakIdx=numPeaks; peakIdx<numPeaks2; peakIdx++)
		{ specDouble[peakIdx]=spec1[peakIdx-numPeaks]; specDouble[peakIdx][0]+=massOffset; }
/*	specDouble.resize(3*spec1.size());   specDouble.parentMass=specDouble.parentMass*3-2*(ctermMass+AAJumps::massHion);
	for(peakIdx=0; peakIdx<numPeaks; peakIdx++) { specDouble[peakIdx]=spec1[peakIdx]; specDouble[peakIdx][0]-=massOffset; }
	for(peakIdx=0; peakIdx<numPeaks; peakIdx++) specDouble[numPeaks+peakIdx]=spec1[peakIdx];
	for(peakIdx=0; peakIdx<numPeaks; peakIdx++) { specDouble[numPeaks2+peakIdx]=spec1[peakIdx]; specDouble[numPeaks2+peakIdx][0]+=massOffset; }
*/
	vector<float> shiftsList;   computeShifts2(spec1, spec2, shiftsList, resolution);
	scoredShifts.resize(shiftsList.size());

	vector<int> idx1, idx2, idxMatched1, idxMatched2;
	unsigned int shiftsIdx=0;
	for(pivot=0; pivot<shiftsList.size(); pivot++) {
		if(shiftsList[pivot]<0) continue;
		scoredShifts[shiftsIdx].set(0,shiftsList[pivot]);
		FindMatchPeaksAll2(specDouble, spec2, shiftsList[pivot], peakTol, idx1, idx2);
		for(peakIdx=0; peakIdx<idx1.size(); peakIdx++)
			scoredShifts[shiftsIdx][0]+=specDouble[idx1[peakIdx]][1]*spec2[idx2[peakIdx]][1];
//			scoredShifts[shiftsIdx][0]+=specDouble[idx1[peakIdx]][1]+spec2[idx2[peakIdx]][1];
//		ScoreOverlap6(specDouble, idx1, spec2, idx2, shiftsList[pivot], peakTol, idxMatched1, idxMatched2, minAAmass);
//		for(peakIdx=0; peakIdx<idxMatched1.size(); peakIdx++)
//			scoredShifts[shiftsIdx][0]+=specDouble[idxMatched1[peakIdx]][1]*spec2[idxMatched2[peakIdx]][1];
//			scoredShifts[shiftsIdx][0]+=specDouble[idxMatched1[peakIdx]][1]+spec2[idxMatched2[peakIdx]][1];
		shiftsIdx++;
	}
}

// Helper global variables for computeShifts (to avoid large resizes every time the function is called)

static unsigned int csMaxVecSize=0;       // Allocated size for the arrays below
static TwoValues<float> *shiftScoresV=0;  // TwoValues to hold match scores per spectrum (spec1=[0],spec2=[1])
static TwoValues<int> *shiftMPcount=0;    // Counts the number of matched peaks per shift (spec1=[0],spec2=[1])
static unsigned int *shiftPairsV=0;       // index of the highest scoring symmetric shift
static char* validShift=0;                // Used to signal whether a shift is valid (in terms of start/end-point AA offsets)

void computeShifts_cleanup() {
	if(csMaxVecSize>0)
		{ delete[] shiftScoresV; delete[] shiftMPcount; delete[] shiftPairsV; delete[] validShift; csMaxVecSize=0; }
}

// Vectors are too slow for this computShifts - this is the "gatekeeper" so it must be as fast as possible
/* static bool computeShiftsFirstUse=true;
static vector<TwoValues<float> > shiftScoresV;   // TwoValues to hold match scores per spectrum (spec1=[0],spec2=[1])
static vector<TwoValues<int> > shiftMPcount;     // Counts the number of matched peaks per shift (spec1=[0],spec2=[1])
static vector<unsigned int> shiftPairsV;         // index of the highest scoring symmetric shift
//static vector<bool> validShift;                  // vector<bool> is way too slow to access
static vector<char> validShift;                  // Used to signal whether a shift is valid (in terms of start/end-point AA offsets)
                                                 //   and whether a shift's maximum match score was above the threshold (thus keeping matched peak lists in shiftMatchedPeaks.
*/

//
//  computeShifts - Computes all eligible shifts of spec2 in relation to spec1.
//
//  peakTol, pmTol - peak and parent mass tolerances (typically 0.5 and 1Da)
//  minRatio - minimum acceptable ratio of matched peaks score to total spectrum score (typically 0.3-0.4)
//  minPeakAreaOvlp  - If >0 require the overlaps to span at least this much area (peak area = largest peak mass-smallest peak mass)
//  validJumps       - If != NULL then enforce that shifts and parent mass differences must match some value in validJumps
//
//  shiftScores - list of shift pair scores, sorted by decreasing shift score
//  shiftPairs  - list of shift pairs (shift,shiftSym), both including shiftOffset, in the same order as shiftScores
//  shiftMatchedPeaks - list of matched peaks in (spec1,spec2) for the every shift in any shiftPairs entry. Directly addressable by round(shift mass offset/InputParams::Resolution)
//  bestCandidateScores - maximal pair of matched peak scores. Matched score in spec1/pos[0], spec2/pos[1]
//  bestCandidateMP - maximal number of matched peaks. Num matched peaks in spec1/pos[0], spec2/pos[1], including best symmetric shift
//
// Returns a pair of values:
//   pos[0]: integer shifts offset from the minimum shift considered
//   pos[1]: integer index of the maximum eligible symmetric shift
//
TwoValues<int> computeShifts(Spectrum &spec1, Spectrum &spec2, float peakTol, float pmTol,
                     float minRatio, float minPeakAreaOvlp, int minNumMatchedPeaks,
                     AAJumps &validJumps, list<float> &shiftScores, list<TwoValues<unsigned int> > &shiftPairs,
//                     vector<vector<TwoValues<int> > > &shiftMatchedPeaks,
                     vector<list<TwoValues<int> > > &shiftMatchedPeaks,
                     TwoValues<float> &bestCandidateScores, TwoValues<int> &bestCandidateMP,
                     float minAbsShift, bool addSymmetric) {
	unsigned int szSpec1=spec1.size(), szSpec2=spec2.size(),
	             szVecs=0;  // Size of the vectors used in computeShifts. Determined by minimum/maximum shifts to consider
	float maxPM = max(spec1.parentMass,spec2.parentMass);

	if(csMaxVecSize==0) {
		szVecs = 1+(unsigned int)ceil(2*(maxPM+peakTol+InputParams::PMTol)/InputParams::Resolution);
	    shiftScoresV = new TwoValues<float>[szVecs];   shiftPairsV = new unsigned int[szVecs];
		shiftMPcount = new TwoValues<int>[szVecs];     validShift = new char[szVecs];
		csMaxVecSize = szVecs;
	}
	if(minPeakAreaOvlp>0) minPeakAreaOvlp = minPeakAreaOvlp*maxPM;
	if(spec2.parentMass<minPeakAreaOvlp or spec1.parentMass<minPeakAreaOvlp) return(TwoValues<int>(0,0)); // Test if required overlap areas are feasible

	int shiftsOffset = (int)ceil((spec2.parentMass+peakTol)/InputParams::Resolution);
	szVecs = 1+shiftsOffset+(int)ceil((spec1.parentMass+peakTol)/InputParams::Resolution);
	if(szVecs>csMaxVecSize) { computeShifts_cleanup();
	    shiftScoresV = new TwoValues<float>[szVecs];   shiftPairsV = new unsigned int[szVecs];
		shiftMPcount = new TwoValues<int>[szVecs];     validShift = new char[szVecs];
		csMaxVecSize = szVecs;
	}
	shiftMatchedPeaks.resize(szVecs);
	unsigned int leftShiftIdx = (unsigned int)(minPeakAreaOvlp<=2*peakTol?0:floor((minPeakAreaOvlp-2*peakTol)/InputParams::Resolution));
	unsigned int rightShiftIdx = szVecs-1-leftShiftIdx;

	// Initializations
	for(unsigned int i=0; i<szVecs; i++) { shiftScoresV[i].set(0,0);   shiftMPcount[i].set(0,0); }
	float totalScore1=0; for(unsigned int i=0;i<szSpec1;i++) totalScore1+=spec1[i][1];
	float totalScore2=0; for(unsigned int i=0;i<szSpec2;i++) totalScore2+=spec2[i][1];
	float minScore1 = minRatio*totalScore1, minScore2 = minRatio*totalScore2;
	shiftScores.clear();   shiftPairs.clear(); for(unsigned int i=0; i<szVecs; i++) shiftMatchedPeaks[i].clear(); //.resize(0);
	int shiftIndex, intPeakTol=(int)round(peakTol/InputParams::Resolution), intPMTol=(int)round(pmTol/InputParams::Resolution);

	// Populate shiftScoresV and shiftMatchedPeaks - spectra should have peaks at 0/19/PM-19/PM to match endpoints to internal peaks
	TwoValues<int> peakPair;
	for(unsigned int idxSpec1=0; idxSpec1<szSpec1; idxSpec1++)
		for(unsigned int idxSpec2=0; idxSpec2<szSpec2; idxSpec2++) {
			shiftIndex = shiftsOffset + (int)round((spec1[idxSpec1][0]-spec2[idxSpec2][0])/InputParams::Resolution);
			peakPair.set((int)idxSpec1,(int)idxSpec2);

			for(int tolIdx=-intPeakTol; tolIdx<=intPeakTol; tolIdx++) {
				int shiftIndexTol = shiftIndex+tolIdx; // float curPenalty=abs(tolIdx/10.0);
				if(shiftIndexTol>=0 and shiftIndexTol<szVecs) shiftMatchedPeaks[shiftIndexTol].push_back(peakPair);
			}
		}

	// Add peak scores to all reachable shift positions
	int lastPeak;                      // Index of last peak in spectrum 1 (per shift)
	char *match2 = new char[szSpec2];   // Matched peaks in spectrum 2 (per shift)
	for(unsigned int shiftIndex=leftShiftIdx; shiftIndex<rightShiftIdx; shiftIndex++) {  // Add the scores of all matched peaks without double-counting errors
		if(shiftMatchedPeaks[shiftIndex].size()>0) {
			lastPeak = -1;
			for(unsigned int j=0; j<szSpec2; j++) match2[j]=0;
			list<TwoValues<int> >::iterator matchStart = shiftMatchedPeaks[shiftIndex].begin(),
										    matchEnd = shiftMatchedPeaks[shiftIndex].end();
			for(list<TwoValues<int> >::iterator matchIter = matchStart; matchIter != matchEnd; matchIter++) {
				if((*matchIter)[0]>lastPeak) { shiftScoresV[shiftIndex][0]+=spec1[(*matchIter)[0]][1];   shiftMPcount[shiftIndex][0]++;   lastPeak = (*matchIter)[0]; }
				if(!match2[(*matchIter)[1]]) { shiftScoresV[shiftIndex][1]+=spec2[(*matchIter)[1]][1];   shiftMPcount[shiftIndex][1]++;   match2[(*matchIter)[1]]=1; }
			}
		}
	}
	delete[] match2;

	// Find and mark valid shifts
	float shiftMass=0, maxMass=validJumps.masses[validJumps.size()-1];
//  	for(shiftIndex=0; shiftIndex<(int)validShift.size(); shiftIndex++) validShift[shiftIndex]=true;
  	for(shiftIndex=0; shiftIndex<leftShiftIdx; shiftIndex++) validShift[shiftIndex]=0;
  	for(shiftIndex=rightShiftIdx+1; shiftIndex<szVecs; shiftIndex++) validShift[shiftIndex]=0;
  	for(shiftIndex=leftShiftIdx; shiftIndex<=(int)rightShiftIdx; shiftIndex++) { validShift[shiftIndex]=0;
	  //		if(not shiftMatchedPeaks[shiftIndex].empty()) {
			shiftMass = fabs((shiftIndex-shiftsOffset)*InputParams::Resolution);
			if(shiftMass>=minAbsShift-peakTol and (shiftMass<=peakTol+0.000001 or shiftMass>maxMass+peakTol-0.000001 or validJumps.isValid(shiftMass,peakTol))) {
				validShift[shiftIndex]=1;
			}
	  //		}
	}

	// Create shiftsTmp with entries observing the minRatio and minPeakAreaOvlp restrictions.
    int middleShift = (int)round((spec1.parentMass-spec2.parentMass)/(2*InputParams::Resolution));
    int middleTimesTwo = (int)round((spec1.parentMass-spec2.parentMass)/InputParams::Resolution);
    int upperShiftLimit = shiftsOffset+middleShift;
    int shiftSym;
    list<TwoValues<float> > shiftsTmp;  // List of all curShift pairs (see line below)
    TwoValues<float> curShift;          // Maximum score of the shift pair (pos[0]), index of the base shift (pos[1])

	bestCandidateScores.set(0,0);     bestCandidateMP.set(0,0);
	for(shiftIndex=leftShiftIdx; shiftIndex<=upperShiftLimit; shiftIndex++) {
		if(!validShift[shiftIndex]) continue;
		shiftSym = shiftsOffset+(middleTimesTwo-(shiftIndex-shiftsOffset));

		// Choose highest scoring symmetric shift within tolerance
		float maxScore=0, score1, score2, curPenalty; int maxScoreIdx=-1;
		int numPeaks1, numPeaks2;
		if (addSymmetric) {
			for(int symIdx=max(0,shiftSym-intPMTol); symIdx<=shiftSym+intPMTol and symIdx<(int)szVecs ; symIdx++) {
				if(!validShift[symIdx]) continue;
				score1 = shiftScoresV[shiftIndex][0]+shiftScoresV[symIdx][0];  // +shiftScoresV[shiftSym][0];  Bug fixed 2006/06/19
				score2 = shiftScoresV[shiftIndex][1]+shiftScoresV[symIdx][1];  // +shiftScoresV[shiftSym][1];
				numPeaks1 = shiftMPcount[shiftIndex][0]+shiftMPcount[symIdx][0];
				numPeaks2 = shiftMPcount[shiftIndex][1]+shiftMPcount[symIdx][1];
				if(score1>=minScore1 and score2>=minScore2 and numPeaks1>=minNumMatchedPeaks and numPeaks2>=minNumMatchedPeaks) {
					curPenalty=round(abs(shiftSym-symIdx)*InputParams::Resolution);
					if(score1+score2-curPenalty>maxScore) { maxScore=score1+score2-curPenalty; maxScoreIdx=symIdx; }

					if((score1+score2)>(bestCandidateScores[0]+bestCandidateScores[1])) bestCandidateScores.set(score1,score2);
					if((numPeaks1+numPeaks2)>(bestCandidateMP[0]+bestCandidateMP[1])) bestCandidateMP.set(numPeaks1,numPeaks2);
				}
			}

		} else {
			score1 = shiftScoresV[shiftIndex][0];
			score2 = shiftScoresV[shiftIndex][1];
			numPeaks1 = shiftMPcount[shiftIndex][0];
			numPeaks2 = shiftMPcount[shiftIndex][1];
			if(score1>=minScore1 and score2>=minScore2 and numPeaks1>=minNumMatchedPeaks and numPeaks2>=minNumMatchedPeaks) {
				if((score1+score2)>(bestCandidateScores[0]+bestCandidateScores[1])) bestCandidateScores.set(score1,score2);
				if((numPeaks1+numPeaks2)>(bestCandidateMP[0]+bestCandidateMP[1])) bestCandidateMP.set(numPeaks1,numPeaks2);
				maxScore=score1+score2;
				maxScoreIdx=0;
			}
		}
		if(maxScoreIdx==-1) continue;

		curShift.set(maxScore,shiftIndex);   shiftsTmp.push_back(curShift);
//		shiftPairsV[shiftIndex] = maxScoreIdx;  // This function should not decide the best symmetric shift - that's for the DP function to decide
		shiftPairsV[shiftIndex] = shiftSym;
	}

	// Sort shiftsTmp, create output structures
	shiftsTmp.sort();
	int maxUsedShiftIndex=0;  for(unsigned int i=0;i<szVecs;i++) validShift[i]=0;
	list<TwoValues<float> >::reverse_iterator iter = shiftsTmp.rbegin();
	for(; iter!=shiftsTmp.rend(); iter++) {
		shiftIndex = (int)round((*iter)[1]);   shiftSym = shiftPairsV[shiftIndex];   if(maxUsedShiftIndex<shiftSym) maxUsedShiftIndex=shiftSym;
		shiftScores.push_back((*iter)[0]);
		shiftPairs.push_back(TwoValues<unsigned int>(shiftIndex,shiftSym));
		validShift[shiftIndex]=1;   for(int shiftSymTol=max(0,shiftSym-intPMTol); shiftSymTol<=shiftSym+intPMTol and shiftSymTol<(int)szVecs; shiftSymTol++) validShift[shiftSymTol]=1;
	}

	// Clear list of matched peaks for unused shifts
	for(shiftIndex=0; shiftIndex<szVecs; shiftIndex++)
		if(!validShift[shiftIndex]) shiftMatchedPeaks[shiftIndex].clear(); //.resize(0);

	return TwoValues<int>(shiftsOffset,maxUsedShiftIndex);
}

float computeBestShift(Spectrum& spec1, Spectrum& spec2, float shiftOffset, int minMatchedPeaks, list<TwoValues<int> >& matched) {
	matched.clear();
	
	map<int, list<TwoValues<int> > > shifts;
	list<TwoValues<int> > shiftMP;
	TwoValues<int> MP;
	map<int, float> shiftScore;
	int intPkTol;
	float bestShift;
	float bestScore = 0;
	intPkTol = 2 * (int)(round(InputParams::PeakTol/InputParams::Resolution) + 0.01);

	for (int i = 0; i < spec1.size(); i++) {
		for (int j = 0; j < spec2.size(); j++) {
			float shift = spec1[i][0] - spec2[j][0] + shiftOffset;
			int intShift = (int)(round(shift/InputParams::Resolution) + 0.01);
			MP[0] = i; MP[1] = j;

			for (int sShift = intShift - intPkTol; sShift <= intShift + intPkTol; sShift ++) {
				float score = spec1[i][1] + spec2[j][1] - (((float)abs(sShift-intShift))*0.001);
				if (shiftScore.count(sShift) > 0) {
					shiftScore[sShift] += score;
					shifts[sShift].push_back(MP);
				} else {
					shiftScore[sShift] = score;
					shiftMP.clear();
					shiftMP.push_back(MP);
					shifts[sShift] = shiftMP;
				}
				
				if (shiftScore[sShift] > bestScore && shifts[sShift].size() >= minMatchedPeaks) {
					bestScore = shiftScore[sShift];
					bestShift = ((float)sShift)*InputParams::Resolution;
					matched = shifts[sShift];
				}
			}
		}
	}
	return bestShift;
}

void computeShiftsRaw(Spectrum& spec1, Spectrum& spec2, float shiftOffset, map<int, list<TwoValues<int> > >& bestShifts, unsigned int minNumMatchedPeaks, bool usePPM) {
	bestShifts.clear();
	map<int, list<TwoValues<int> > > shifts;
	list<TwoValues<int> > shiftMP;
	TwoValues<int> MP;
	map<int, float> shiftScore;
	int intPkTol;
	if (!usePPM) {
		intPkTol = 2 * (int)(round(InputParams::PeakTol/InputParams::Resolution) + 0.01);
	}
	
	//cout << "finding shifts ... "; cout.flush();
	for (int i = 0; i < spec1.size(); i++) {
		for (int j = 0; j < spec2.size(); j++) {
			float shift = spec1[i][0] - spec2[j][0] + shiftOffset;
			int intShift = (int)(round(shift/InputParams::Resolution) + 0.01);
			MP[0] = i; MP[1] = j;

			if (usePPM) {
				float shiftol = (InputParams::PPM * shift)/1000000.0;
				float resol = getResolution(shiftol);
				intPkTol = (int)(round(shiftol/resol) + 0.01);
			}

			for (int sShift = intShift - intPkTol; sShift <= intShift + intPkTol; sShift ++) {
				float score = spec1[i][1] + spec2[j][1] - (((float)abs(sShift-intShift))*0.001);
				if (shiftScore.count(sShift) > 0) {
					shiftScore[sShift] += score;
					shifts[sShift].push_back(MP);
				} else {
					shiftScore[sShift] = score;
					shiftMP.clear();
					shiftMP.push_back(MP);
					shifts[sShift] = shiftMP;
				}
			}
		}
	}

	/*
	cout << "\n\nShifts\n";
	int mx = 0;
	for (map<int, float>::iterator shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt ++) {
		int shift = shiftScoreIt->first;
		float score = shiftScoreIt->second;
		cout << shift << " - " << score << " : ";
		if (shifts[shift].size() > mx) mx = shifts[shift].size();
		for (list<TwoValues<int> >::iterator it = shifts[shift].begin(); it != shifts[shift].end(); it++) {
			cout << spec1[(*it)[0]][0] << "," << spec2[(*it)[1]][0] << "(" << (*it)[0] << "," << (*it)[1] << "); ";
		}
		cout << "\n";
	}

	cout << "\nMAX = " << mx << "\n\n";
	*/

	// only keep shifts that are at the center of their resolution distribution
	float lastScore = -1.0;
	int lastShift;
	bool increasing = true, decreasing = false, peaking = false;

	list<int> prevShifts;
	//cout << "finished, picking top shifts ... "; cout.flush();
	for (map<int, float>::iterator shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt ++) {
		int shift = shiftScoreIt->first;
		float score = shiftScoreIt->second;
		
		if (shifts[shift].size() < minNumMatchedPeaks) {increasing = true; decreasing = false; peaking = false; continue;}
		
		if (lastScore > 0 && score > lastScore && lastShift + 1 == shift) {increasing = true; decreasing = false; peaking = false;}
		else if (lastScore > 0 && score < lastScore && lastShift + 1 == shift) {
			if (increasing || peaking) {prevShifts.push_back(lastShift);}
			decreasing = true; increasing = false; peaking = false;
		}
		else if (lastScore > 0 && score == lastScore && lastShift + 1 == shift) {increasing = false; decreasing = false; peaking = true; prevShifts.push_back(lastShift);}
		else {
			if (increasing || peaking) {prevShifts.push_back(lastShift);}
			increasing = true; decreasing = false; peaking = false;
		}

		if (! peaking && prevShifts.size() > 0) {
			float bestShift = 0;
			for (list<int>::iterator lit = prevShifts.begin(); lit != prevShifts.end(); lit++) {
				bestShift += *lit;
			}
			bestShift /= (float)prevShifts.size();
			int bestShiftInt = (int)(round(bestShift)+0.01);
			bestShifts[bestShiftInt] = shifts[bestShiftInt];
			prevShifts.clear();
		}

		lastScore = score;
		lastShift = shift;
	}

	/*
	cout << "\n\nBest Shifts\n";
	mx = 0;
	for (map<int, list<TwoValues<int> > >::iterator shiftScoreIt = bestShifts.begin(); shiftScoreIt != bestShifts.end(); shiftScoreIt ++) {
		int shift = shiftScoreIt->first;
		float score = shiftScore[shift];
		cout << shift << " - " << score << " : ";
		if (bestShifts[shift].size() > mx) mx = bestShifts[shift].size();
		for (list<TwoValues<int> >::iterator it = bestShifts[shift].begin(); it != bestShifts[shift].end(); it++) {
			cout << spec1[(*it)[0]][0] << "," << spec2[(*it)[1]][0] << "(" << (*it)[0] << "," << (*it)[1] << "); ";
		}
		cout << "\n";
	}

	cout << "\nMAX = " << mx << "\n\n";
	cout.flush();
	*/
	//cout << "finished, "; cout.flush();
}

// Helper global variables for computeShifts2
static bool computeShifts2firstUse=true;
void computeShifts2(Spectrum &spec1, Spectrum &spec2, vector<float> &shiftsList, float resolution) {
	if(computeShifts2firstUse) {computeShifts2firstUse=false;}
	shiftsList.resize(spec1.size()*spec2.size()+1);
	unsigned int shiftsIdx=0, spec1idx, spec2idx;
	shiftsList[shiftsIdx++]=0;
//	for(spec1idx=0; spec1idx<spec1.size(); spec1idx++) shiftsList[shiftsIdx++]=spec1[spec1idx][0]; // One peak is always matched
//	for(spec1idx=0; spec1idx<spec1.size(); spec1idx++) shiftsList[shiftsIdx++]=spec1[spec1idx][0]-spec2.parentMass; // One peak is always matched
	for(spec1idx=0; spec1idx<spec1.size(); spec1idx++)
		for(spec2idx=0; spec2idx<spec2.size(); spec2idx++)
			shiftsList[shiftsIdx++]=(spec1[spec1idx][0]-spec2[spec2idx][0]);

	// Remove duplicate shift entries
	sort(shiftsList.begin(), shiftsList.end());
	unsigned int shiftsIdxUnique=0;
	for(shiftsIdx=1; shiftsIdx<shiftsList.size(); shiftsIdx++)
		if(fabs(shiftsList[shiftsIdx]-shiftsList[shiftsIdxUnique])>resolution) {
			shiftsIdxUnique++;
			if(shiftsIdx>shiftsIdxUnique) shiftsList[shiftsIdxUnique]=shiftsList[shiftsIdx];
		}
	shiftsList.resize(shiftsIdxUnique+1);
}

//
//  ScoreOverlap6 - Front-end function for ScoreOverlap6
//
//  minIPdist - minimum allowed mass distance between 2 consecutive peaks (defaults to AAJumps::minAAmass)
//  idxMatched       - Indices of the matched PRMs (sparse set), col 0 for spec1 and col 1 for spec2
//
float ScoreOverlap6(Spectrum &spec1, Spectrum &spec2, float shift, float tolerance,
					 vector<int> &idxMatched1, vector<int> &idxMatched2,
					 float minIPdist, float *offsetPenalty){
	vector<int> idx1all, idx2all;
	FindMatchPeaksAll(spec1, spec2, shift, tolerance, idx1all, idx2all);
	return ScoreOverlap6mp(spec1,idx1all,spec2,idx2all,shift,tolerance,idxMatched1,idxMatched2,minIPdist,offsetPenalty);
}

//
//  ScoreOverlap6mp - like ScoreOverlapE.m with distance function 6. Min distance between PRMs is 57 - 2*tolerance
//
//  minIPdist - minimum allowed mass distance between 2 consecutive peaks (defaults to AAJumps::minAAmass)
//  idx1all, idx2all - as returned by FindMatchPeaksAll
//  idxMatched       - Indices of the matched PRMs (sparse set), col 0 for spec1 and col 1 for spec2
//
float ScoreOverlap6mp(Spectrum &spec1, vector<int> idx1all, Spectrum &spec2, vector<int> idx2all, float shift,
					 float tolerance, vector<int> &idxMatched1, vector<int> &idxMatched2,
					 float minIPdist, float *offsetPenalty){
    vector< TwoValues<float> > values(min(idx1all.size(),idx2all.size())+1);  // Keeps the values for the dynamic programming recursion: predecessor (col 0) and predecessor score (col 1)
                                                          // +1 because first line is (0,0) for DP initialization
    int forbiddenIdx;  // Index of the first forbidden PRM (because it is too close to the current PRM
    int maxIdx;        // Index of the best PRM that the current PRM can connect to
    float bestMatch;   // Best path score so far
    int bestMatchIdx;  // Index of the last PRM in the best path so far
    int   i,j;         // Iterator vars

    idxMatched1.resize(0);   idxMatched2.resize(0);

    values[0][0]=0;   values[0][1]=0;
    maxIdx = 0;       forbiddenIdx = 1;
    bestMatch = 0;    bestMatchIdx = 0;
    for (i=1; i<(int)values.size(); i++) {
        while (forbiddenIdx<i and
               spec1[idx1all[forbiddenIdx-1]][0] <= (spec1[idx1all[i-1]][0]-minIPdist+2*tolerance) and
               spec2[idx2all[forbiddenIdx-1]][0] <= (spec2[idx2all[i-1]][0]-minIPdist+2*tolerance)) {
            // This is executed only when forbidden is advanced
            if (values[forbiddenIdx][1] > values[maxIdx][1]) { maxIdx = forbiddenIdx; }
            forbiddenIdx++;
        }
        values[i][0] = maxIdx;
        values[i][1] = values[maxIdx][1] + spec1[idx1all[i-1]][1] + spec2[idx2all[i-1]][1];
        if(offsetPenalty) values[i][1]-= abs(spec1[idx1all[i-1]][0] - (spec2[idx2all[i-1]][0]+shift));
        if (values[i][1]>bestMatch) { bestMatch=values[i][1];  bestMatchIdx=i; }  // Keep track of where the best path ends
    }

    list<int> bestPath;
    while (bestMatchIdx>0) {
    	if(offsetPenalty) (*offsetPenalty)+= abs(spec1[idx1all[bestMatchIdx-1]][0] - (spec2[idx2all[bestMatchIdx-1]][0]+shift));
    	bestPath.push_back(bestMatchIdx-1);  bestMatchIdx=(int)values[bestMatchIdx][0];
    }

    // ******************************************
    //    Populate final lists of matched PRMs
    // ******************************************
    unsigned int bestPathLength = bestPath.size();
    idxMatched1.resize(bestPathLength);
    idxMatched2.resize(bestPathLength);

    for(unsigned int idxPath=0; idxPath<bestPathLength; idxPath++) {
        idxMatched1[idxPath] = idx1all[bestPath.back()];
        idxMatched2[idxPath] = idx2all[bestPath.back()];
        bestPath.pop_back();
    }

    return bestMatch;
}

//
//  ScoreOverlap7 - like ScoreOverlapE.m with distance function 7 but gets matching PRMs from findMatchPeaksAll instead of findMatchPeaks2a. Min distance between PRMs is 57 - 2*tolerance
//
//  idx1all, idx2all - as returned by FindMatchPeaksAll
//  idxMatched1,2    - Indices of the matched PRMs (sparse set) for spec1 and for spec2
//  symmetryOffset   - Masses of symmetric peaks should add up to sum(peptide masses)+18+symmetryOffset
//
float ScoreOverlap7(Spectrum &spec1, vector<int> idx1all, Spectrum &spec2, vector<int> idx2all, float shift, float tolerance,
                    vector<int> &idxMatched1, vector<int> &idxMatched2, float symmetryOffset){
	Spectrum tmpSpec;   vector<int> idxMatched;   float score;

// Use FindMatchPeaksAll results to construct a single spectrum with all possible peak matches
	tmpSpec.parentMass = (spec1.parentMass+spec2.parentMass)/2;
	tmpSpec.peakList.resize(idx1all.size());
	tmpSpec.idDist = spec1.idDist;

//cerr<<"tmpSpec: (parent mass = "<<tmpSpec.parentMass<<")\n";
	for(unsigned int i=0; i<tmpSpec.size(); i++) {
		tmpSpec[i].set((spec1[idx1all[i]][0]+spec2[idx2all[i]][0]+shift)/2,spec1[idx1all[i]][1]+spec2[idx2all[i]][1]);
//cerr<<i<<": "<<spec1[idx1all[i]][0]<<"+"<<spec2[idx2all[i]][0]+shift<<" -> ["<<tmpSpec[i][0]<<","<<tmpSpec[i][1]<<"]\n";
	}


// Use getMaxSparseSet to select which peak matches to keep and determine match score
//	score = getMaxSparseSet(tmpSpec, tolerance, 0, idxMatched);
	score = getMaxSparseSet(tmpSpec, tolerance, symmetryOffset, idxMatched, true);

/*cerr << "tmpSpec's matched peaks:\n";
for(int i=0;i<idxMatched.size();i++)
	cerr <<idxMatched[i]<<"\t"<<tmpSpec[idxMatched[i]][0]<<"\t"<<tmpSpec[idxMatched[i]][1]<<"\n";
*/
	idxMatched1.resize(idxMatched.size());   idxMatched2.resize(idxMatched.size());
	for(unsigned int i=0; i<idxMatched.size(); i++)
		{ idxMatched1[i]=idx1all[idxMatched[i]]; idxMatched2[i]=idx2all[idxMatched[i]]; }

	return score;
}

//
//  getMaxSparseSet - like ScoreOverlap7 but takes only one spectrum as input
//
//  pmOffset=0 for PRM spectra and pmOffset=2 for MS/MS spectra.
//  idxMatched - Indices of the matched PRMs (sparse set)
//  includeSymmetric - set to true if the sparse set should also include the symmetric
//                      peaks of the matched peaks (e.g. ScoreOverlap7)
//
//  NOTE: Make sure that spec.idDist is adequately set.
//
float getMaxSparseSet(Spectrum &spec, float tolerance, float pmOffset, vector<int> &idxMatched, bool includeSymmetric){
	vector<vector<TwoValues<int> > > prevPair,     // Previous pair for every possible pair
	                                   bestPair;     // Best pair so far
	vector<vector<float> > scores,         // Scores for every possible path
	                       bestScores;     // Best path score so far
	vector<float> peakScores;    // Adjusted peak scores (including the score of the symmetric peaks)
	vector<int> peakPairs;     // Indices of the other paired peak (if any) used in peakScores
	vector<TwoValues<int> > prefixPRMs;  // PRMs with mass<=aaMass/2 and sorted by increasing mass
	                                       //  Col 1 is PRM index, col 2 is index of closest PRM >=57-2*tolerance Da away
	vector<TwoValues<int> > suffixPRMs;  // PRMs with mass>aaMass/2 and sorted by increasing distance to aaMass
	TwoValues<int> globalBestPair(0,0);
	float aaMass = spec.parentMass+(pmOffset-spec.idDist)*AAJumps::massHion,
		  minInterPeakDist = 57*spec.idDist-2*tolerance,
		  globalBestScore=0;
	int i,j,k,p,idxPref,idxSuff;

	if(spec.size()==0) { idxMatched.resize(0); return 0; }

	// Make sure that the spectrum has PRMs at zero and at aaMass
	short addZero=0, addPM=0;
	Spectrum oldSpec;           // oldSpec is used to keep a copy of the input spectrum
	                            //   whenever addZero or addPM are >0
	if(spec[0][0]>tolerance) addZero++;
	if(spec[spec.size()-1][0]<aaMass-AAJumps::massH2O*spec.idDist-tolerance) addPM++;
	if(addZero+addPM>0) {
		oldSpec = spec;
		spec.peakList.resize(spec.size()+addZero+addPM);
		if (addZero) { for(i=spec.size()-1-addPM; i>=1; i--) spec[i]=spec[i-1]; spec[0].set(0,0); }
		if (addPM) spec[spec.size()-1].set(aaMass-AAJumps::massH2O*spec.idDist,0);
	}

//for(unsigned int peakIdx=0; peakIdx<spec.size(); peakIdx++)
//	cerr<<peakIdx<<": "<<spec[peakIdx][0]<<", "<<spec[peakIdx][1]<<endl;

   	// Populate prefixPRMs and suffixPRMs - peak index (col.1), predecessor index (col.2)
	prefixPRMs.resize(spec.size());
	for(i=0; i<(int)spec.size() && spec[i][0]<=aaMass/2; i++) {
		for(p=i-1; p>=0 && spec[i][0]-spec[p][0]<minInterPeakDist; p--);
		prefixPRMs[i].set(i,max((int)p,0));
	}
	prefixPRMs.resize(i);   suffixPRMs.resize(spec.size()-prefixPRMs.size());
	suffixPRMs[suffixPRMs.size()-1].set(spec.size()-1,0);
	for (j=spec.size()-1, k=0; j>=0 && spec[j][0]>aaMass/2; j--) {
		suffixPRMs[k][0]=j;
		for(p=0; suffixPRMs[p][0]>j && spec[suffixPRMs[p][0]][0]-spec[j][0]>minInterPeakDist; p++);
		suffixPRMs[k++].set(j,max(p-1,0));
	}

	// Resize and initialize all the DP variables
	prevPair.resize(prefixPRMs.size());   bestPair.resize(prefixPRMs.size());
	scores.resize(prefixPRMs.size());     bestScores.resize(prefixPRMs.size());
	for(int i=0; i<(int)prefixPRMs.size(); i++) {
		prevPair[i].resize(suffixPRMs.size());   bestPair[i].resize(suffixPRMs.size());
		scores[i].resize(suffixPRMs.size());     bestScores[i].resize(suffixPRMs.size());
		for(j=0; j<suffixPRMs.size(); j++) { scores[i][j]=0; bestScores[i][j]=0; }
	}

	// Consolidate the scores of symmetric peaks
	peakScores.resize(spec.size());   peakPairs.resize(spec.size());
	for(k=0;k<spec.size();k++) { peakScores[k]=spec[k][1]; peakPairs[k]=-1; }
	vector<vector<float> > pairs; vector<vector<int> > pairsIdx;
	spec.getPairs(pmOffset, tolerance, pairs, pairsIdx);
	for(k=0;k<pairsIdx.size();k++) {
		if(pairsIdx[k][0]<0 || pairsIdx[k][1]<0) continue;
		if(peakScores[pairsIdx[k][0]]<pairs[k][2]) { peakScores[pairsIdx[k][0]]=pairs[k][2]; peakPairs[pairsIdx[k][0]]=pairsIdx[k][1]; }
		if(peakScores[pairsIdx[k][1]]<pairs[k][2]) { peakScores[pairsIdx[k][1]]=pairs[k][2]; peakPairs[pairsIdx[k][1]]=pairsIdx[k][0]; }
	}


	//
	// DP to compute maximum sparse set
	//
	scores[0][0]=peakScores[prefixPRMs[0][0]]+peakScores[suffixPRMs[0][0]];
	bestScores[0][0] = scores[0][0];  bestPair[0][0].set(0,0);
	globalBestScore = scores[0][0];   globalBestPair.set(0,0);
	for(i=0; i<prefixPRMs.size(); i++)
		for(j=0; j<suffixPRMs.size(); j++) {
			if ((i==0 && j==0) || spec[suffixPRMs[j][0]][0]-spec[prefixPRMs[i][0]][0]<minInterPeakDist) continue;

			// Set default values of best scores/pairs for position [i][j]
			if(i>0) { bestScores[i][j]=bestScores[i-1][j];  bestPair[i][j]=bestPair[i-1][j]; }
			if(j>0 && bestScores[i][j-1]>bestScores[i][j]) { bestScores[i][j]=bestScores[i][j-1];  bestPair[i][j]=bestPair[i][j-1]; }

			idxPref = prefixPRMs[i][0];   idxSuff = suffixPRMs[j][0];
			// j ranges over suffixes whose masses differ from spec[i][0]
			if (fabs(spec[idxPref][0]+spec[idxSuff][0]-aaMass)>2*tolerance) {
//cerr<<" --- "<<spec[idxPref][0]+spec[idxSuff][0]<<", "<<aaMass<<endl;
				if (spec[idxPref][0]>aaMass-spec[idxSuff][0]) {  // last jump was on the prefix side
					scores[i][j] = peakScores[idxPref] + bestScores[prefixPRMs[i][1]][j];
					prevPair[i][j] = bestPair[prefixPRMs[i][1]][j];
//cerr<<"["<<i<<","<<j<<"] jumping from ["<<prefixPRMs[i][1]<<","<<j<<"], score "<<scores[i][j]<<"\n";
				} else {    // last jump was on the suffix side
					scores[i][j] = peakScores[idxSuff] + bestScores[i][suffixPRMs[j][1]];
					prevPair[i][j] = bestPair[i][suffixPRMs[j][1]];
//cerr<<"["<<i<<","<<j<<"] jumping from ["<<i<<","<<suffixPRMs[j][1]<<"], score "<<scores[i][j]<<"\n";
				}
			} else {  // still consider these pairs but don't increase the score
				scores[i][j] = bestScores[prefixPRMs[i][1]][j];
				prevPair[i][j] = bestPair[prefixPRMs[i][1]][j];
				if(scores[i][j] < bestScores[i][suffixPRMs[j][1]]) {
					scores[i][j] = bestScores[i][suffixPRMs[j][1]];
					prevPair[i][j] = bestPair[i][suffixPRMs[j][1]];
				}
//cerr<<"["<<i<<","<<j<<"] (same!) jumping from ["<<i<<","<<suffixPRMs[j][1]<<"] or ["<<prefixPRMs[i][1]<<","<<j<<"], score "<<scores[i][j]<<"\n";
			}

			if(scores[i][j]>bestScores[i][j]) { bestScores[i][j] = scores[i][j];  bestPair[i][j].set(i,j); }
			if(bestScores[i][j]>globalBestScore) { globalBestScore=bestScores[i][j];  globalBestPair.set(i,j); }
		}

//cerr<<"Global best score is "<<globalBestScore<<", pair ["<<globalBestPair[0]<<","<<globalBestPair[1]<<"]\n";

	// Construct idxMatched
	TwoValues<int> tmpPair;     int curPeakIdx;
	vector<bool> idxMatchedBool;   // Boolean vector used to mark matched peaks
	idxMatchedBool.resize(spec.size());   for(int i=0; i<spec.size(); i++) idxMatchedBool[i]=false;
	if(globalBestPair[0]>0) {
		idxMatchedBool[prefixPRMs[globalBestPair[0]][0]-addZero]=true;
		if(includeSymmetric and peakPairs[prefixPRMs[globalBestPair[0]][0]]>0) idxMatchedBool[peakPairs[prefixPRMs[globalBestPair[0]][0]]-addZero]=true;
	}
	if(globalBestPair[1]>0) {
		idxMatchedBool[suffixPRMs[globalBestPair[1]][0]-addZero]=true;
		if(includeSymmetric and peakPairs[suffixPRMs[globalBestPair[1]][0]]>0) idxMatchedBool[peakPairs[suffixPRMs[globalBestPair[1]][0]]-addZero]=true;
	}
	while (globalBestPair[0]>0 || globalBestPair[1]>0) {
		tmpPair = prevPair[globalBestPair[0]][globalBestPair[1]];
		curPeakIdx=-1;
		if(tmpPair[0]!=globalBestPair[0])
			{ if(tmpPair[0]>0) curPeakIdx=prefixPRMs[tmpPair[0]][0]; }
		else if(tmpPair[1]>0) curPeakIdx=suffixPRMs[tmpPair[1]][0];
		if (curPeakIdx>0) {
			idxMatchedBool[curPeakIdx-addZero]=true;
			if (includeSymmetric and peakPairs[curPeakIdx]>0) idxMatchedBool[peakPairs[curPeakIdx]-addZero]=true;
//cerr << "Marked ["<<curPeakIdx-addZero<<","<<peakPairs[curPeakIdx]-addZero<<"]\n";
		}
		globalBestPair = tmpPair;
	}
	if(addZero==0) idxMatchedBool[prefixPRMs[0][0]]=true;          // If these were already in the spectrum
	if(addPM==0) idxMatchedBool[suffixPRMs[0][0]-addZero]=true;    //  then their scores are part of the match

	// Copy matches from idxMatchedBool to idxMatched
	idxMatched.resize(spec.size()); k=0;
	for(int i=0; i<spec.size(); i++) if(idxMatchedBool[i]) idxMatched[k++]=i;
	idxMatched.resize(k);

	if (addZero+addPM>0) spec=oldSpec;  // Reverse the addition of peaks at masses zero and parentMass

	return globalBestScore;
}
