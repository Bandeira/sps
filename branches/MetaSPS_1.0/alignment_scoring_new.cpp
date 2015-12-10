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
//  NOTE: idx1 and idx2 should be pre-allocated to have enough storage and avoid resizing operations in the middle
//
short FindMatchPeaksAll(Spectrum &spec1, Spectrum &spec2, float shift, float tolerance, vector<short> &idx1, vector<short> &idx2){
    short i,j;            // Iterators over the peaks indices
    short low=0,high=0;   // Index bounds of the peaks in spec2 the lie within tolerance of current peak in spec1

    idx1.resize(0);   idx2.resize(0);
//    while (low<(int)spec2.size() && spec2[low][0]<=57.0-tolerance) low++;  // Don't match peaks under 57-tolerance Daltons - THIS IS NOT THE PLACE TO CHECK THIS: if these peaks are not supposed to be matched then remove these peaks from the spectra before calling this function
    for(i=0; i<(int)spec1.size(); i++) {
//        if (spec1[i][0] >= 57.0-tolerance) {
            while ((low<(int)spec2.size()) && (spec1[i][0]-tolerance-0.0001) > (spec2[low][0]+shift)) low++;
            while ((high<(int)spec2.size()) && (spec1[i][0]+tolerance+0.0001) >= (spec2[high][0]+shift)) high++;  // high is index of first unreachable peak
            for (j=low; j<high; j++) { idx1.push_back(i); idx2.push_back(j); }
//        }
    }
    return idx1.size();
}

// Helper global variables for computeShifts (to avoid large resizes every time the function is called)

static bool computeShiftsFirstUse=true;
static vector<TwoValues<float> > shiftScoresV;
static vector<unsigned int> shiftPairsV;
static vector<TwoValues<short> > maxPeakUsed;  // Used to prevent having the same peak scores repeatedly added to the same shift positions (because of tolerance window overlaps)
static vector<bool> validShift, usedShift;                
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
//  shiftMatchedPeaks - list of matched peaks in (spec1,spec2) for the every shift in any shiftPairs entry
//

TwoValues<int> computeShifts(Spectrum &spec1, Spectrum &spec2, float peakTol, float minRatio,
                     float minPeakAreaOvlp, AAJumps &validJumps,
                     list<float> &shiftScores, list<TwoValues<unsigned int> > &shiftPairs, 
                     vector<vector<TwoValues<short> > > &shiftMatchedPeaks, TwoValues<float> &bestCandidateScores) {
	if(computeShiftsFirstUse) {
		unsigned int szVecs = 1+(unsigned int)ceil(2*(InputParams::MaxShift+InputParams::PeakTol+InputParams::PMTol)/InputParams::Resolution);
		shiftScoresV.resize(szVecs);   shiftPairsV.resize(szVecs);   maxPeakUsed.resize(szVecs);   validShift.resize(szVecs);   usedShift.resize(szVecs);
		computeShiftsFirstUse=false; 
	}
    if(minPeakAreaOvlp>0 and (spec2.parentMass<minPeakAreaOvlp*spec1.parentMass or spec1.parentMass<minPeakAreaOvlp*spec2.parentMass)) return(TwoValues<int>(0,0)); // Test if required overlap areas are feasible

//	int shiftsOffset = (int)ceil((spec2.parentMass+4*InputParams::PeakTol)/InputParams::Resolution);
//	shiftScoresV.resize(shiftsOffset+(int)ceil((spec1.parentMass+4*InputParams::PeakTol)/InputParams::Resolution));        // TwoValues to hold match scores per spectrum (spec1=[0],spec2=[1])
	int shiftsOffset = (int)ceil((min(spec2.parentMass*(1-minPeakAreaOvlp),InputParams::MaxShift)+InputParams::PeakTol)/InputParams::Resolution);
	shiftScoresV.resize(shiftsOffset+(int)ceil(min(spec1.parentMass*(1-minPeakAreaOvlp),InputParams::MaxShift)/InputParams::Resolution));        // TwoValues to hold match scores per spectrum (spec1=[0],spec2=[1])
	shiftPairsV.resize(shiftScoresV.size());                        // index of the highest scoring symmetric shift
	shiftMatchedPeaks.resize(shiftScoresV.size()); // Directly addressable by round(shift mass offset/InputParams::Resolution)
	validShift.resize(shiftScoresV.size());   usedShift.resize(shiftScoresV.size());	
	
	// Initializations
	bestCandidateScores.set(0,0);
	for(unsigned int i=0; i<shiftScoresV.size(); i++) shiftScoresV[i].set(0,0);
	float totalScore1=0; for(unsigned int i=0;i<spec1.size();i++) totalScore1+=spec1[i][1];
	float totalScore2=0; for(unsigned int i=0;i<spec2.size();i++) totalScore2+=spec2[i][1];
	float minScore1 = minRatio*totalScore1, minScore2 = minRatio*totalScore2;
	shiftScores.clear();   shiftPairs.clear(); for(unsigned int i=0; i<shiftMatchedPeaks.size(); i++) shiftMatchedPeaks[i].resize(0);

	// Populate shiftScoresV and shiftMatchedPeaks - spectra should have peaks at 0/19/PM-19/PM to match endpoints to internal peaks
	int shiftIndex, intPeakTol=(int)round(peakTol/InputParams::Resolution);

	TwoValues<short> peakPair;

	// Add peak scores (and offset penalties) from spec1 to all reachable shift positions
	for(unsigned int i=0; i<maxPeakUsed.size();i++) maxPeakUsed[i].set(-1,-1);
	for(unsigned int idxSpec1=0; idxSpec1<spec1.size(); idxSpec1++)
		for(unsigned int idxSpec2=0; idxSpec2<spec2.size(); idxSpec2++) {
			shiftIndex = shiftsOffset + (int)round((spec1[idxSpec1][0]-spec2[idxSpec2][0])/InputParams::Resolution);
			peakPair.set((short)idxSpec1,(short)idxSpec2);
			
//			if(shiftIndex+intPeakTol>=shiftScoresV.size()) { shiftScoresV.resize(shiftIndex+intPeakTol+100); shiftPairsV.resize(shiftIndex+intPeakTol+100); maxPeakUsed.resize(shiftIndex+intPeakTol+100);}
			
			for(int tolIdx=-intPeakTol; tolIdx<=intPeakTol; tolIdx++) {
				int shiftIndexTol = shiftIndex+tolIdx; // float curPenalty=abs(tolIdx/10.0);
				if(shiftIndexTol>=0 and shiftIndexTol<shiftMatchedPeaks.size()) {
					shiftMatchedPeaks[shiftIndexTol].push_back(peakPair); 
	
					if(maxPeakUsed[shiftIndexTol][0]<(short)idxSpec1) 
						{ shiftScoresV[shiftIndexTol][0]+=spec1[idxSpec1][1]; maxPeakUsed[shiftIndexTol][0]=(short)idxSpec1; }
					if(maxPeakUsed[shiftIndexTol][1]<(short)idxSpec2) 
						{ shiftScoresV[shiftIndexTol][1]+=spec2[idxSpec2][1]; maxPeakUsed[shiftIndexTol][1]=(short)idxSpec2; }
				}
			}
		}

	float shiftMass=0, maxMass=validJumps.masses[validJumps.size()-1];
	for(shiftIndex=0; shiftIndex<(int)validShift.size(); shiftIndex++) { validShift[shiftIndex]=false;
		if(shiftMatchedPeaks[shiftIndex].size()>0) {
			shiftMass = fabs((shiftIndex-shiftsOffset)*InputParams::Resolution);
			if(fabs(shiftMass)<=peakTol or fabs(shiftMass)>maxMass or validJumps.isValid(shiftMass,peakTol)) validShift[shiftIndex]=true;
		}
	}
	
	// Create shiftsTmp with entries observing the minRatio and minPeakAreaOvlp restrictions.
/*    int middleShift = (int)round(10*(spec1.parentMass-spec2.parentMass)/2);
    int upperShiftLimit = shiftsOffset+middleShift, lowerShiftLimit=0;
    if(minPeakAreaOvlp>0) lowerShiftLimit = shiftsOffset - (int)round(10*(spec2.parentMass - max(minPeakAreaOvlp*spec1.parentMass,minPeakAreaOvlp*spec2.parentMass)));
    int middleTimesTwo = (int)round(10*(spec1.parentMass-spec2.parentMass)); */
//cerr<<"middleShift="<<middleShift<<", upperShiftLimit="<<upperShiftLimit<<", lowerShiftLimit="<<lowerShiftLimit<<", middleTimesTwo="<<middleTimesTwo<<endl;
    int middleShift = (int)round((spec1.parentMass-spec2.parentMass)/(2*InputParams::Resolution));
    int upperShiftLimit = shiftsOffset+middleShift, lowerShiftLimit=0;
    if(minPeakAreaOvlp>0) lowerShiftLimit = shiftsOffset - (int)round((spec2.parentMass - max(minPeakAreaOvlp*spec1.parentMass,minPeakAreaOvlp*spec2.parentMass))/InputParams::Resolution);
    int middleTimesTwo = (int)round((spec1.parentMass-spec2.parentMass)/InputParams::Resolution);
//cerr<<"middleShift2="<<middleShift2<<", upperShiftLimit2="<<upperShiftLimit2<<", lowerShiftLimit2="<<lowerShiftLimit2<<", middleTimesTwo2="<<middleTimesTwo2<<endl;
    int shiftSym, validCount=0;
    list<TwoValues<float> > shiftsTmp;   TwoValues<float> curShift;
	for(shiftIndex=lowerShiftLimit; shiftIndex<upperShiftLimit; shiftIndex++) {
		if(!validShift[shiftIndex]) continue; validCount++;
		shiftSym = shiftsOffset+(middleTimesTwo-(shiftIndex-shiftsOffset));

		// Choose highest scoring symmetric shift within tolerance
		float maxScore=0, score1, score2, curPenalty; int maxScoreIdx=-1;
		for(int symIdx=max(0,shiftSym-intPeakTol); symIdx<=shiftSym+intPeakTol and symIdx<shiftScoresV.size() ; symIdx++) {
			score1 = shiftScoresV[shiftIndex][0]+shiftScoresV[shiftSym][0];  // +shiftScoresV[symIdx][0];  Bug fixed 2006/06/19
			score2 = shiftScoresV[shiftIndex][1]+shiftScoresV[shiftSym][1];  // +shiftScoresV[symIdx][1];  // 
			if(validShift[symIdx] and (score1+score2)>(bestCandidateScores[0]+bestCandidateScores[1])) bestCandidateScores.set(score1,score2);
			if(score1<minScore1 or score2<=minScore2 or !validShift[symIdx]) continue;
			else { curPenalty=round(abs(shiftSym-(int)symIdx)*InputParams::Resolution); if(score1+score2-curPenalty>maxScore) { maxScore=score1+score2-curPenalty; maxScoreIdx=symIdx; } } 
		}
		if(maxScoreIdx==-1) continue;

		curShift.set(maxScore,shiftIndex);   shiftsTmp.push_back(curShift);
		shiftPairsV[shiftIndex] = maxScoreIdx;
	}
//cerr<<"["<<validCount<<"/"<<shiftsTmp.size()<<"]";cerr.flush();
//	int maxUsedShiftIndex = aux3(intPeakTol, shiftMatchedPeaks, shiftsTmp, shiftScores, shiftPairs);

	// Sort shiftsTmp, create output structures
	shiftsTmp.sort();
	int maxUsedShiftIndex=0;  for(unsigned int i=0;i<usedShift.size();i++) usedShift[i]=false;
	list<TwoValues<float> >::reverse_iterator iter = shiftsTmp.rbegin();
	for(; iter!=shiftsTmp.rend(); iter++) {
		shiftIndex = (int)round((*iter)[1]);   shiftSym = shiftPairsV[shiftIndex];   if(maxUsedShiftIndex<shiftSym) maxUsedShiftIndex=shiftSym;
		shiftScores.push_back((*iter)[0]);
		shiftPairs.push_back(TwoValues<unsigned int>(shiftIndex,shiftSym));
		usedShift[shiftIndex]=true;   for(int shiftSymTol=max(0,shiftSym-intPeakTol); shiftSymTol<=shiftSym+intPeakTol and shiftSymTol<usedShift.size(); shiftSymTol++) usedShift[shiftSymTol]=true;
	}
	
	// Clear list of matched peaks for unused shifts
	for(unsigned int i=0; i<shiftMatchedPeaks.size(); i++)
		if(!usedShift[i]) shiftMatchedPeaks[i].resize(0);
//		if(!usedShift[i]) shiftMatchedPeaks[i].clear();

	return TwoValues<int>(shiftsOffset,maxUsedShiftIndex);
}

// Helper global variables for computeShifts2
static bool computeShifts2firstUse=true;
void computeShifts2(Spectrum &spec1, Spectrum &spec2, vector<float> &shiftsList) {
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
		if(fabs(shiftsList[shiftsIdx]-shiftsList[shiftsIdxUnique])>0.01) {
			shiftsIdxUnique++;
			if(shiftsIdx>shiftsIdxUnique) shiftsList[shiftsIdxUnique]=shiftsList[shiftsIdx];
		}
	shiftsList.resize(shiftsIdxUnique+1);
}

//
//  ScoreOverlap6 - like ScoreOverlapE.m with distance function 6. Min distance between PRMs is 57 - 2*tolerance
//
//  idx1all, idx2all - as returned by FindMatchPeaksAll
//  idxMatched       - Indices of the matched PRMs (sparse set), col 0 for spec1 and col 1 for spec2
//
float ScoreOverlap6(Spectrum &spec1, vector<short> idx1all, Spectrum &spec2, vector<short> idx2all, float shift, 
					 float tolerance, vector<short> &idxMatched1, vector<short> &idxMatched2, bool matchEndpoints, float *offsetPenalty){
    
    vector< TwoValues<float> > values(idx1all.size()+1);  // Keeps the values for the dynamic programming recursion: predecessor (col 0) and predecessor score (col 1)
                                                          // +1 because first line is (0,0) for DP initialization
    short forbiddenIdx;  // Index of the first forbidden PRM (because it is too close to the current PRM
    short maxIdx;        // Index of the best PRM that the current PRM can connect to
    float bestMatch;     // Best path score so far  
    short bestMatchIdx;  // Index of the last PRM in the best path so far
    int   i,j;           // Iterator vars
    
    idxMatched1.resize(0);   idxMatched2.resize(0);
    
    values[0][0]=0;   values[0][1]=0;
    maxIdx = 0;       forbiddenIdx = 1;
    bestMatch = 0;    bestMatchIdx = 0;
    for (i=1; i<(int)values.size(); i++) {
        while (spec1[idx1all[forbiddenIdx-1]][0] <= (spec1[idx1all[i-1]][0]-57+2*tolerance) &&
               spec2[idx2all[forbiddenIdx-1]][0] <= (spec2[idx2all[i-1]][0]-57+2*tolerance)) {
            // This is executed only when forbidden is advanced
            if (values[forbiddenIdx][1] > values[maxIdx][1]) { maxIdx = forbiddenIdx; }
            forbiddenIdx++;
        }
        values[i][0] = maxIdx;
        values[i][1] = values[maxIdx][1] + spec1[idx1all[i-1]][1] + spec2[idx2all[i-1]][1];
        if(offsetPenalty) values[i][1]-= abs(spec1[idx1all[i-1]][0] - (spec2[idx2all[i-1]][0]+shift));        	
        if (values[i][1]>bestMatch) { bestMatch=values[i][1];  bestMatchIdx=i; }  // Keep track of where the best path ends
    }

    list<short> bestPath;
    while (bestMatchIdx>0) { 
    	if(offsetPenalty) (*offsetPenalty)+= abs(spec1[idx1all[bestMatchIdx-1]][0] - (spec2[idx2all[bestMatchIdx-1]][0]+shift));
    	bestPath.push_back(bestMatchIdx-1);  bestMatchIdx=(short)values[bestMatchIdx][0]; 
    }
    
    // *********************************************
    //    Find the scores of the endpoint matches
    // *********************************************
    //
    // Option 1 - Spectra have no additional strange peaks at 0,19,PM-19,PM Da and these matches are checked here
    //            Problem: how to incorporate with sparse subset DP ?  Not likely to be a problem if spectra have no peaks <57 nor >PM+19-57
    // Option 2 - Spectra have the additional strange peaks
    //            Problem: Peaks can be matched incoherently, e.g. by matching both (0 and PM+19) or (19 and PM)
    
    // Option 1 implemented. Evaluate all 4 possible PRM matches in a single pass
    float toMatch[4][4];  // 4 points to be matched (0,19,PM-19,PM), each has mass in spec 1 (col 1), mass in spec 2 (col 2), match score (col 3), match PRM idx (col 4)
	float matchedSum = 0;   // Keeps track of whether any PRMs were matched at all (to endpoint PRMs)
	if(!matchEndpoints) { for(int i=0;i<4;i++) for(int j=0;j<4;j++) toMatch[i][j]=0; }
	else {
	    if (shift>=57-tolerance) {
	        toMatch[0][0] = shift;       toMatch[0][1] = 0;
	        toMatch[1][0] = shift+19;    toMatch[1][1] = 19;
	    } else if (shift<=57+tolerance) {
	        toMatch[0][0] = 0;       toMatch[0][1] = -shift;
	        toMatch[1][0] = 19;      toMatch[1][1] = -shift+19;
	    } else { toMatch[0][0]=0; toMatch[0][1]=0; toMatch[1][0]=0; toMatch[1][1]=0; }  // Starts match, don't look for matching PRMs
	    
	    if (spec1.parentMass>shift+spec2.parentMass+57-2*tolerance) {   // 2*tolerance to account for increased parent mass error tolerance
	        toMatch[2][0] = shift+spec2.parentMass-19;   toMatch[2][1] = spec2.parentMass-19;
	        toMatch[3][0] = shift+spec2.parentMass;      toMatch[3][1] = spec2.parentMass;
	    } else if (spec1.parentMass<shift+spec2.parentMass-57+2*tolerance) {
	        toMatch[2][0] = spec1.parentMass-19;   toMatch[2][1] = spec1.parentMass-shift-19;
	        toMatch[3][0] = spec1.parentMass;      toMatch[3][1] = spec1.parentMass-shift;
	    } else { toMatch[2][0]=0; toMatch[2][1]=0; toMatch[3][0]=0; toMatch[3][1]=0; }  // Ends match, don't look for matching PRMs
	    
	    i=0; j=0;               // initialize the iterators for spec1 (i) and spec2 (j)
	    for (int curMatch=0; curMatch<4; curMatch++) {
	        toMatch[curMatch][2]=0;   toMatch[curMatch][3]=0;
	        if (toMatch[curMatch][0]>=57-tolerance && toMatch[curMatch][0]<=spec1.parentMass-57+tolerance) {
	            while (i<spec1.size() && spec1[i][0]<toMatch[curMatch][0]-tolerance) i++;
	            while (i<spec1.size() && spec1[i][0]<=toMatch[curMatch][0]+tolerance) 
	                { if (spec1[i][1]>toMatch[curMatch][2]) { toMatch[curMatch][2]=spec1[i][1]; toMatch[curMatch][3]=i; } i++; }
	            matchedSum = matchedSum+toMatch[curMatch][3];
	            continue;
	        }
	        if (toMatch[curMatch][1]>=57-tolerance && toMatch[curMatch][1]<=spec2.parentMass-57+tolerance) {
	            while (j<spec2.size() && spec2[j][0]<toMatch[curMatch][1]-tolerance-0.0001) j++;
	            while (j<spec2.size() && spec2[j][0]<=toMatch[curMatch][1]+tolerance+0.0001) 
	                { if (spec2[j][1]>toMatch[curMatch][2]) { toMatch[curMatch][2]=spec2[j][1]; toMatch[curMatch][3]=j; } j++; }
	            matchedSum = matchedSum+toMatch[curMatch][3];
	        }
	    }
	    
	    if (matchedSum>0) {
	        // Choose the best matching pair: (0,PM-19) or (19,PM)
	        if (toMatch[0][2]+toMatch[2][2] >= toMatch[1][2]+toMatch[3][2]) { for(i=0;i<4;i++) toMatch[1][i]=toMatch[2][i]; }
	        else for(i=0;i<4;i++) { toMatch[0][i]=toMatch[1][i]; toMatch[1][i]=toMatch[3][i]; } // Copy both lines to maintain PRM order
	        
	        matchedSum = toMatch[0][2]+toMatch[1][2];
	
	        for(int curMatch=0; curMatch<2 ; curMatch++) {
	            if (toMatch[curMatch][0]>=57-tolerance && toMatch[curMatch][0]<=spec1.parentMass-57+tolerance) { toMatch[curMatch][0]=1; }  // insert in idxMatched1
	            if (toMatch[curMatch][1]>=57-tolerance && toMatch[curMatch][1]<=spec2.parentMass-57+tolerance) { toMatch[curMatch][0]=2; }  // insert in idxMatched2
	        }
	    }
	}
	
    // ******************************************
    //    Populate final lists of matched PRMs
    // ******************************************
    int bestPathLength = bestPath.size();
    idxMatched1.resize(bestPathLength+(toMatch[0][0]==1?1:0)+(toMatch[1][0]==1?1:0)+100);    
    idxMatched2.resize(bestPathLength+(toMatch[0][0]==2?1:0)+(toMatch[1][0]==2?1:0)+100);

    i=0,j=0; int k=0; // idxMatched iterators, k iterates through toMatch(0:1)
    for(int idxPath=0; idxPath<bestPathLength; idxPath++) {
        if(k<2 && toMatch[k][0]==1 && toMatch[k][3]<idx1all[bestPath.back()] && (i==0 || toMatch[k][3]>idxMatched1[i-1])) { idxMatched1[i++]=(short)toMatch[k][3]; k++; /* cerr<<"Added "<<idxMatched1[i-1]<<" (endpoint in spec1)\n"; */ }
        if(k<2 && toMatch[k][0]==2 && toMatch[k][3]<idx2all[bestPath.back()] && (j==0 || toMatch[k][3]>idxMatched2[j-1])) { idxMatched2[j++]=(short)toMatch[k][3]; k++; /* cerr<<"Added "<<idxMatched2[j-1]<<" (endpoint in spec2)\n"; */ }
        idxMatched1[i++] = idx1all[bestPath.back()];
        idxMatched2[j++] = idx2all[bestPath.back()];
// cerr<<"Added "<<idxMatched1[i-1]<<","<<idxMatched2[j-1]<<endl; 
        bestPath.pop_back();
    }
    idxMatched1.resize(i);   idxMatched2.resize(j);
    
    return bestMatch + matchedSum;
}

//
//  ScoreOverlap7 - like ScoreOverlapE.m with distance function 7 but gets matching PRMs from findMatchPeaksAll instead of findMatchPeaks2a. Min distance between PRMs is 57 - 2*tolerance
//
//  idx1all, idx2all - as returned by FindMatchPeaksAll
//  idxMatched1,2    - Indices of the matched PRMs (sparse set) for spec1 and for spec2
//
float ScoreOverlap7(Spectrum &spec1, vector<short> idx1all, Spectrum &spec2, vector<short> idx2all, float shift, float tolerance, 
                    vector<short> &idxMatched1, vector<short> &idxMatched2){
	Spectrum tmpSpec;   vector<short> idxMatched;   float score;
	
// Use FindMatchPeaksAll results to construct a single spectrum with all possible peak matches
	tmpSpec.parentMass = (spec1.parentMass+spec2.parentMass)/2;
	tmpSpec.peakList.resize(idx1all.size());
	tmpSpec.idDist = spec1.idDist;

	for(unsigned int i=0; i<tmpSpec.size(); i++) 
		tmpSpec[i].set((spec1[idx1all[i]][0]+spec2[idx2all[i]][0]+shift)/2,spec1[idx1all[i]][1]+spec2[idx2all[i]][1]);

// Use getMaxSparseSet to select which peak matches to keep and determine match score
	score = getMaxSparseSet(tmpSpec, tolerance, 0, idxMatched);
	
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
//
//  NOTE: Make sure that spec.idDist is adequately set.
//
float getMaxSparseSet(Spectrum &spec, float tolerance, float pmOffset, vector<short> &idxMatched){
	vector<vector<TwoValues<short> > > prevPair,     // Previous pair for every possible pair
	                                   bestPair;     // Best pair so far
	vector<vector<float> > scores,      // Scores for every possible path
	                       bestScores;  // Best path score so far
	vector<float> peakScores;    // Adjusted peak scores (including the score of the symmetric peaks)
	vector<short> peakPairs;     // Indices of the other paired peak (if any) used in peakScores
	vector<TwoValues<short> > prefixPRMs;  // PRMs with mass<=aaMass/2 and sorted by increasing mass
	                                       //  Col 1 is PRM index, col 2 is index of closest PRM >=57-2*tolerance Da away
	vector<TwoValues<short> > suffixPRMs;  // PRMs with mass>aaMass/2 and sorted by increasing distance to aaMass
	TwoValues<short> globalBestPair(0,0);
	float aaMass = spec.parentMass+pmOffset-spec.idDist, 
		  minInterPeakDist = 57*spec.idDist-2*tolerance,
		  globalBestScore=0;
	short i,j,k,p,idxPref,idxSuff;

	if(spec.size()==0) { idxMatched.resize(0); return 0; }

	// Make sure that the spectrum has PRMs at zero and at aaMass
	short addZero=0, addPM=0;  
	Spectrum oldSpec;           // oldSpec is used to keep a copy of the input spectrum
	                            //   whenever addZero or addPM are >0
	if(spec[0][0]>tolerance) addZero++;
	if(spec[spec.size()-1][0]<aaMass-18*spec.idDist-tolerance) addPM++;
	if(addZero+addPM>0) {
		oldSpec = spec;
		spec.peakList.resize(spec.size()+addZero+addPM);
		if (addZero) { for(i=spec.size()-1-addPM; i>=1; i--) spec[i]=spec[i-1]; spec[0].set(0,0); }
		if (addPM) spec[spec.size()-1].set(aaMass-18*spec.idDist,0);
	}

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
	vector<vector<float> > pairs; vector<vector<short> > pairsIdx;
	spec.getPairs(pmOffset, tolerance, pairs, pairsIdx);
	for(k=0;k<pairsIdx.size();k++) {
		if(pairsIdx[k][0]<0 || pairsIdx[k][1]<0) continue;
		if(peakScores[pairsIdx[k][0]]<pairs[k][2]) { peakScores[pairsIdx[k][0]]=pairs[k][2]; peakPairs[pairsIdx[k][0]]=pairsIdx[k][1]; }
		if(peakScores[pairsIdx[k][1]]<pairs[k][2]) { peakScores[pairsIdx[k][1]]=pairs[k][2]; peakPairs[pairsIdx[k][1]]=pairsIdx[k][0]; }
	}

//cerr << "Prefixes:\n";
//for(k=0;k<prefixPRMs.size();k++) cerr<<k<<": ["<<prefixPRMs[k][0]<<","<<prefixPRMs[k][1]<<"], mass = "<<spec[prefixPRMs[k][0]][0]<<", score = "<<peakScores[prefixPRMs[k][0]]<<"\n";
//cerr << "Suffixes:\n";
//for(k=0;k<suffixPRMs.size();k++) cerr<<k<<": ["<<suffixPRMs[k][0]<<","<<suffixPRMs[k][1]<<"], mass = "<<spec[suffixPRMs[k][0]][0]<<", score = "<<peakScores[suffixPRMs[k][0]]<<"\n";

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
	TwoValues<short> tmpPair;     short curPeakIdx;
	vector<bool> idxMatchedBool;   // Boolean vector used to mark matched peaks
	idxMatchedBool.resize(spec.size());   for(int i=0; i<spec.size(); i++) idxMatchedBool[i]=false;
	if(globalBestPair[0]>0) { 
		idxMatchedBool[prefixPRMs[globalBestPair[0]][0]-addZero]=true;
		if(peakPairs[prefixPRMs[globalBestPair[0]][0]]>0) idxMatchedBool[peakPairs[prefixPRMs[globalBestPair[0]][0]]-addZero]=true;
	}
	if(globalBestPair[1]>0) {
		idxMatchedBool[suffixPRMs[globalBestPair[1]][0]-addZero]=true;
		if(peakPairs[suffixPRMs[globalBestPair[1]][0]]>0) idxMatchedBool[peakPairs[suffixPRMs[globalBestPair[1]][0]]-addZero]=true;
	}
	while (globalBestPair[0]>0 || globalBestPair[1]>0) {
		tmpPair = prevPair[globalBestPair[0]][globalBestPair[1]];
		curPeakIdx=-1;
		if(tmpPair[0]!=globalBestPair[0]) 
			{ if(tmpPair[0]>0) curPeakIdx=prefixPRMs[tmpPair[0]][0]; }
		else if(tmpPair[1]>0) curPeakIdx=suffixPRMs[tmpPair[1]][0];
		if (curPeakIdx>0) {
			idxMatchedBool[curPeakIdx-addZero]=true;
			if (peakPairs[curPeakIdx]>0) idxMatchedBool[peakPairs[curPeakIdx]-addZero]=true;
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
