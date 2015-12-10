#include "batch.h"
#include "aminoacid.h"
#include "denovo.h"
#include "alignment_scoring.h"
#include "inputParams.h"

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdio>

pthread_mutex_t cout_mutex = PTHREAD_MUTEX_INITIALIZER;

//
// specSet        - Set of spectra to align
// baseSpectraIdx - Indices of the spectra to align to others on this run. Note that spectra only align to spectra with higher indices.
// aaDiff         - Difference between matched spectra's parent masses must be at most the sum of aaDiffs amino acid masses
//                    or -1 for computation of all shifts between minShift and maxShift
// minShift, maxShift - Min and max parent mass differences between spectra (to compute ASP shifts)
// pmTol          - Parent mass error tolerance
// peakTol        - Peak mass error tolerance
// minMatchRatio  - minimum ratio of matched score / max score (in both spectra) to retain a match
// results        - Variable to hold results
// resolution     - Resolution for parent mass histogram (to determine which spectra to match)
//
void getPairAlignsASP(SpecSet &specSet, vector<int> &baseSpectraIdx, short aaDiff, float minShift, float maxShift,
									float pmTol, float peakTol, float minMatchRatio, list<Results_ASP> &results,
									list<TwoValues<float> > &ratios, vector<TwoValues<float> > &means, vector<float> &varTerms,
									vector<vector<TTag> > &tags, unsigned int tagsMatchFlank, unsigned int tagsMatchCount,
									float resolution, float symmetryOffset){
	int j, massIdx;
	float pMass1;
	float shiftIncrement = 2*resolution;  // Step size used to increment shifts from -pmTol to pmTol
	list<int>::iterator curPair;
	float histResolution = 0.1;           // Resolution of the parent mass histogram - only used for binning, independent of overall peak resolution settings
	vector<list<int> > hist;   specSet.getMassesHistogram(hist,histResolution);
	TwoValues<float> curRatios;
	AAJumps jumps(-1);
	Results_ASP curResult;

	// Determine parent mass differences for valid matches
	if (aaDiff>=0) {
		jumps.getjumps(aaDiff,resolution);
		jumps.forceJump(0);
		jumps.forceTolerance(pmTol,resolution);
		jumps.forceDoubleSided();
		jumps.forceUnique(resolution);
	}

	// Determine totalScores per spectrum for later computation of match ratios
	vector<float> totalSpecScores;   totalSpecScores.resize(specSet.size());
	for(unsigned int i=0; i<specSet.size(); i++) { totalSpecScores[i]=0;
		for(unsigned int p=0; p<specSet[i].size(); p++) totalSpecScores[i]+=specSet[i][p][1];
	}

	// Initialize means vector
	means.resize(specSet.size());   varTerms.resize(specSet.size());
	for(unsigned int i=0; i<means.size(); i++) { means[i].set(0,0); varTerms[i]=0; }

	// Iterate through pairs using FindMatchPeaksAll and ScoreOverlap6
//    vector<int> idx1,idx2;                             idx1.reserve(1500);               idx2.reserve(1500);
	vector<int> idxMatched1_zero, idxMatched1_other;   idxMatched1_zero.reserve(1500);   idxMatched1_other.reserve(1500);
	vector<int> idxMatched2_zero, idxMatched2_other;   idxMatched2_zero.reserve(1500);   idxMatched2_other.reserve(1500);
	vector<int> idxUnion;                              idxUnion.reserve(1500);
	vector<int>::iterator idxUnionEnd, idxUnionStart;
	int spec1, spec2,
		maxShiftInt = (int)ceil(maxShift/histResolution);   minShift = max((float)0,minShift);
	float pmDiff, score1, score2, shift, bestScore1, bestScore2, bestShift;
//    clock_t startTime = clock();    double curTime, totTime=0;
	list<int> candidates;
	list<int>::iterator curCandidate;
	for(unsigned int i=0; i<baseSpectraIdx.size(); i++) {
		cout << "Starting with spectrum " << baseSpectraIdx[i] << "..."; cout.flush();

		// Find candidate spectra to compare to
		spec1 = baseSpectraIdx[i];
		pMass1 = specSet[spec1].parentMass;   int massInt = (int)round(pMass1/histResolution);
		candidates.clear();  // Make sure candidates is empty from previous run
	if(specSet[spec1].size()==0 or specSet[spec1].parentMass<400) { cout<<"(spectrum is empty or parent mass too low)\n"; continue; }
		// Find candidates in [minShift,maxShift]
		int upperLimit = max(0,(int)ceil((pMass1-minShift)/histResolution));  // Look at interval below current parent mass
		for(massIdx=(int)floor((pMass1-min(pMass1,maxShift))/histResolution); massIdx<upperLimit; massIdx++){
			for(curPair = hist[massIdx].begin(); curPair!=hist[massIdx].end(); curPair++) {
				if (*curPair>spec1 and fabs(pMass1-specSet[*curPair].parentMass)>=minShift-2*pmTol and fabs(pMass1-specSet[*curPair].parentMass)<=maxShift+2*pmTol)
					candidates.push_back(*curPair);
			}
		}
		upperLimit = min((int)hist.size(),massInt+maxShiftInt);      // Look at interval above current parent mass
		for(massIdx=(int)floor((pMass1+minShift)/histResolution); massIdx<upperLimit; massIdx++){
			for(curPair = hist[massIdx].begin(); curPair!=hist[massIdx].end(); curPair++) {
				if (*curPair>spec1 and fabs(pMass1-specSet[*curPair].parentMass)>=minShift-2*pmTol and fabs(pMass1-specSet[*curPair].parentMass)<=maxShift+2*pmTol)
					candidates.push_back(*curPair);
			}
		}

		if(aaDiff>0) {
			for(j=0, massIdx=(int)round((pMass1+jumps[j])/histResolution); j<jumps.size(); j++) if((int)round((pMass1+jumps[j])/histResolution)>=0) break;  // Skips masses<0
			for(; j<jumps.size(); j++){   if(fabs(jumps[j])>=minShift-2*pmTol and fabs(jumps[j])<=maxShift+2*pmTol) continue;  // Already added above
				massIdx=(int)round((pMass1+jumps[j])/histResolution); if(massIdx>=hist.size()) break;
				for(curPair = hist[massIdx].begin(); curPair!=hist[massIdx].end(); curPair++)
					if (*curPair>spec1 and fabs(pMass1+jumps[j]-specSet[*curPair].parentMass)<=2*pmTol) { candidates.push_back(*curPair); }
			}
		}
		candidates.sort();   candidates.unique();   // The same spectrum may occur in multiple bins - make sure the code below only aligns once to each eligible spectrum
		cout << " num. candidates is " << candidates.size(); cout.flush();

		if(tags.size()>0) {
			vector<unsigned int> candV(candidates.size());   unsigned int cIdx=0;
			for(curCandidate=candidates.begin(); curCandidate!=candidates.end(); curCandidate++, cIdx++) candV[cIdx]=(unsigned int)(*curCandidate);
			IntersectTags(tags,spec1,candV,tagsMatchCount,tagsMatchFlank);
			candidates.clear(); for(cIdx=0; cIdx<candV.size(); cIdx++) candidates.push_back((int)candV[cIdx]);
			cout<<" ("<<candidates.size()<<" after tag filtering)";
		}
		cout << "..."; cout.flush();

		curCandidate = candidates.begin();
		while(curCandidate!=candidates.end()) {
			spec2 = *curCandidate;   curCandidate++;
			if (spec1>=specSet.size() || spec2>=specSet.size()) { cerr << "ERROR: Alignment pair ("<<spec1<<","<<spec2<<") is invalid - PRM spectra set has only "<<specSet.size()<<" spectra!\n"; exit(-1); }
//cerr<<spec2<<":"; cerr.flush();
			pmDiff = specSet[spec1].parentMass-specSet[spec2].parentMass;
			bestScore1 = -99999; bestScore2 = -99999;
			if (abs(pmDiff)>pmTol+0.00001) {
//	            FindMatchPeaksAll(specSet[spec1], specSet[spec2], 0, peakTol, idx1, idx2);
//cerr<<"X"; cerr.flush();
				ScoreOverlap6(specSet[spec1], specSet[spec2], 0, peakTol, idxMatched1_zero, idxMatched2_zero);
//cerr<<"X"; cerr.flush();

				for(shift = pmDiff-pmTol; shift<=pmDiff+pmTol+0.00001; shift+=shiftIncrement) {
//	                FindMatchPeaksAll(specSet[spec1], specSet[spec2], shift, peakTol, idx1, idx2);
//cerr<<"<"; cerr.flush();
					ScoreOverlap6(specSet[spec1], specSet[spec2], shift, peakTol, idxMatched1_other, idxMatched2_other);
//cerr<<">"; cerr.flush();

					idxUnionEnd = set_union(idxMatched1_zero.begin(),idxMatched1_zero.end(),idxMatched1_other.begin(),idxMatched1_other.end(),idxUnion.begin());
					for(idxUnionStart = idxUnion.begin(), score1=0; idxUnionStart!=idxUnionEnd; idxUnionStart++) score1+=specSet[spec1][(*idxUnionStart)][1];

					idxUnionEnd = set_union(idxMatched2_zero.begin(),idxMatched2_zero.end(),idxMatched2_other.begin(),idxMatched2_other.end(),idxUnion.begin());
					for(idxUnionStart = idxUnion.begin(), score2=0; idxUnionStart!=idxUnionEnd; idxUnionStart++) score2+=specSet[spec2][(*idxUnionStart)][1];

					if ((score1+score2>bestScore1+bestScore2) ||
						(score1+score2==bestScore1+bestScore2 && abs(shift-pmDiff)<abs(bestShift-pmDiff)))
						{ bestScore1=score1; bestScore2=score2; bestShift=shift; }
				}
			} else {
				vector<int> idx1all, idx2all, idxMatched1, idxMatched2;

				for(shift=-pmTol; shift<=pmTol+0.00001; shift+=shiftIncrement) {
//cerr<<"{"; cerr.flush();
					FindMatchPeaksAll(specSet[spec1], specSet[spec2], shift, peakTol, idx1all, idx2all);
//cerr<<"}"; cerr.flush();
					ScoreOverlap7(specSet[spec1], idx1all, specSet[spec2], idx2all, shift, peakTol, idxMatched1, idxMatched2, symmetryOffset);

					score1=0; score2=0;
					for(int i=0; i<idxMatched1.size(); i++) score1+=specSet[spec1][idxMatched1[i]][1];
					for(int i=0; i<idxMatched2.size(); i++) score2+=specSet[spec2][idxMatched2[i]][1];
					if ((score1+score2>bestScore1+bestScore2) ||
						(score1+score2==bestScore1+bestScore2 && abs(shift)<abs(bestShift)))
						{ bestScore1=score1; bestScore2=score2; bestShift=shift;
/*
cerr<<"Matching peaks with shift "<<bestShift<<", scores ("<<score1<<","<<score2<<"):\nSpectrum "<<spec1<<":\n";
for(unsigned int i=0; i<idxMatched1.size(); i++) cerr<<specSet[spec1][idxMatched1[i]][0]<<"\t"<<specSet[spec1][idxMatched1[i]][1]<<endl;
cerr<<"Spectrum "<<spec2<<":\n";
for(unsigned int i=0; i<idxMatched2.size(); i++) cerr<<specSet[spec2][idxMatched2[i]][0]<<"\t"<<specSet[spec2][idxMatched2[i]][1]<<endl;
*/
						}
				}
			}
//cerr<<"|"; cerr.flush();

			curResult.spec1 = spec1;   curResult.spec2 = spec2;
			curResult.score1 = bestScore1;   curResult.score2 = bestScore2;   curResult.shift1 = bestShift;

			means[spec1][0]=(means[spec1][0]*means[spec1][1]+curResult.score1)/(means[spec1][1]+1);
			means[spec2][0]=(means[spec2][0]*means[spec2][1]+curResult.score2)/(means[spec2][1]+1);
			means[spec1][1]++;    means[spec2][1]++;
			float updateRatio = (means[spec1][1]-1)/means[spec1][1];  // Update variance terms
			varTerms[spec1] = updateRatio*varTerms[spec1]+(curResult.score1*curResult.score1)/means[spec1][1];
			updateRatio = (means[spec2][1]-1)/means[spec2][1];
			varTerms[spec2] = updateRatio*varTerms[spec2]+(curResult.score2*curResult.score2)/means[spec2][1];

			// Compute match ratios
			curRatios[0] = curResult.score1/totalSpecScores[spec1];
			curRatios[1] = curResult.score2/totalSpecScores[spec2];
			if(curRatios[0]>=minMatchRatio && curRatios[1]>=minMatchRatio)
				{ results.push_back(curResult); ratios.push_back(curRatios); }

//cerr<<"|"; cerr.flush();
		}  // while over all candidates for spec1
		cout << " done. Current number of results in memory: "<<results.size()<<"\n";
	} // for i over baseSpectraIdx
	cout << "Done with the alignment!\n";
}


void getPairAlignsPA(SpecSet &specSet, float pmTol, float peakTol, vector<Results_PA> &results){
	int pair, massOffset, baseShift, shift, shiftSym, middleShift, middleTimesTwo, spec1, spec2, i, j, k;
	int intPeakTolMin=-(int)round(peakTol*10*2), intPeakTolMax=(int)round(10*peakTol*2);
	int intPmTolMin=-(int)round(pmTol*10), intPmTolMax=(int)round(10*pmTol);
	float realShift1, realShift2, maxMass;
	vector<bool> shifts;
//    vector<int> idx1,idx2;                        idx1.reserve(500);              idx2.reserve(500);
	vector<int> idxMatched1,idxMatched2;          idxMatched1.reserve(500);       idxMatched2.reserve(500);
	vector<int> idxMatchedSym1,idxMatchedSym2;    idxMatchedSym1.reserve(500);    idxMatchedSym2.reserve(500);
	vector<int> idxUnion;                         idxUnion.reserve(500);
	vector<int>::iterator idxUnionEnd, idxUnionStart;

cerr<<"WARNING: getPairAlignsPA() in batch.cpp is not ready for high-accuracy spectra!\n";

	if (results.size()==0) {
		// Fill results.spec1 and results.spec2 to compute ALL pairwise alignments
		results.resize((specSet.size()*(specSet.size()-1))/2);   i=0;
		for(spec1=0; spec1<specSet.size(); spec1++)
			for(spec2=spec1+1; spec2<specSet.size(); spec2++)
				{ results[i].spec1=spec1;   results[i].spec2=spec2;   i++; }
	}

	// Pre-allocate vector to keep track of valid shifts - implicit 0.1 Da resolution for shifts
	for(i=0, maxMass=0; i<=specSet.size(); i++) if (specSet[i].parentMass>maxMass) maxMass=specSet[i].parentMass;
	massOffset = (int)ceil(maxMass*10);
	shifts.resize(2*massOffset+1);

	for(pair=0; pair<=results.size(); pair++) {
		spec1 = results[pair].spec1;   spec2 = results[pair].spec2;

		// Step 1: Get the list of possible shifts
		for(shift=0; shift<shifts.size(); shift++) shifts[shift]=0;
		for(i=0; i<specSet[spec1].size(); i++)
			for(j=0; j<specSet[spec2].size(); j++) {
				baseShift = massOffset+(int)round(10*(specSet[spec1][i][0]-specSet[spec2][j][0]));
				for(k=intPeakTolMin; k<intPeakTolMax; k++) shifts[baseShift+k]=1;
			}


// THERE IS A LOT OF RECOMPUTATION GOING ON HERE!
// Every symmetric shift is being computed as many times as it is used - once should be enough if the memory requirements stay acceptable

		// Step 2: Compute the summed scores for all pairs of symmetric shifts (use pmTol) + choose the best pair
		float score1, bestScore1, score2, bestScore2, bestShift1, bestShift2;
		middleShift = massOffset+(int)floor(5*(specSet[spec1].parentMass-specSet[spec2].parentMass));
		middleTimesTwo = (int)round(10*(specSet[spec1].parentMass-specSet[spec2].parentMass));

cerr << "Spec1: " << spec1 << ", parent mass = " << specSet[spec1].parentMass << endl;
cerr << "Spec2: " << spec2 << ", parent mass = " << specSet[spec2].parentMass << endl;
cerr << "middleShift = " << middleShift-massOffset << ", middleTimeTwo = " << middleTimesTwo << endl;
int first=0;

		for(shift=0; shift<middleShift; shift++) {
			if (shifts[shift]==0) continue;
			realShift1=(shift-massOffset)/10;

//            FindMatchPeaksAll(specSet[spec1], specSet[spec2], realShift1, peakTol, idx1, idx2);
			ScoreOverlap6(specSet[spec1], specSet[spec2], realShift1, peakTol, idxMatched1, idxMatched2);

			shiftSym = massOffset+(middleTimesTwo-(shift-massOffset));
			bestScore1=0; bestScore2=0; bestShift1=0; bestShift2=0;

if (first==0) {
	cerr << "shift = " << massOffset-shift << ", realShift1 = " << realShift1 << endl;
}

			for(k=intPmTolMin; k<intPmTolMax; k++){
				if (shifts[shiftSym+k]==0) continue;
				realShift2=(shiftSym-massOffset)/10;

//                FindMatchPeaksAll(specSet[spec1], specSet[spec2], realShift2, peakTol, idx1, idx2);
				ScoreOverlap6(specSet[spec1], specSet[spec2], realShift2, peakTol, idxMatchedSym1, idxMatchedSym2);

				// Union of matched peaks
				idxUnionEnd = set_union(idxMatched1.begin(),idxMatched1.end(),idxMatchedSym1.begin(),idxMatchedSym1.end(),idxUnion.begin());
				for(idxUnionStart = idxUnion.begin(), score1=0; idxUnionStart!=idxUnionEnd; idxUnionStart++) score1+=specSet[spec1][(*idxUnionStart)][1];

				idxUnionEnd = set_union(idxMatched2.begin(),idxMatched2.end(),idxMatchedSym2.begin(),idxMatchedSym2.end(),idxUnion.begin());
				for(idxUnionStart = idxUnion.begin(), score2=0; idxUnionStart!=idxUnionEnd; idxUnionStart++) score2+=specSet[spec2][(*idxUnionStart)][1];

if (first==0) {
	cerr << "-- shiftSym = " << massOffset-shiftSym << ", realShift2 = " << realShift2 << ", score1 = " << score1 << ", score2 = " << score2 << endl;
}

				if ((score1+score2>bestScore1+bestScore2) ||
					(score1+score2==bestScore1+bestScore2 && abs(shift+shiftSym-middleTimesTwo-2*massOffset)<abs(bestShift1+bestShift2-middleTimesTwo-2*massOffset)))
					{ bestScore1=score1; bestScore2=score2; bestShift1=realShift1; bestShift2=realShift2;}
			}
if (first==0) {
	cerr << "bestShift1 = " << bestShift1 << ", bestShift2 = " << bestShift2 << ", bestScore1 = " << bestScore1 << ", bestScore2 = " << bestScore2 << endl;
	first = 1;
}
		}

cerr << "bestShift1 = " << bestShift1 << ", bestShift2 = " << bestShift2 << ", bestScore1 = " << bestScore1 << ", bestScore2 = " << bestScore2 << endl;

		results[pair].score1 = bestScore1;    results[pair].score2 = bestScore2;
		results[pair].shift1 = bestShift1;    results[pair].shift2 = bestShift2;
	}
}

void copy(MatepairContigAlignParams& from, MatepairContigAlignParams& to) {
	to.contigs = from.contigs;
	to.connectors = from.connectors;
	to.minRatio = from.minRatio;
	to.outputFile = from.outputFile;
	to.minNumMatchedPeaks = from.minNumMatchedPeaks;
	to.minAbsShift = from.minAbsShift;
	to.allowedJumps = from.allowedJumps;
	to.precCharge = from.precCharge;
	to.checkConsPeaks = from.checkConsPeaks;
}


bool alignMatepairContig(int mpIdx, int contigIdx, MatepairContigAlignParams* params, Results_OCC& result) {
	float minRatio = params->minRatio;
	int minNumMatchedPeaks = params->minNumMatchedPeaks;
	float minAbsShift = params->minAbsShift;
	AAJumps allowedJumps = *(params->allowedJumps);
	bool checkConsPeaks = params->checkConsPeaks;
	
	//cout << "minRatio: " << minRatio << ", minNumMatchedPeaks: " << minNumMatchedPeaks << ", minAbsShift: " << minAbsShift << ", checkConsPeaks: " << checkConsPeaks << endl;
	
	float peakTol = InputParams::PeakTol, pmTol = InputParams::PMTol;
	int intPMTol = (int)(round(pmTol/InputParams::Resolution)+0.01);
	int intPeakTol = (int)(round(peakTol/InputParams::Resolution)+0.01);
	int massH2OInt = (int)(round(AAJumps::massH2O/InputParams::Resolution)+0.01);
	
	Spectrum matepair = (*(params->connectors))[mpIdx];
	Spectrum contig = (*(params->contigs))[contigIdx];
	
	Results_OCC curResult;
	
	map<int, list<TwoValues<int> > > shiftMP;
	map<int, list<TwoValues<int> > > shiftMP_r;
	map<int, list<TwoValues<int> > >::iterator shiftMPIt;
	list<TwoValues<int> >::iterator mpIt;
	
	if (contig.size() == 0 || matepair.size() == 0) return false;
	
	vector<bool> o_contig(contig.size());
	for (int i = 0; i < o_contig.size(); i++) o_contig[i] = false;
	
	vector<bool> o_matepair(matepair.size());
	for (int i = 0; i < o_matepair.size(); i++) o_matepair[i] = false;
	
	Spectrum r_contig = contig;
	r_contig.reverse(0.0 - AAJumps::massH2O, 0);
	
	int contigPMInt = (int)(round(contig.parentMass/InputParams::Resolution)+0.01);
	int mpPMInt = (int)(round(matepair.parentMass/InputParams::Resolution)+0.01);
	
	computeShiftsRaw(matepair, contig, 0, shiftMP, 2, false);
	
	computeShiftsRaw(matepair, r_contig, 0, shiftMP_r, 2, false);
	
	map<int, int> seenShifts;

	float maxScore = 0.0;
	int peaks, o_peaks, c_peaks, cons_peaks, max_cons_peaks;
	float s1, s2, ratio;
	float o_m_intensity, o_t_intensity, o_a_intensity, c_m_intensity, c_t_intensity, c_a_intensity;
	float f_shift, r_shift, f_min, f_max, f_c_min, f_c_max, r_min, r_max, y_mass;
	float mass, intensity;
	
	curResult.spec1 = mpIdx;
	curResult.spec2 = contigIdx;
	
	for (shiftMPIt = shiftMP.begin(); shiftMPIt != shiftMP.end(); shiftMPIt ++) {
		int shiftIndex = shiftMPIt->first;
		int shiftSym = mpPMInt - contigPMInt - shiftIndex + massH2OInt;
		
		f_shift = ((float)shiftIndex)*InputParams::Resolution;
		f_min = f_shift - pmTol;
		f_max = f_shift + contig.parentMass + pmTol;
		f_c_min = 0.0 - f_shift - pmTol;
		f_c_max = matepair.parentMass - f_shift + pmTol;
		
		for (int symIdx = shiftSym-intPMTol; symIdx <= shiftSym+intPMTol; symIdx ++) {
			seenShifts[symIdx] = shiftIndex;
			float r_shift = ((float)symIdx)*InputParams::Resolution;
			
			r_min = r_shift - pmTol;
			r_max = r_shift + contig.parentMass + pmTol;

			o_m_intensity = 0.0; o_t_intensity = 0.0; o_a_intensity = 0.0; c_m_intensity = 0.0; c_t_intensity = 0.0; c_a_intensity = 0.0;
			o_peaks = 0; c_peaks = 0; cons_peaks = 0; max_cons_peaks = 0;
			
			for (mpIt = (shiftMPIt->second).begin(); mpIt != (shiftMPIt->second).end(); mpIt ++) {
				o_matepair[(*mpIt)[0]] = true;
				o_contig[(*mpIt)[1]] = true;
			}
			
			if (shiftMP_r.count(symIdx) > 0) {
				for (mpIt = shiftMP_r[symIdx].begin(); mpIt != shiftMP_r[symIdx].end(); mpIt ++) {
					o_matepair[(*mpIt)[0]] = true;
					o_contig[contig.size() - (*mpIt)[1] - 1] = true;
				}
			}
			
			f_max = max(f_max, r_max);
			f_min = min(f_min, r_min);
			
			for(int i = 0; i < matepair.size(); i++ ) {
				mass = matepair[i][0]; intensity = matepair[i][1];
				if (mass > f_min && mass < f_max) {o_a_intensity += intensity;}
				o_t_intensity += intensity;
				if (o_matepair[i]) {
					o_m_intensity += intensity;
					o_peaks ++;
				}
				o_matepair[i] = false;
			}
			
			for (int i = 0; i < contig.size(); i++) {
				mass = contig[i][0]; intensity = contig[i][1];
				
				if (mass >= f_c_min && mass <= f_c_max) {c_a_intensity += intensity;}
				
				c_t_intensity += intensity;
				if (o_contig[i]) {
					c_m_intensity += intensity;
					c_peaks ++;
					cons_peaks ++;
					if (cons_peaks > max_cons_peaks) max_cons_peaks = cons_peaks;
				} else cons_peaks = 0;
				
				o_contig[i] = false;
			}
			
			//s1 = o_m_intensity/o_t_intensity;
			//s2 = c_m_intensity/c_t_intensity;
			s1 = o_m_intensity/o_a_intensity;
			s2 = c_m_intensity/c_a_intensity;
			//s1 = o_m_intensity/((o_t_intensity/orbcal.parentMass) * (min(f_max, orbcal.parentMass) - max((float)0, f_min)));
			//s2 = c_m_intensity/((c_t_intensity/contig.parentMass) * (min(f_c_max, contig.parentMass) - max((float)0, f_c_min)));
			ratio = min(s1, s2) - 0.0001*(float)abs(shiftSym-symIdx);
			peaks = (checkConsPeaks) ? max_cons_peaks : min(o_peaks, c_peaks);
			
			if (ratio > maxScore && peaks >= minNumMatchedPeaks && ratio >= minRatio) {
				maxScore = ratio;
				
				curResult.score1 = s1;
				curResult.score2 = s2;
				curResult.shift1 = ((float)shiftIndex)*InputParams::Resolution;
				curResult.shift2 = ((float)symIdx)*InputParams::Resolution;
				curResult.matchedPeaks = peaks;
			}
		}
	}
	
	for (shiftMPIt = shiftMP_r.begin(); shiftMPIt != shiftMP_r.end(); shiftMPIt ++) {
		int shiftIndex = shiftMPIt->first;
		int shiftSym = mpPMInt - contigPMInt - shiftIndex + massH2OInt;
		
		float r_shift = ((float)shiftIndex)*InputParams::Resolution;
		r_min = r_shift - pmTol;
		r_max = r_shift + contig.parentMass + pmTol;
		
		
		for (int symIdx = shiftSym-intPMTol; symIdx <= shiftSym+intPMTol; symIdx ++) {
			if (seenShifts.count(symIdx) > 0 && seenShifts[symIdx] == shiftIndex) continue;
			float f_shift = ((float)symIdx)*InputParams::Resolution;
			
			f_min = f_shift - pmTol;
			f_max = f_shift + contig.parentMass + pmTol;
			f_c_min = 0.0 - f_shift - pmTol;
			f_c_max = matepair.parentMass - f_shift + pmTol;
			
			o_m_intensity = 0.0; o_t_intensity = 0.0; o_a_intensity = 0.0; c_m_intensity = 0.0; c_t_intensity = 0.0; c_a_intensity = 0.0;
			o_peaks = 0; c_peaks = 0; cons_peaks = 0; max_cons_peaks = 0;
			
			for (mpIt = (shiftMPIt->second).begin(); mpIt != (shiftMPIt->second).end(); mpIt ++) {
				o_matepair[(*mpIt)[0]] = true;
				o_contig[contig.size() - (*mpIt)[1] - 1] = true;
			}
			if (shiftMP.count(symIdx) > 0) {
				for (mpIt = shiftMP_r[symIdx].begin(); mpIt != shiftMP_r[symIdx].end(); mpIt ++) {
					o_matepair[(*mpIt)[0]] = true;
					o_contig[(*mpIt)[1]] = true;
				}
			}
			f_max = max(f_max, r_max);
			f_min = min(f_min, r_min);
			
			for(int i = 0; i < matepair.size(); i++ ) {
				mass = matepair[i][0]; intensity = matepair[i][1];
				if (mass > f_min && mass < f_max) {o_a_intensity += intensity;}
				o_t_intensity += intensity;
				if (o_matepair[i]) {
					o_m_intensity += intensity;
					o_peaks ++;
				}
			}
			for (int i = 0; i < contig.size(); i++) {
				mass = contig[i][0]; intensity = contig[i][1];
				
				if (mass >= f_c_min && mass <= f_c_max) {c_a_intensity += intensity;}
				
				c_t_intensity += intensity;
				if (o_contig[i]) {
					c_m_intensity += intensity;
					c_peaks ++;
					cons_peaks ++;
					if (cons_peaks > max_cons_peaks) max_cons_peaks = cons_peaks;
				} else cons_peaks = 0;
			}
			
			//s1 = o_m_intensity/o_t_intensity;
			//s2 = c_m_intensity/c_t_intensity;
			s1 = o_m_intensity/o_a_intensity;
			s2 = c_m_intensity/c_a_intensity;
			//s1 = o_m_intensity/((o_t_intensity/orbcal.parentMass) * (min(f_max, orbcal.parentMass) - max((float)0, f_min)));
			//s2 = c_m_intensity/((c_t_intensity/contig.parentMass) * (min(f_c_max, contig.parentMass) - max((float)0, f_c_min)));
			ratio = min(s1, s2) - 0.0001*(float)abs(shiftSym-symIdx);
			peaks = (checkConsPeaks) ? max_cons_peaks : min(o_peaks, c_peaks);
			
			if (ratio > maxScore && peaks >= minNumMatchedPeaks && ratio >= minRatio) {
				maxScore = ratio;
				
				curResult.score1 = s1;
				curResult.score2 = s2;
				curResult.shift1 = ((float)symIdx)*InputParams::Resolution;
				curResult.shift2 = ((float)shiftIndex)*InputParams::Resolution;
				curResult.matchedPeaks = peaks;
			}
		}
	}
	
	if (maxScore > 0.0) {
		result = curResult;
		return true;
	}
	return false;
}

void getMatepairContigAlignment(MatepairContigAlignParams* inparams, unsigned int num_threads) {
	inparams->results.clear();
	SpecSet* contigs = inparams->contigs;
	SpecSet* matepairs = inparams->connectors;
	short precCharge = inparams->precCharge;
	
	if (num_threads == 0) {
		cerr << "ERROR: num_threads must be >= 1\n";
		return;
	}
	
	float o_min_intensity = 0, c_min_intensity = 0;
	for (unsigned int orbit = 0; orbit < matepairs->size(); orbit ++) {
		for (int i = 0; i < (*matepairs)[orbit].size(); i ++) {
			if ((*matepairs)[orbit][i][1] < o_min_intensity) o_min_intensity = (*matepairs)[orbit][i][1];
		}
	}
	for (unsigned int contit = 0; contit < contigs->size(); contit ++) {
		if ((*contigs)[contit].size() == 0) continue;
		(*contigs)[contit].parentMass = (*contigs)[contit].peakList.back()[0] + AAJumps::massMH;
		for (int i = 0; i < (*contigs)[contit].size(); i ++) {
			if ((*contigs)[contit][i][1] < c_min_intensity) c_min_intensity = (*contigs)[contit][i][1];
		}
	}

	if (o_min_intensity < 0) {
		printf("Minimum matepair intensity: %.3f", o_min_intensity);
		printf(", correcting ... "); fflush(stdout);
		o_min_intensity = 0.0 - o_min_intensity;
		for (int orbit = 0; orbit < matepairs->size(); orbit ++) {
			for (int i = 0; i < (*matepairs)[orbit].size(); i ++) {
				(*matepairs)[orbit][i][1] += o_min_intensity;
			}
		}
		printf("finished.\n"); fflush(stdout);
	}
	if (c_min_intensity < 0) {
		printf("Minimum contig intensity: %.3f", c_min_intensity); fflush(stdout);
		printf(", correcting ... "); fflush(stdout);
		c_min_intensity = 0.0 - c_min_intensity;
		for (int contit = 0; contit < contigs->size(); contit ++) {
			for (int i = 0; i < (*contigs)[contit].size(); i ++) {
				(*contigs)[contit][i][1] += c_min_intensity;
			}
		}
		printf("finished.\n"); fflush(stdout);
	}
	
	unsigned long num_peaks = 0;
	for (unsigned int i = 0; i < matepairs->size(); i++) {
		if ((*matepairs)[i].parentCharge >= precCharge) num_peaks += (*matepairs)[i].size();
	}
	
	unsigned long num_peaks_per_thread = num_peaks / ( (unsigned long) num_threads );
	
	
	list<list<int> > mpsToAlign;
	list<int> mps;
	list<int> allMps;
	
	printf("Spawning %lu threads ...\n\n", (unsigned long)num_threads);
	
	unsigned long peaks_counted = 0;
	unsigned int num_used = 0;
	for (int i = 0; i < matepairs->size(); i++) {
		if ((*matepairs)[i].parentCharge >= precCharge && (*matepairs)[i].size() > 0) {
			mps.push_back(i);
			allMps.push_back(i);
			peaks_counted += (*matepairs)[i].size();
			
			if (peaks_counted >= num_peaks_per_thread && (num_used < num_threads - 1)) {
				num_used ++;
				printf("Thread %lu will align %lu matepairs (%lu peaks) to all contigs\n", (unsigned long)num_used, (unsigned long)mps.size(), peaks_counted);
				mpsToAlign.push_back(mps);
				mps.clear();
				peaks_counted = 0;
			}
		}
	}
	
	num_used ++;
	printf("Thread %lu will align %lu matepairs (%lu peaks) to all contigs\n\n", (unsigned long)num_used, (unsigned long)mps.size(), peaks_counted);
	mpsToAlign.push_back(mps);
	
	if (mpsToAlign.size() != num_threads) {
		cerr << "ERROR: Did not split matepair indicies properly for each thread, split for " << mpsToAlign.size() << " threads instead of " << num_threads <<  " (" << num_peaks_per_thread << " peaks per thread)" << endl;
		return;
	}
	
	MatepairContigAlignParams params;
	copy(*inparams, params);
	
	if (num_threads == 1) {
		params.mpsToAlign = allMps;
		ostringstream outF;
		outF << "_temp_mp_aligns_0.txt";
		params.outputFile = outF.str();
		getMatepairContigAlignment_thread( (void*) &params );
		
		bool success = transferTempBuffer(params.outputFile.c_str(), inparams->results);
		
		return;
	}
	
	pthread_t* threads = (pthread_t*) calloc( num_threads, sizeof(pthread_t) );
	
	vector<MatepairContigAlignParams*> paramPtrs(num_threads);
	
	
	unsigned int i = 0;
	for (list<list<int> >::iterator idxIt = mpsToAlign.begin(); idxIt != mpsToAlign.end(); idxIt++) {
		paramPtrs[i] = new MatepairContigAlignParams;
		copy(params, *(paramPtrs[i]));
		
		ostringstream outF;
		outF << "_temp_mp_aligns_" << i << ".txt";
		(paramPtrs[i])->outputFile = outF.str();
		
		paramPtrs[i]->mpsToAlign = *idxIt;
		
		pthread_create(&(threads[i]), NULL, getMatepairContigAlignment_thread, (void*) paramPtrs[i]);
		
		i++;
	}
	
	for (unsigned int i = 0; i < num_threads; i++) {
		pthread_join(threads[i], NULL);
		
		//(inparams->results).splice((inparams->results).end(), paramPtrs[i]->results);
		
		bool success = transferTempBuffer(paramPtrs[i]->outputFile.c_str(), inparams->results);
		
		delete paramPtrs[i];
	}
	
	free( threads );
}

void* getMatepairContigAlignment_thread(void* vparams) {
	MatepairContigAlignParams* params = (MatepairContigAlignParams*) vparams;
	
	Results_OCC curResult;
	
	FILE* tempBuf = fopen(params->outputFile.c_str(), "w");
	
	params->results.clear();
	int mpIdx, contigIdx;
	
	float totalCount = 0, totalWait = 0;
	
	clock_t init, final;
	int alignCount;

	for (list<int>::iterator mpIt = params->mpsToAlign.begin(); mpIt !=  params->mpsToAlign.end(); mpIt ++) {
		mpIdx = *mpIt;
		init=clock();
		alignCount = 0;
		
		for (int contigIdx = 0; contigIdx < params->contigs->size(); contigIdx++) {
			if ((*params->contigs)[contigIdx].size() < params->minNumMatchedPeaks) continue;
			
			if ( alignMatepairContig(mpIdx, contigIdx, params, curResult) )
			{
				curResult.output(tempBuf);
				alignCount ++;
			}
			
		}
		
		final=clock()-init;
		float wait = (float)final / ((float)CLOCKS_PER_SEC);
		
		totalWait += wait;
		totalCount += 1.0;
		pthread_mutex_lock(&cout_mutex);
		cout << "Processed matepair " << mpIdx << " in " << wait << " sec, found " << alignCount << " new alignments (taking " << totalWait/totalCount << " sec per matepair) - %" << 100.0*(totalCount / ((float)params->mpsToAlign.size())) << " done\n";
		pthread_mutex_unlock(&cout_mutex);
	}
	fclose(tempBuf);
}

bool transferTempBuffer( const char* buf_filename, list<Results_OCC>& results) {
	
	vector<string> strV(7);
	vector<float> floatV(7);
	
	Results_OCC curResult;
	BufferedLineReader readBuf;
	
	if (! readBuf.Load(buf_filename) )
	{
		cerr << "ERROR: Failed to load " << buf_filename << endl;
		return false;
	}
	
	for (unsigned int j = 0; j <= readBuf.size(); j++)
	{
		const char* next_line = readBuf.getline(j);
		
		if (next_line == NULL || strlen(next_line) == 0) break;
		
		bool res = splitText(next_line, strV, " ");
		
		for (unsigned int k = 0; k < strV.size(); k++)
		{
			floatV[k] = (float) atof( strV[k].c_str() );
		}
		curResult.load(floatV);
		results.push_back(curResult);
	}
	readBuf.reset();
	
	if( remove( buf_filename ) != 0 )
	{
		cerr << "Error deleting temporary buffer file " <<buf_filename << endl;
	}
	
	return true;
}

bool alignContigs(int idx1, int idx2, ContigContigAlignParams* params, Results_OCC& result) {
	Spectrum contig1 = (*(params->contigs))[idx1];
	Spectrum contig2 = (*(params->contigs))[idx2];
	float minScore = params->minScore;
	int minNumMatchedPeaks = params->minNumMatchedPeaks;
	float minAbsShift = params->minAbsShift;
	
	if (contig1.size() == 0 || contig2.size() == 0) return false;
	
	Spectrum contig2_r = contig2;
	contig2_r.reverse(0.0 - AAJumps::massH2O, 0);
	
	vector<bool> o_contig1(contig1.size());
	vector<bool> o_contig2(contig2.size());
	for (int i = 0; i < o_contig1.size(); i++) o_contig1[i] = false;
	for (int i = 0; i < o_contig2.size(); i++) o_contig2[i] = false;
	
	list<TwoValues<int> > matched;
	list<TwoValues<int> > matched_r;
	
	//map<int, list<TwoValues<int> > > shiftMP;
	//map<int, list<TwoValues<int> > > shiftMP_r;
	//map<int, list<TwoValues<int> > >::iterator shiftMPIt;
	list<TwoValues<int> >::iterator mpIt;
	float maxScore = 0.0;
	float score1, score2, overlapScore;
	float fIndex, rIndex;
	int mPeaks;
	int shiftIndex;
	float s1, s2, ratio;
	float o_m_intensity, o_t_intensity, o_a_intensity, c_m_intensity, c_t_intensity, c_a_intensity;
	float f_shift, r_shift, f_min, f_max, f_c_min, f_c_max, r_min, r_max, r_c_min, r_c_max;
	float mass, intensity;
	int o_peaks, c_peaks;
	int peaks;
	float peakTol = InputParams::PeakTol, pmTol = InputParams::PMTol;
	int intPMTol = (int)(round(pmTol/InputParams::Resolution)+0.01);
	int intPeakTol = (int)(round(peakTol/InputParams::Resolution)+0.01);
	
	float shiftF = computeBestShift(contig1, contig2, 0, minNumMatchedPeaks, matched);
	
	f_min = shiftF - pmTol;
	f_max = shiftF + contig2.parentMass + pmTol;
	f_c_min = 0.0 - shiftF - pmTol;
	f_c_max = contig1.parentMass - shiftF + pmTol;
	
	o_m_intensity = 0.0, o_t_intensity = 0.0, o_a_intensity = 0.0, c_m_intensity = 0.0, c_t_intensity = 0.0, c_a_intensity = 0.0;
	o_peaks = 0, c_peaks = 0;
	for (mpIt = matched.begin(); mpIt != matched.end(); mpIt ++) {
		o_contig1[(*mpIt)[0]] = true;
		o_contig2[(*mpIt)[1]] = true;
	}
	
	for(int i = 0; i < contig1.size(); i++ ) {
		mass = contig1[i][0]; intensity = contig1[i][1];
		if (mass >= f_min && mass <= f_max) {o_a_intensity += intensity;}
		o_t_intensity += intensity;
		if (o_contig1[i]) {
			o_m_intensity += intensity;
			o_peaks ++;
		}
		o_contig1[i] = false;
	}
	
	for (int i = 0; i < contig2.size(); i++) {
		mass = contig2[i][0]; intensity = contig2[i][1];
		if (mass >= f_c_min && mass <= f_c_max) {c_a_intensity += intensity;}
		
		c_t_intensity += intensity;
		if (o_contig2[i]) {
			c_m_intensity += intensity;
			c_peaks ++;
		}
		o_contig2[i] = false;
	}
	
	//s1 = o_m_intensity/o_t_intensity;
	//s2 = c_m_intensity/c_t_intensity;
	
	s1 = o_m_intensity/o_a_intensity;
	s2 = c_m_intensity/c_a_intensity;
	
	ratio = min(s1, s2);
	overlapScore = ratio * ((float) minNumMatchedPeaks);
	peaks = min(o_peaks, c_peaks);
	if (overlapScore > maxScore && peaks >= minNumMatchedPeaks && overlapScore >= minScore) {
		maxScore = overlapScore;
		fIndex = shiftF;
		rIndex = 0.0;
		score1 = s1;
		score2 = s2;
		mPeaks = peaks;
	}
	
	float shiftR = computeBestShift(contig1, contig2_r, 0, minNumMatchedPeaks, matched_r);
	
	r_min = shiftR - pmTol;
	r_max = shiftR + contig2.parentMass + pmTol;
	
	r_c_min = 0.0 - shiftR - pmTol;
	r_c_max = contig1.parentMass - shiftR + pmTol;
	
	o_m_intensity = 0.0, o_t_intensity = 0.0, o_a_intensity = 0.0, c_m_intensity = 0.0, c_t_intensity = 0.0, c_a_intensity = 0.0;
	o_peaks = 0, c_peaks = 0;
	for (mpIt = matched_r.begin(); mpIt != matched_r.end(); mpIt ++) {
		o_contig1[(*mpIt)[0]] = true;
		o_contig2[(*mpIt)[1]] = true;
	}

	for(int i = 0; i < contig1.size(); i++ ) {
		mass = contig1[i][0]; intensity = contig1[i][1];
		if (mass > r_min && mass < r_max ) {o_a_intensity += intensity;}
		o_t_intensity += intensity;
		if (o_contig1[i]) {
			o_m_intensity += intensity;
			o_peaks ++;
		}
		o_contig1[i] = false;
	}
	for (int i = 0; i < contig2_r.size(); i++) {
		mass = contig2_r[i][0]; intensity = contig2_r[i][1];
		
		if (mass >= r_c_min && mass <= r_c_max) {c_a_intensity += intensity;}
		
		c_t_intensity += intensity;
		if (o_contig2[i]) {
			c_m_intensity += intensity;
			c_peaks ++;
		}
		o_contig2[i] = false;
	}
	
	//s1 = o_m_intensity/o_t_intensity;
	//s2 = c_m_intensity/c_t_intensity;
	
	s1 = o_m_intensity/o_a_intensity;
	s2 = c_m_intensity/c_a_intensity;
	
	ratio = min(s1, s2);
	overlapScore = ratio * ((float) minNumMatchedPeaks);
	peaks = min(o_peaks, c_peaks);
	if (overlapScore > maxScore && peaks >= minNumMatchedPeaks && overlapScore >= minScore) {
		maxScore = overlapScore;
		fIndex = shiftR;
		rIndex = 1.0;
		score1 = s1;
		score2 = s2;
		mPeaks = peaks;
	}
	
	
	if (maxScore > 0) {
		result.spec1 = idx1;
		result.spec2 = idx2;
		result.score1 = score1;
		result.score2 = score2;
		result.shift1 = fIndex;
		result.shift2 = rIndex;
		result.matchedPeaks = mPeaks;
		return true;
	}
	return false;
}


void getContigContigAlignment(ContigContigAlignParams* inparams, unsigned int num_threads) {
	inparams->results.clear();
	SpecSet* contigs = inparams->contigs;
	float minScore = inparams->minScore;
	int minNumMatchedPeaks = inparams->minNumMatchedPeaks;
	float minAbsShift = inparams->minAbsShift;
	AAJumps allowedJumps = *(inparams->allowedJumps);
	
	if (num_threads == 0) {
		cerr << "ERROR: num_threads must be >= 1\n";
		return;
	}
	
	float c_min_intensity = 0;
	for (int contit = 0; contit < contigs->size(); contit ++) {
		if ((*contigs)[contit].size() == 0) continue;
		(*contigs)[contit].parentMass = (*contigs)[contit].peakList.back()[0] + AAJumps::massMH;
		for (int i = 0; i < (*contigs)[contit].size(); i ++) {
			if ((*contigs)[contit][i][1] < c_min_intensity) c_min_intensity = (*contigs)[contit][i][1];
		}
	}

	if (c_min_intensity < 0) {
		printf("Minimum contig intensity: %.3f", c_min_intensity); fflush(stdout);
		printf(", correcting ... "); fflush(stdout);
		c_min_intensity = 0.0 - c_min_intensity;
		for (int contit = 0; contit < contigs->size(); contit ++) {
			for (int i = 0; i < (*contigs)[contit].size(); i ++) {
				(*contigs)[contit][i][1] += c_min_intensity;
			}
		}
		printf("finished.\n"); fflush(stdout);
	}
	
	map<int, set<int> >* seen = new map<int, set<int> >;
	set<int>* nodes = new set<int>;
	
	list<TwoValues<int> > allPairs;
	TwoValues<int> pair;
	
	for (int contit = 0; contit < contigs->size(); contit ++) {
		if ((*contigs)[contit].size() == 0) continue;
		
		for (int contit2 = 0; contit2 < contigs->size(); contit2 ++) {
			if (contit2 == contit) continue;
			if (seen->count(contit) > 0 && (*seen)[contit].count(contit2) > 0) continue;
			if (seen->count(contit) == 0) {
				nodes->clear(); nodes->insert(contit2);
				(*seen)[contit] = *nodes;
			} else (*seen)[contit].insert(contit2);
			if (seen->count(contit2) == 0) {
				nodes->clear(); nodes->insert(contit);
				(*seen)[contit2] = *nodes;
			} else (*seen)[contit2].insert(contit);
			
			pair[0] = contit;
			pair[1] = contit2;
			allPairs.push_back(pair);
		}
	}
	
	delete seen;
	delete nodes;
	
	ContigContigAlignParams params;
	copy(*inparams, params);
	
	if (num_threads == 1) {
		params.pairsToCompute = allPairs;
		getContigContigAlignment_thread( (void*) &params );
		inparams->results = params.results;
		return;
	}
	
	unsigned int pairsPerThread = allPairs.size() / num_threads;
	
	printf("Spawning %lu threads ...\n\n", (unsigned long)num_threads);
	
	pthread_t* threads = (pthread_t*) calloc( num_threads, sizeof(pthread_t) );
	
	vector<ContigContigAlignParams*> paramPtrs(num_threads);
	
	for (unsigned int i = 0; i < num_threads; i++) {
		paramPtrs[i] = new ContigContigAlignParams;
		copy(params, *(paramPtrs[i]));
		
		list<TwoValues<int> >::iterator startIt = allPairs.begin();
		list<TwoValues<int> >::iterator endIt = allPairs.begin();
		for ( unsigned int j = 0; j < pairsPerThread && endIt != allPairs.end(); j++) {
			++endIt;
		}
		
		(paramPtrs[i]->pairsToCompute).splice((paramPtrs[i]->pairsToCompute).begin(), allPairs, startIt, endIt);
		
		pthread_create(&(threads[i]), NULL, getContigContigAlignment_thread, (void*) paramPtrs[i]);
	}
	
	for (unsigned int i = 0; i < num_threads; i++) {
		pthread_join(threads[i], NULL);
		
		(inparams->results).splice((inparams->results).end(), paramPtrs[i]->results);
		
		delete paramPtrs[i];
	}
	
	free( threads );
}

void copy(ContigContigAlignParams& from, ContigContigAlignParams& to) {
	to.contigs = from.contigs;
	to.minScore = from.minScore;
	to.minNumMatchedPeaks = from.minNumMatchedPeaks;
	to.minAbsShift = from.minAbsShift;
	to.allowedJumps = from.allowedJumps;
}

void* getContigContigAlignment_thread(void* vparams) {
	ContigContigAlignParams* params = (ContigContigAlignParams*) vparams;
	
	Results_OCC curResult;
	
	//clock_t init, final;
	params->results.clear();
	int idx1, idx2;
	
	for (list<TwoValues<int> >::iterator contit = params->pairsToCompute.begin(); contit !=  params->pairsToCompute.end(); contit ++) {
		idx1 = (*contit)[0];
		idx2 = (*contit)[1];
		
		if ( alignContigs(idx1, idx2, params, curResult) ) {
			params->results.push_back(curResult);
			
			pthread_mutex_lock(&cout_mutex);
			curResult.output(cout, ' ');
			pthread_mutex_unlock(&cout_mutex);
		}
	}
}

inline void mergeVectors(vector<int> &v1, vector<int> &v2, vector<int> &vOut) {
	unsigned int i1=0, i2=0, iOut=0;
	vOut.resize(v1.size()+v2.size());
	while(i1<v1.size() or i2<v2.size()) {
		if(i1==v1.size()) { vOut[iOut++]=v2[i2++]; continue; }
		if(i2==v2.size()) { vOut[iOut++]=v1[i1++]; continue; }
		if(v1[i1]==v2[i2]) { vOut[iOut++]=v1[i1++]; i2++; continue; }
		if(v1[i1]<v2[i2]) vOut[iOut++]=v1[i1++]; else vOut[iOut++]=v2[i2++];
	}
	vOut.resize(iOut);
}

// aux_updateMeanVariance - Updates mean/variance values to also include newValue.
//
// IMPORTANT NOTE: sampleSize in NOT UPDATED because multiple different means/variances may
//                    be updated on the same sample (e.g. num peaks, matched intensities) and
//                    sampleSize should only be updated once.
inline void aux_updateMeanVariance(float newValue, float sampleSize, float &mean, float &variance) {
	float updateRatio = sampleSize/(sampleSize+1);   sampleSize++;
	mean = mean*updateRatio + newValue/sampleSize;
	variance = updateRatio*variance+(newValue*newValue)/sampleSize;
}

//
//  getPairAlignsPA2 - Like getPairAlignsPA but fixes some problems and adds some extensions
//    - Each shift is only computed once
//    - Shift pairs have scores precomputed by computeShifts and are passed to ScoreOverlap6 only if
//        their potential DP score is higher than the best achieved DP score for the pair
//    - Only returns pairs that meet the minRatio/minPeakAreaOvlp criteria
//
//  alignStats - Alignment statistics per returned pair (each is minimum between the 2 aligned spectra):
//                  pos 0: num matched peaks (before DP/merging)
//                  pos 1: matched intensity (before DP/merging)
//                  pos 2: num matched peaks (after DP/merging)
//                  pos 3: matched intensity (after DP/merging)
//  specStats  - Alignment statistics per aligned spectrum: (intensity not considered because intensity distributions should vary per spectrum)
//                  pos 0: num alignments
//                  pos 1/2: mean/variance for num matched peaks (before DP/merging)
//                  pos 3/4: mean/variance num matched peaks (after DP/merging)
//                  pos 5/6: mean/variance for matched intensity (before DP/merging)
//                  pos 7/8: mean/variance matched intensity (after DP/merging)

void getPairAlignsPA2(SpecSet &specSet, unsigned int startIdx, unsigned int endIdx,
						float peakTol, float pmTol, float minRatio, float minPeakAreaOvlp,
						short minNumMatchedPeaks, AAJumps &allowedJumps, float minAbsShift,
						list<Results_PA> &results, list<TwoValues<float> > &ratios, list<TwoValues<int> > &numMatchedPeaks,
						vector<TwoValues<float> > &means, vector<float> &varTerms,
						list<vector<float> > &alignStats, vector<vector<float> > &specStats){
	unsigned int spec1, spec2;

	// Variables and initializations
	TwoValues<int> res;
	list<float> shiftScores;    // Maximum scores for eligible pairs of shifts
	list<TwoValues<unsigned int> > shiftPairs;     // List of eligible pairs of shifts (in the same order as shiftScores)
//	vector<vector<TwoValues<int> > > shiftMatchedPeaks;  // Lists of matched peaks per shift between the two spectra (indexed using
	vector<list<TwoValues<int> > > shiftMatchedPeaks;  // Lists of matched peaks per shift between the two spectra (indexed using
												   //  a shift offset returned by computeShifts)
	vector<TwoValues<float> > shiftDPscores;       // Match score per shift after applying dynamic programming
	vector<float> shiftDPpenalties;                // Peak mass mismatch penalties (in DP)
	vector<TwoValues<vector<int> > > shiftDPmatchedPeaks;  // Peaks matched by the DP algorithm
	vector<TwoValues<float> > minMaxMatchScores(specSet.size());  // Minimum acceptable / Maximum possible match scores per spectrum in the dataset
	Results_PA curResult;
	vector<float> curAlignStats(4);

	// Initialize stats vectors
	means.resize(specSet.size());       varTerms.resize(specSet.size());     specStats.resize(specSet.size());
	for(unsigned int i=0; i<means.size(); i++) {
		means[i].set(0,0); varTerms[i]=0;
		specStats[i].resize(9); for(unsigned int j=0; j<9; j++) specStats[i][j]=0;
	}

	float maxParentMass = 0;  unsigned int maxNumPeaks=0;
	int intPeakTol = (int)round(peakTol/InputParams::Resolution), intPMTol = (int)round(pmTol/InputParams::Resolution);
	for(unsigned int i=0; i<specSet.size(); i++) {
		if(specSet[i].parentMass>maxParentMass) maxParentMass=specSet[i].parentMass;
		if(specSet[i].size()>maxNumPeaks) maxNumPeaks=specSet[i].size();
		minMaxMatchScores[i].set(0,0); for(unsigned int j=0; j<specSet[i].size();j++) minMaxMatchScores[i][1]+=specSet[i][j][1];
		minMaxMatchScores[i][0]=minRatio*minMaxMatchScores[i][1];
	}
//	shiftMatchedPeaks.resize(1+(int)round(2*(min(maxParentMass,InputParams::MaxShift)+8*peakTol)/InputParams::Resolution)); //for(unsigned int i=0;i<shiftMatchedPeaks.size();i++) shiftMatchedPeaks[i].reserve(10);
	shiftMatchedPeaks.resize(1+(int)ceil(2*(min(maxParentMass,InputParams::MaxShift)+peakTol+pmTol)/InputParams::Resolution)); //for(unsigned int i=0;i<shiftMatchedPeaks.size();i++) shiftMatchedPeaks[i].reserve(10);
	shiftDPscores.reserve(shiftMatchedPeaks.capacity());   shiftDPpenalties.reserve(shiftMatchedPeaks.capacity());
	shiftDPmatchedPeaks.reserve(shiftMatchedPeaks.capacity());
	vector<int> idx1, idx2; idx1.reserve(maxNumPeaks*(2*intPeakTol+1)+1); idx2.reserve(idx1.capacity());  // Temporary variables used as input to ScoreOverlap6
	vector<int> idxMerged1, idxMerged2; idxMerged1.reserve(idx1.capacity()); idxMerged2.reserve(idx1.capacity());  // Temporary variables used to merge symmetric results of ScoreOverlap6

	// Compute the pairwise alignments
	time_t curTime=time(0); double totTime=0.0;  int numResults=0;
	float bestShiftScore;
	TwoValues<int> bestShiftPair;
	TwoValues<int> bestNumMatchedPeaks;
	TwoValues<float> bestShiftScores, bestCandidateScores;
	list<int> shiftsToCompute;
	for(spec1=startIdx; spec1<=endIdx; spec1++) {
		cout << "Processing spectrum "<<spec1<<"... "; cout.flush();
		for(spec2=spec1+1; spec2<specSet.size(); spec2++) {
			for(unsigned int i=0; i<curAlignStats.size(); i++) curAlignStats[i]=0;

			// Step 1: Get the list of possible shifts
			res = computeShifts(specSet[spec1], specSet[spec2], peakTol, pmTol, minRatio, minPeakAreaOvlp, minNumMatchedPeaks, allowedJumps, shiftScores, shiftPairs, shiftMatchedPeaks, bestCandidateScores, bestNumMatchedPeaks, minAbsShift);
			if(res[1]==0) continue;  // No eligible shifts found; could happen when maximum overlap is less than minPeakAreaOvlp

			// Matched intensity before DP/merging
			aux_updateMeanVariance(bestCandidateScores[0],specStats[spec1][0],specStats[spec1][5],specStats[spec1][6]);
			aux_updateMeanVariance(bestCandidateScores[1],specStats[spec2][0],specStats[spec2][5],specStats[spec2][6]);
			// Num matched peaks before DP/merging
			aux_updateMeanVariance(bestNumMatchedPeaks[0],specStats[spec1][0],specStats[spec1][1],specStats[spec1][2]);
			aux_updateMeanVariance(bestNumMatchedPeaks[1],specStats[spec2][0],specStats[spec2][1],specStats[spec2][2]);
			// Pre-DP alignment stats, in case this alignment gets reported
			curAlignStats[0] = min(bestNumMatchedPeaks[0],bestNumMatchedPeaks[1]);
			curAlignStats[1] = min(bestCandidateScores[0]/minMaxMatchScores[spec1][1],bestCandidateScores[1]/minMaxMatchScores[spec2][1]);

/*					means[spec1][0]=(means[spec1][0]*means[spec1][1]+bestCandidateScores[0])/(means[spec1][1]+1);  // Include accepted pairs in means
					means[spec2][0]=(means[spec2][0]*means[spec2][1]+bestCandidateScores[1])/(means[spec2][1]+1);
					means[spec1][1]++;    means[spec2][1]++;
					float updateRatio = (means[spec1][1]-1)/means[spec1][1];  // Update variance terms
					varTerms[spec1] = updateRatio*varTerms[spec1]+(bestCandidateScores[0]*bestCandidateScores[0])/means[spec1][1];
					updateRatio = (means[spec2][1]-1)/means[spec2][1];
					varTerms[spec2] = updateRatio*varTerms[spec2]+(bestCandidateScores[1]*bestCandidateScores[1])/means[spec2][1];
*/

//cerr<<"shiftDPscores resized to "<<shiftMatchedPeaks.size()<<endl;
			shiftDPscores.resize(shiftMatchedPeaks.size()); for(unsigned int i=0; i<shiftDPscores.size(); i++) shiftDPscores[i].set(-1.0,-1.0);
			shiftDPpenalties.resize(shiftMatchedPeaks.size()); for(unsigned int i=0; i<shiftDPpenalties.size(); i++) shiftDPpenalties[i]=0.0;
			shiftDPmatchedPeaks.resize(shiftMatchedPeaks.size());
			for(unsigned int i=0; i<shiftDPmatchedPeaks.size(); i++) { shiftDPmatchedPeaks[i][0].resize(0);  shiftDPmatchedPeaks[i][1].resize(0); }

			bestShiftScore = 0;   bestShiftPair.set(-1,-1);   bestShiftScores.set(0,0);
			bestCandidateScores.set(0,0);   bestNumMatchedPeaks.set(0,0); // Reuse to store best DP match score / num matched peaks (used to update means/vars)
			list<float>::iterator nextShiftScore = shiftScores.begin();
			list<TwoValues<unsigned int> >::iterator nextShiftPair = shiftPairs.begin();
			while(nextShiftScore!=shiftScores.end() and *nextShiftScore > bestShiftScore) {
				int shiftIndex=(*nextShiftPair)[0], shiftSym=(*nextShiftPair)[1];
//cerr<<"Shift pair "<<shiftIndex<<"/"<<shiftSym<<endl;
				if(shiftDPscores[shiftIndex][0]+shiftDPscores[shiftIndex][1]<0.0001 and shiftMatchedPeaks[shiftIndex].size()>0) shiftsToCompute.push_back(shiftIndex);
				for(int tolIdx=max(0,shiftSym-intPMTol); tolIdx<=shiftSym+intPMTol and tolIdx<(int)shiftDPscores.size(); tolIdx++)
					if(shiftDPscores[tolIdx][0]+shiftDPscores[tolIdx][1]<0.0001 and shiftMatchedPeaks[tolIdx].size()>0) shiftsToCompute.push_back(tolIdx);

				for(list<int>::iterator curShift = shiftsToCompute.begin(); curShift != shiftsToCompute.end(); curShift++) {
					float shiftMass = ((*curShift)-res[0])*InputParams::Resolution;

					idx1.resize(shiftMatchedPeaks[*curShift].size());   idx2.resize(shiftMatchedPeaks[*curShift].size());  int pivot=0;
					for(list<TwoValues<int> >::iterator iterPeaks=shiftMatchedPeaks[*curShift].begin(); iterPeaks!=shiftMatchedPeaks[*curShift].end(); iterPeaks++)
						{ idx1[pivot]=(*iterPeaks)[0]; idx2[pivot]=(*iterPeaks)[1]; pivot++; }

/*cerr<<"shift = "<<*curShift<<", idx1=(";
for(unsigned int dbgIdx=0; dbgIdx<idx1.size(); dbgIdx++) cerr<<idx1[dbgIdx]<<",";
cerr<<"), idx2=(";
for(unsigned int dbgIdx=0; dbgIdx<idx2.size(); dbgIdx++) cerr<<idx2[dbgIdx]<<",";
cerr<<")\n";*/
					ScoreOverlap6mp(specSet[spec1], idx1, specSet[spec2], idx2, shiftMass, peakTol, shiftDPmatchedPeaks[*curShift][0], shiftDPmatchedPeaks[*curShift][1], AAJumps::minAAmass, &shiftDPpenalties[*curShift]);
//					ScoreOverlap6(specSet[spec1], idx1, specSet[spec2], idx2, shiftMass, peakTol, shiftDPmatchedPeaks[*curShift][0], shiftDPmatchedPeaks[*curShift][1], false);

/*cerr<<"shift = "<<*curShift<<", matched1=(";
for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[*curShift][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[*curShift][0][dbgIdx]<<",";
cerr<<"), matched2=(";
for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[*curShift][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[*curShift][1][dbgIdx]<<",";
cerr<<")\n";*/
					shiftDPscores[*curShift][0]=0; for(unsigned int i=0;i<shiftDPmatchedPeaks[*curShift][0].size();i++) shiftDPscores[*curShift][0]+=specSet[spec1][shiftDPmatchedPeaks[*curShift][0][i]][1];
					shiftDPscores[*curShift][1]=0; for(unsigned int i=0;i<shiftDPmatchedPeaks[*curShift][1].size();i++) shiftDPscores[*curShift][1]+=specSet[spec2][shiftDPmatchedPeaks[*curShift][1][i]][1];
				}
				shiftsToCompute.clear();

//cerr<<"--- repeat shift pair "<<shiftIndex<<"/"<<shiftSym<<endl;
				for(int tolIdx=max(0,shiftSym-intPMTol); tolIdx<=shiftSym+intPMTol and tolIdx<(int)shiftDPmatchedPeaks.size(); tolIdx++) {
//					if(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0]+shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]>bestCandidateScores[0]+bestCandidateScores[1])
//						bestCandidateScores.set(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0],shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]);
//					if(shiftDPscores[shiftIndex][0]+shiftDPscores[tolIdx][0]<minMaxMatchScores[spec1][0] or  // Sole purpose is to avoid unnecessary computations
//					   shiftDPscores[shiftIndex][1]+shiftDPscores[tolIdx][1]<minMaxMatchScores[spec2][0]) continue;
					if(shiftMatchedPeaks[tolIdx].size()==0) continue;

					// merge lists of matched peaks - spec1
					mergeVectors(shiftDPmatchedPeaks[shiftIndex][0], shiftDPmatchedPeaks[tolIdx][0], idxMerged1);
					mergeVectors(shiftDPmatchedPeaks[shiftIndex][1], shiftDPmatchedPeaks[tolIdx][1], idxMerged2);

/*cerr<<"shifts = "<<shiftIndex<<"/"<<tolIdx<<":\n  idxMerged1: (";
for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[shiftIndex][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[shiftIndex][0][dbgIdx]<<",";
cerr<<") U (";
for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[tolIdx][0].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[tolIdx][0][dbgIdx]<<",";
cerr<<") -> (";
for(unsigned int dbgIdx=0; dbgIdx<idxMerged1.size(); dbgIdx++) cerr<<idxMerged1[dbgIdx]<<",";
cerr<<")\n  idxMerged2: (";
for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[shiftIndex][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[shiftIndex][1][dbgIdx]<<",";
cerr<<") U (";
for(unsigned int dbgIdx=0; dbgIdx<shiftDPmatchedPeaks[tolIdx][1].size(); dbgIdx++) cerr<<shiftDPmatchedPeaks[tolIdx][1][dbgIdx]<<",";
cerr<<") -> (";
for(unsigned int dbgIdx=0; dbgIdx<idxMerged2.size(); dbgIdx++) cerr<<idxMerged2[dbgIdx]<<",";
cerr<<")\n";*/
					float score1=0;   for(unsigned int i=0;i<idxMerged1.size();i++) score1+=specSet[spec1][idxMerged1[i]][1];
					float score2=0;   for(unsigned int i=0;i<idxMerged2.size();i++) score2+=specSet[spec2][idxMerged2[i]][1];
// Moved (and changed) from above on 2006/09/11
					if(bestShiftScore<0.0001) {  // These statements ensure that we get values to update means/variances
						if(score1+score2>bestCandidateScores[0]+bestCandidateScores[1])	bestCandidateScores.set(score1,score2); // Used to update means/vars even if this pair is not returned
						if(idxMerged1.size()+idxMerged2.size()>bestNumMatchedPeaks[0]+bestNumMatchedPeaks[1]) bestNumMatchedPeaks.set(idxMerged1.size(),idxMerged2.size());
					}

					if(score1<minMaxMatchScores[spec1][0] or score2<minMaxMatchScores[spec2][0]) continue;
					if(score1+score2-shiftDPpenalties[shiftIndex]-shiftDPpenalties[tolIdx]>bestShiftScore) {
						bestShiftScore=score1+score2-shiftDPpenalties[shiftIndex]-shiftDPpenalties[tolIdx];
						bestShiftPair.set(shiftIndex,tolIdx);
						bestShiftScores.set(score1,score2);
						bestNumMatchedPeaks.set(idxMerged1.size(),idxMerged2.size());
					}
				}
				nextShiftScore++;   nextShiftPair++;
			}

			// Build results for best match
			TwoValues<float> curRatios;
			if(bestShiftScore>0) {
				curResult.spec1 = spec1;                  curResult.spec2 = spec2;
				curResult.score1 = bestShiftScores[0];    curResult.score2 = bestShiftScores[1];
				curResult.shift1 = (bestShiftPair[0]-res[0])*InputParams::Resolution;
				curResult.shift2 = (bestShiftPair[1]-res[0])*InputParams::Resolution;
				results.push_back(curResult);
				curRatios.set(curResult.score1/minMaxMatchScores[spec1][1],curResult.score2/minMaxMatchScores[spec2][1]);
				ratios.push_back(curRatios);
				numMatchedPeaks.push_back(bestNumMatchedPeaks);
				bestCandidateScores = bestShiftScores;

				// Post-DP alignment stats since this alignment is reported
				curAlignStats[2] = min(bestNumMatchedPeaks[0],bestNumMatchedPeaks[1]);
				curAlignStats[3] = min(bestCandidateScores[0]/minMaxMatchScores[spec1][1],bestCandidateScores[1]/minMaxMatchScores[spec2][1]);
				alignStats.push_back(curAlignStats);
			}

			// Regular mean/vars output to text files
			aux_updateMeanVariance(bestCandidateScores[0],means[spec1][1],means[spec1][0],varTerms[spec1]);
			aux_updateMeanVariance(bestCandidateScores[1],means[spec2][1],means[spec2][0],varTerms[spec2]);
			// Matched intensity after DP/merging
			aux_updateMeanVariance(bestCandidateScores[0],specStats[spec1][0],specStats[spec1][7],specStats[spec1][8]);
			aux_updateMeanVariance(bestCandidateScores[1],specStats[spec2][0],specStats[spec2][7],specStats[spec2][8]);
			// Num matched peaks after DP/merging
//cerr<<"Num matched peaks: "<<bestNumMatchedPeaks[0]<<"/"<<bestNumMatchedPeaks[1]<<", sampleSize = "<<specStats[spec1][0]<<"/"<<specStats[spec2][0];
			aux_updateMeanVariance(bestNumMatchedPeaks[0],specStats[spec1][0],specStats[spec1][3],specStats[spec1][4]);
			aux_updateMeanVariance(bestNumMatchedPeaks[1],specStats[spec2][0],specStats[spec2][3],specStats[spec2][4]);
//cerr<<", updated means = "<<specStats[spec1][3]<<"/"<<specStats[spec2][3]<<endl;
			specStats[spec1][0]++;   means[spec1][1]++;
			specStats[spec2][0]++;   means[spec2][1]++;

/*			means[spec1][0]=(means[spec1][0]*means[spec1][1]+bestCandidateScores[0])/(means[spec1][1]+1);  // Include accepted pairs in means
			means[spec2][0]=(means[spec2][0]*means[spec2][1]+bestCandidateScores[1])/(means[spec2][1]+1);
			means[spec1][1]++;    means[spec2][1]++;
			float updateRatio = (means[spec1][1]-1)/means[spec1][1];  // Update variance terms
			varTerms[spec1] = updateRatio*varTerms[spec1]+(bestCandidateScores[0]*bestCandidateScores[0])/means[spec1][1];
			updateRatio = (means[spec2][1]-1)/means[spec2][1];
			varTerms[spec2] = updateRatio*varTerms[spec2]+(bestCandidateScores[1]*bestCandidateScores[1])/means[spec2][1];
*/	    }
		time_t newTime = time(0); double ellapsed=difftime(newTime,curTime); totTime+=ellapsed; curTime=newTime;  int N=endIdx,CUR=spec1-startIdx+1;
		cout<< ellapsed << " secs, number of new pair aligns = "<<results.size()-numResults<<" / "<<results.size()<<". Average time = "<<totTime/CUR<<", ETA = "<< (N-CUR-1)*(N-CUR-2)*totTime/(CUR*(N-CUR)+(CUR-1)*CUR/2.0) << endl; cout.flush();
		numResults = results.size();
	}
	computeShifts_cleanup();
}

//
// getPairAlignsPAext_aux - Computes alignment of one MA against one PRM spectrum
//
// Shifts are computed at the same resolution as that of the input spectra, i.e.
//   shifts lambda = round(prmSpec[p][0]-cSpec[p'][0]). Convert
//   input spectra to 0.1 resolution before calling this function to get shifts
//   with 0.1 Da resolution.
//
//   WARNING: result.score2 (for consensus spectrum) is from both cSpec and
//      cSpecRev - compute ratios by result.score2/(2*totalScoreInCSpec)\n";
//
void getPairAlignsPAext_aux(Spectrum &prmSpec, Spectrum &cSpec, Spectrum &cSpecRev,
							short aaDiff, float pmTol, float peakTol,
							Results_PA &result){
	vector<TwoValues<float> > preShiftScores;     // Max. shift scores before the DP: consensus spec (col 0), reversed spec (col 1)
	vector<TwoValues<Spectrum*> > preShiftPeaks;  // Matched peaks' scores before the DP
	vector<TwoValues<float> > finalShiftScores;   // Shift scores after the DP
	int shift, shiftSym;
	Spectrum specToScore;      // Filtered version of prmSpec that actually gets scored
	vector<int> idxMatched;  // List of matched peaks in specToScore (after DP)

cerr<<"WARNING: batch.cpp::getPairAlignsPAext_aux() was not revised to handle high precision data!\n";

	// Initialize specToScore
	specToScore.parentMass=prmSpec.parentMass;	specToScore.peakList.reserve(prmSpec.size());

// These steps are common to every alignment and should therefore be computed outside
//  this function and reused accross all alignments

	int peakTolInt = (int)round(peakTol);
	int pmTolInt = (int)round(pmTol);
	int middleShift = (int)round((prmSpec.parentMass-cSpec.parentMass)/2);
	int middleTimesTwo = (int)round(prmSpec.parentMass-cSpec.parentMass);

	// Step 1: Compute possible shifts
	int shiftsOffset = (int)ceil(cSpec.parentMass);
	float peakDiff;   // Tracks mass difference between matching peaks
	preShiftScores.resize(shiftsOffset+(int)ceil(prmSpec.parentMass)+1);
	preShiftPeaks.resize(shiftsOffset+(int)ceil(prmSpec.parentMass)+1);
	finalShiftScores.resize(shiftsOffset+(int)ceil(prmSpec.parentMass)+1);
	int highShift;   float score;
	for(int i=0; i<cSpec.size(); i++) {
		for(int j=0; j<prmSpec.size(); j++) {
			highShift = shiftsOffset+(int)prmSpec[j][0]-(int)cSpec[i][0]+peakTolInt;   score = cSpec[i][1]+prmSpec[j][1];
			peakDiff = -peakTolInt;
			for(int p=shiftsOffset+(int)prmSpec[j][0]-(int)cSpec[i][0]-peakTolInt; p<=highShift; p++, peakDiff++)
				preShiftScores[p][0]+=score-fabs(peakDiff)/(double)peakTolInt;

			// correspondingly for the peak matches with the reversed spectrum
			if(i<cSpecRev.size()) {
				highShift = shiftsOffset+(int)prmSpec[j][0]-(int)cSpecRev[i][0]+peakTolInt;	  score = cSpecRev[i][1]+prmSpec[j][1];
				peakDiff = -peakTolInt;
				for(int p=shiftsOffset+(int)prmSpec[j][0]-(int)cSpecRev[i][0]-peakTolInt; p<=highShift; p++, peakDiff++)
					preShiftScores[p][1]+=score-fabs(peakDiff)/(double)peakTolInt;
			}
		}
	}
//cerr<<"Shifts with cSpec:\n";
//for(int i=0; i<preShiftScores.size(); i++) if(preShiftScores[i][0]>0) {if(i>0 && preShiftScores[i-1][0]==0)cerr<<"\n"; cerr<<"["<<(i-shiftsOffset)<<","<<preShiftScores[i][0]<<"],"; } cerr<<endl;
//cerr<<"Shifts with cSpecRev:\n";
//for(int i=0; i<preShiftScores.size(); i++) if(preShiftScores[i][1]>0) {if(i>0 && preShiftScores[i-1][1]==0)cerr<<"\n"; cerr<<"["<<(i-shiftsOffset)<<","<<preShiftScores[i][1]<<"],"; } cerr<<endl;

	// Step 2: Find the the most promising pair of shifts

// TODO Here:
//   - Filter out shift pairs whose maximum score doesn't meet the minRatio criterion
//   - Create sorted list of shift pairs (by decreasing maximum possible score)
//      -> use to stop computing matches when maximum possible score < found DP score

	float preDPscore=0;  TwoValues<int> preDPshifts;
	int tolMax=pmTolInt, tolMin=-tolMax, t;
	for(shift=0; shift<preShiftScores.size(); shift++) {
		shiftSym = shiftsOffset+(middleTimesTwo-(shift-shiftsOffset));
		for(t=max(0,shiftSym+tolMin); t<shiftSym+tolMax && t<preShiftScores.size(); t++)
			if(preShiftScores[shift][0]+preShiftScores[t][1]>preDPscore) {
				preDPscore = preShiftScores[shift][0]+preShiftScores[t][1];
				preDPshifts.set(shift,t);
//cerr<<"Updated preDPshifts to ["<<(shift-shiftsOffset)<<","<<(t-shiftsOffset)<<"], score = "<<preDPscore<<endl;
			}
	}

	//
	// Step 3: Initialize the best match with the most promising pair of shifts
	//
	float score1, score2; // spec 1 is prmSpec, spec 2 is consensus spec
	vector<float> filteredCSpec,    // Scores of the cSpec peaks used in the first filtered spectrum
				  filteredCSpecRev, // Scores of the cSpecRev peaks used in the second filtered spectrum
				  matchCSpec;       // Merged scores of the consensus peaks used when matching prmSpec

	// Get the filtered versions of prmSpec
	preShiftPeaks[preDPshifts[0]][0] = new Spectrum();
	preShiftPeaks[preDPshifts[1]][1] = new Spectrum();

	float endpointScores=0;
	prmSpec.filterAdd(cSpec, (preDPshifts[0]-shiftsOffset), peakTolInt, (*preShiftPeaks[preDPshifts[0]][0]), filteredCSpec, endpointScores);
//cerr<<"filteredCSpec:\n"; for(int i=0; i<filteredCSpec.size(); i++) cerr<<(*preShiftPeaks[preDPshifts[0]][0])[i][0]<<"\t"<<(*preShiftPeaks[preDPshifts[0]][0])[i][1]<<"\t["<<filteredCSpec[i]<<"]\n";
	prmSpec.filterAdd(cSpecRev, (preDPshifts[1]-shiftsOffset), peakTolInt, (*preShiftPeaks[preDPshifts[1]][1]), filteredCSpecRev, endpointScores);
//cerr<<"filteredCSpecRev:\n"; for(int i=0; i<filteredCSpecRev.size(); i++) cerr<<(*preShiftPeaks[preDPshifts[1]][1])[i][0]<<"\t"<<(*preShiftPeaks[preDPshifts[1]][1])[i][1]<<"\t["<<filteredCSpecRev[i]<<"]\n";
	// Merge the filtered versions to construct specToScore
	preShiftPeaks[preDPshifts[0]][0]->merge(preShiftPeaks[preDPshifts[1]][1]->peakList,specToScore,filteredCSpec,filteredCSpecRev,matchCSpec);
//cerr<<"specToScore:\n"; for(int i=0; i<specToScore.size(); i++) cerr<<specToScore[i][0]<<"\t"<<specToScore[i][1]<<"\n";
	// Find matched peaks
	getMaxSparseSet(specToScore, peakTol, 0, idxMatched);
//cerr<<"idxMatched:"; for(int i=0; i<idxMatched.size(); i++) cerr<<idxMatched[i]<<", "; cerr<<"\n";
	// Set the result variables
	result.shift1=(preDPshifts[0]-shiftsOffset), result.shift2=(preDPshifts[1]-shiftsOffset);
	result.score1=0; result.score2=0;
	for(int i=0; i<idxMatched.size(); i++) {
		result.score1+=specToScore[idxMatched[i]][1];
		result.score2+=matchCSpec[idxMatched[i]];
	}
	result.score1-=result.score2;  // result.score1 includes result.score2 because the peak scores were summed into specToScore

	//
	// Step 4: Go over all remaining shifts with preShiftScores > bestScore1
	//
	Spectrum prevMergedSpec;
	for(shift=0; shift<preShiftScores.size(); shift++) {
		preShiftPeaks[shift][0] = new Spectrum();
		shiftSym = shiftsOffset+(middleTimesTwo-(shift-shiftsOffset));
		for(t=max(0,shiftSym+tolMin); t<shiftSym+tolMax && t<preShiftScores.size(); t++)
			if(preShiftScores[shift][0]+preShiftScores[t][1]>0  &&
			   preShiftScores[shift][0]+preShiftScores[t][1]>result.score1+result.score2) {
//				preDPscore = preShiftScores[shift][0]+preShiftScores[t][1];
//				preDPshifts.set(shift,t);
				if(preShiftPeaks[t][1]==0) preShiftPeaks[t][1] = new Spectrum();

//	Do something about this prevMergedSpectrum issue when shift is the same? Maybe just ignore?
//				if((int)preShiftPeaks[shift][0]==0){
//					preShiftPeaks[shift][0] = new Spectrum();
//				} else
//					prmSpec.filterAdd(cSpec, (shift-shiftsOffset), peakTolInt, mergeSpec, filteredCSpec, endpointScores);

				prmSpec.filterAdd(cSpec, (shift-shiftsOffset), peakTolInt, (*preShiftPeaks[shift][0]), filteredCSpec, endpointScores);
				prmSpec.filterAdd(cSpecRev, (t-shiftsOffset), peakTolInt, (*preShiftPeaks[t][1]), filteredCSpecRev, endpointScores);
				preShiftPeaks[shift][0]->merge(preShiftPeaks[t][1]->peakList,specToScore,filteredCSpec,filteredCSpecRev,matchCSpec);
				// Check if peakLists are different between current and previous shifts
				// If they're the same then the outcome must be the same and the remaining code can be skipped
				if(shift>0 && preShiftPeaks[shift-1][0]>0 && preShiftPeaks[shift][0]->compare(*preShiftPeaks[shift-1][0])) continue;
//cerr << "Computing shift pair ["<<shift-shiftsOffset<<","<<t-shiftsOffset<<"]\n";
//cerr<<"specToScore:\n"; for(int i=0; i<specToScore.size(); i++) cerr<<specToScore[i][0]<<"\t"<<specToScore[i][1]<<"\n";
				getMaxSparseSet(specToScore, peakTol, 0, idxMatched);

				score1=0; score2=0;
				for(int i=0; i<idxMatched.size(); i++) {
					score1+=specToScore[idxMatched[i]][1];
					score2+=matchCSpec[idxMatched[i]];
				}
				score1-=score2;  // score1 includes score2 because the peak scores were summed into specToScore

				if(score1+score2>result.score1+result.score2) {
//cerr << "Best solution changed to ["<<shift-shiftsOffset<<","<<t-shiftsOffset<<"], scores ("<<score1<<","<<score2<<")\n";
					result.shift1=(shift-shiftsOffset);   result.shift2=(t-shiftsOffset);
					result.score1=score1;                 result.score2=score2;
				}
			}
	}

	// Dispose of temporary memory structures
	for(shift=0; shift<preShiftPeaks.size(); shift++) {
		if(preShiftPeaks[shift][0]!=0) delete preShiftPeaks[shift][0];
		if(preShiftPeaks[shift][1]!=0) delete preShiftPeaks[shift][1];
	}
}


//
// getPairAlignsPAext - Extend consensus multiple alignments by aligning to single PRM spectra
//
void getPairAlignsPAext(Clusters &clusters, SpecSet &specSet,
						int startIdx, int endIdx, short aaDiff, float pmTol, float peakTol, float minMatchRatio,
						list<Results_PA> &results, list<TwoValues<float> > &ratios,
						vector<TwoValues<float> > &meansClst, vector<TwoValues<float> > &meansSpecs){
	int specIdx, clstIdx;
	Results_PA curResult;   TwoValues<float> curRatios;
	float totalConsensusScore;  // Variables used in checking minMatchRatio
	vector<float> totalSpecScores;
	clock_t curTime, totTime=0;  int alignsCount=0;

	startIdx=max(0,startIdx);   endIdx=min((unsigned int)endIdx,specSet.size()-1);   totalSpecScores.resize(endIdx+1);
//	results.resize((endIdx-startIdx+1)*clusters.size());
	for(specIdx=startIdx; specIdx<=endIdx; specIdx++) { totalSpecScores[specIdx]=0;
		for(int p=0; p<specSet[specIdx].size(); p++) totalSpecScores[specIdx]+=specSet[specIdx][p][1];
//cerr<<"totalSpecScores["<<specIdx<<"] = "<<totalSpecScores[specIdx]<<endl;
	}

	meansClst.resize(clusters.size());   for(int i=0; i<meansClst.size(); i++) meansClst[i].set(0,0);
	meansSpecs.resize(specSet.size());   for(int i=0; i<meansSpecs.size(); i++) meansSpecs[i].set(0,0);
	Spectrum cSpec;   Spectrum cSpecRev;
	for(clstIdx=0; clstIdx<clusters.size(); clstIdx++) {
		cout << "Processing cluster " << clstIdx+1 << "/" << clusters.size() << "... "; curTime=clock(); cout.flush();
		clusters.getSpecIfB(clstIdx,cSpec);	    clusters.getSpecIfY(clstIdx,cSpecRev);
		totalConsensusScore = 0;
		for(int p=0; p<cSpec.size(); p++) totalConsensusScore+=cSpec[p][1];
		for(int p=0; p<cSpecRev.size(); p++) totalConsensusScore+=cSpecRev[p][1];
//cerr<<"totalConsensusScore = "<<totalConsensusScore<<endl;
		for(specIdx=startIdx; specIdx<=endIdx; specIdx++) {
			curResult.spec1=specIdx;
			curResult.spec2=clstIdx;

			getPairAlignsPAext_aux(specSet[specIdx],cSpec,cSpecRev,aaDiff,pmTol,peakTol,curResult);

			meansSpecs[specIdx][0]=(meansSpecs[specIdx][0]*meansSpecs[specIdx][1]+curResult.score1)/(meansSpecs[specIdx][1]+1);
			meansSpecs[specIdx][1]++;
			meansClst[clstIdx][0]=(meansClst[clstIdx][0]*meansClst[clstIdx][1]+curResult.score2)/(meansClst[clstIdx][1]+1);
			meansClst[clstIdx][1]++;

			if(minMatchRatio>0) {
				curRatios[0] = curResult.score1/totalSpecScores[specIdx];
				curRatios[1] = curResult.score2/totalConsensusScore;
//cerr<<"Ratios: " << curRatios[0] << ", "<< curRatios[1] <<endl;
				if(curRatios[0]>=minMatchRatio && curRatios[1]>=minMatchRatio)
					{ results.push_back(curResult); ratios.push_back(curRatios); }
			} else { results.push_back(curResult); ratios.push_back(curRatios); }
		}

		curTime=clock()-curTime; totTime+=curTime;
		cout<< curTime/CLOCKS_PER_SEC << " secs, number of new pair aligns = "<<results.size()-alignsCount<<". Average time = "<<totTime/((clstIdx+1)*CLOCKS_PER_SEC)<<", ETA = "<< (clusters.size()-clstIdx-1)/(clstIdx+1)*(totTime/CLOCKS_PER_SEC) << endl; cout.flush();
		alignsCount = results.size();
	}

}


bool Load_pabatch(const char *filename, vector<PairAlign> &aligns) {
	ifstream input(filename);
	if (!input) { cerr << "ERROR: cannot open " << filename << "\n";   aligns.resize(0);  return false; }
	int numPairs, numShifts, i, j;

	input >> numPairs;    aligns.resize(numPairs);
	for(i=0; i<numPairs ; i++) {
		input >> aligns[i].spec1 >> aligns[i].spec2 >> numShifts;
		aligns[i].shifts.resize(numShifts);   aligns[i].results.resize(numShifts);
		for(j=0 ; j<numShifts ; j++)  input >> aligns[i].shifts[j];
	}

	input.close();
	return true;
}

bool Load_results(const char *filename, vector<Results_ASP> &results) {
	ifstream input(filename);
	if (!input) { cerr << "ERROR: cannot open " << filename << "\n";   results.resize(0);  return false; }
	int numResults, i;
	char separator;

	input >> numResults;    results.resize(numResults);
	for(i=0; i<numResults ; i++) {
		input >> results[i].spec1 >> separator >> results[i].spec2 >> separator >> results[i].shift1 >> separator >> results[i].score1 >> separator >> results[i].score2;
		results[i].spec1--; results[i].spec2--;  // Convert indices to c++ values
	}

	input.close();
	return true;
}

bool Load_results(const char *filename, vector<Results_PA> &results) {
	ifstream input(filename);
	if (!input) { cerr << "ERROR: cannot open " << filename << "\n";   results.resize(0);  return false; }
	int numResults, i;
	char separator;

	input >> numResults;    results.resize(numResults);
	for(i=0; i<numResults ; i++) {
		input >> results[i].spec1 >> separator >> results[i].spec2 >> separator >> results[i].shift1 >> separator >> results[i].shift2 >> separator >> results[i].score1 >> separator >> results[i].score2;
		results[i].spec1--; results[i].spec2--;  // Convert indices to c++ values
	}

	input.close();
	return true;
}

/*
bool Load_resultsASPbin(char *filename, vector<Results_ASP> &results) {
	FILE *fp;   unsigned int numEntries;
	fp = fopen(filename,"r");
	if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return false; }

	fread(&numEntries,sizeof(unsigned int),1,fp);  // Number of entries in the file
	results.resize(numEntries);

	float data[5];   int read_result;
	for(unsigned int resIdx=0; resIdx<numEntries; resIdx++) {
		read_result = fread(data,sizeof(float),5,fp);   if(read_result!=5) { cerr<<"ERROR reading "<<filename<<"!\n"; results.resize(0); return false; }
		results[resIdx].spec1=(int)data[0]-1;   results[resIdx].spec2=(int)data[1]-1;  // +1 for matlab arrays starting at 1
		results[resIdx].shift1=data[2];    results[resIdx].score1=data[3];   results[resIdx].score2=data[4];
	}

	fclose(fp);
	return true;
}

bool Load_resultsPAbin(char *filename, vector<Results_PA> &results) {
	FILE *fp;   unsigned int numEntries;
	fp = fopen(filename,"r");
	if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return false; }

	fread(&numEntries,sizeof(unsigned int),1,fp);  // Number of entries in the file
	results.resize(numEntries);

	float data[6];   int read_result;
	for(unsigned int resIdx=0; resIdx<numEntries; resIdx++) {
		read_result = fread(data,sizeof(float),6,fp);   if(read_result!=6) { cerr<<"ERROR reading "<<filename<<"!\n"; results.resize(0); return false; }
		results[resIdx].spec1=(int)data[0]-1;   results[resIdx].spec2=(int)data[1]-1;  // +1 for matlab arrays starting at 1
		results[resIdx].shift1=data[2];   results[resIdx].shift2=data[3];
		results[resIdx].score1=data[4];   results[resIdx].score2=data[5];
	}

	fclose(fp);
	return true;
}
*/

bool Load_resultsPA(const char *filename, vector<Results_PA> &results) {
	ifstream input(filename);
	if (!input) { cerr << "ERROR: cannot open " << filename << "\n";   results.resize(0);  return false; }
	int numResults, i, j;
	char separator;

	input >> numResults;    results.resize(numResults);
	for(i=0; i<numResults ; i++) {
		input >> results[i].spec1 >> separator >> results[i].spec2 >> separator >> results[i].shift1 >> separator >> results[i].shift2 >> separator >> results[i].score1 >> separator >> results[i].score2;
		results[i].spec1--; results[i].spec2--;  // Convert indices to c++ values
	}

	input.close();
	return true;
}

template<class T> bool Load_results_bin_multiple(const char *filename, vector<T> &results) {
	unsigned int totalCount=0, numPairs, fIdx, pIdx;
	FILE *fp;
	BufferedLineReader blr;
	if(blr.Load(filename)<=0) { cerr<<"ERROR reading "<<filename<<"!\n"; results.resize(0); return false; }

	for(fIdx=0; fIdx<blr.size(); fIdx++) {
		if(strlen(blr.getline(fIdx))==0) continue;
		fp = fopen(blr.getline(fIdx),"r");  if (fp==0) { cerr << "ERROR: cannot open " << blr.getline(fIdx) << "\n";   return false; }
		if (fread(&numPairs,sizeof(unsigned int),1,fp) != 1 && !feof(fp)) {
			cerr<<"ERROR reading "<<filename<<"!\n"; results.resize(0);
			return false;
		}  // Number of pairs in the file
		fclose(fp);
		totalCount+=numPairs;
	}

	results.resize(totalCount);   totalCount=0;
	vector<T> curPairs;
	for(fIdx=0; fIdx<blr.size(); fIdx++) {
		curPairs.resize(0);
		if (strlen(blr.getline(fIdx))>0 and not Load_results_bin(blr.getline(fIdx),curPairs)) { cerr<<"ERROR reading "<<blr.getline(fIdx)<<"!\n"; results.resize(0); return false; }
		for(pIdx=0; pIdx<curPairs.size(); pIdx++) results[totalCount+pIdx]=curPairs[pIdx];
		totalCount+=curPairs.size();
	}

	return true;
}

void instantiator(){
	const char *cP;
	vector<Results_SP> rSP; vector<Results_ASP> rASP; vector<Results_PA> rPA;
	list<Results_SP> lSP; list<Results_ASP> lASP; list<Results_PA> lPA;

	Load_results_bin(cP,rSP); Load_results_bin(cP,rASP); Load_results_bin(cP,rPA);
	Load_results_bin_multiple(cP,rSP); Load_results_bin_multiple(cP,rASP); Load_results_bin_multiple(cP,rPA);
	Save_results_bin(cP,0,rSP.begin());   Save_results_bin(cP,0,lSP.begin());
	Save_results_bin(cP,0,rASP.begin());  Save_results_bin(cP,0,lASP.begin());
	Save_results_bin(cP,0,rPA.begin());   Save_results_bin(cP,0,lPA.begin());
}

template<class T> bool Load_results_bin(const char *filename, vector<T> &results) {
	FILE *fp;   unsigned int numEntries, vCount, resIdx, vIdx;
	fp = fopen(filename,"r");
	if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return false; }
	if (fread(&numEntries,sizeof(unsigned int),1,fp) != 1 && !feof(fp)) {
		cerr << "ERROR: cannot open " << filename << "\n";
		return false;
	}  // Number of entries in the file
	results.resize(numEntries); if(numEntries==0) { fclose(fp); return true; }
	vCount=results[0].loadSz();
	float *data = (float *)malloc(vCount*sizeof(float));   int read_result;
	vector<float> dataV(vCount);
	for(unsigned int resIdx=0; resIdx<numEntries; resIdx++) {
		read_result = fread(data,sizeof(float),vCount,fp);   if(read_result!=vCount) { cerr<<"ERROR reading "<<filename<<"!\n"; results.resize(0); return false; }
		for(vIdx=0;vIdx<2;vIdx++) dataV[vIdx]=data[vIdx]-1;   // Spectrum indices are 1-based
		for(vIdx=2;vIdx<vCount;vIdx++) dataV[vIdx]=data[vIdx];
		results[resIdx].load(dataV);
	}
	free(data); fclose(fp); return true;
}

/* commented out for some other reason..
bool Load_resultsASPbin_multiple(char *filename, vector<Results_ASP> &results) {
	unsigned int totalCount=0, numPairs, fIdx, pIdx;
	FILE *fp;
	BufferedLineReader blr;
	if(blr.Load(filename)<=0) { cerr<<"ERROR reading "<<filename<<"!\n"; results.resize(0); return false; }

	for(fIdx=0; fIdx<blr.size(); fIdx++) {
		if(strlen(blr.getline(fIdx))==0) continue;
		fp = fopen(blr.getline(fIdx),"r");  if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return false; }
		fread(&numPairs,sizeof(unsigned int),1,fp);  // Number of pairs in the file
		fclose(fp);
		totalCount+=numPairs;
	}

	results.resize(totalCount);   totalCount=0;
	vector<Results_ASP> curPairs;
	for(fIdx=0; fIdx<blr.size(); fIdx++) {
		curPairs.resize(0);
		if (strlen(blr.getline(fIdx))>0 and not Load_resultsASPbin(blr.getline(fIdx),curPairs)) { cerr<<"ERROR reading "<<blr.getline(fIdx)<<"!\n"; results.resize(0); return false; }
		for(pIdx=0; pIdx<curPairs.size(); pIdx++) results[totalCount+pIdx]=curPairs[pIdx];
		totalCount+=curPairs.size();
	}

	return true;
}

//
//  Save_resultsPAbin - saves results_PA as blocks of 6 floats, plus one entry count at the start
//
bool Save_resultsPAbin(char *filename, list<Results_PA> &results) {
	FILE *fp;   unsigned int numEntries = results.size();
	float data[6];
	unsigned int i,p;

	fp = fopen(filename,"w");
	if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return false; }

	fwrite(&numEntries,sizeof(int),1,fp);  // Number of entries in the file

	for(list<Results_PA>::iterator iter= results.begin();iter!=results.end();iter++) {
		data[0]=(*iter).spec1+1;  data[1]=(*iter).spec2+1;    // +1 for matlab arrays starting at 1
		data[2]=(*iter).shift1;   data[3]=(*iter).shift2;
		data[4]=(*iter).score1;   data[5]=(*iter).score2;
		fwrite(data,sizeof(float),6,fp);
	}

	fclose(fp);
	return true;
}

//
//  Save_resultsASPbin - saves results ASP as blocks of 5 floats, plus one entry count at the start
//
bool Save_resultsASPbin(char *filename, list<Results_ASP> &results) {
	FILE *fp;   unsigned int numEntries = results.size();
	float data[5];
	unsigned int i,p;

	fp = fopen(filename,"w");
	if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return false; }

	fwrite(&numEntries,sizeof(unsigned int),1,fp);  // Number of entries in the file

	for(list<Results_ASP>::iterator iter= results.begin();iter!=results.end();iter++) {
		data[0]=(*iter).spec1+1;   data[1]=(*iter).spec2+1;  // +1 for matlab arrays starting at 1
		data[2]=(*iter).shift1;    data[3]=(*iter).score1;   data[4]=(*iter).score2;
		fwrite(data,sizeof(float),5,fp);
	}

	fclose(fp);
	return true;
}


bool Save_resultsASPbin(char *filename, vector<Results_ASP> &results) {
	FILE *fp;   unsigned int numEntries = results.size();
	float data[5];
	unsigned int i,p;

	fp = fopen(filename,"w");
	if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return false; }

	fwrite(&numEntries,sizeof(unsigned int),1,fp);  // Number of entries in the file

	for(i=0;i<numEntries;i++) {
		data[0]=results[i].spec1+1;  data[1]=results[i].spec2+1;   // +1 for matlab arrays starting at 1
		data[2]=results[i].shift1;   data[3]=results[i].score1;    data[4]=results[i].score2;
		fwrite(data,sizeof(float),5,fp);
	}

	fclose(fp);
	return true;
}
*/

//
// Save_results - save the alignment results in text format
// Line 1: Number_of_pairs
// Per pair, first line:       spec1 spec2 number_of_shift_lines
// Per pair, following lines:  shift_value shift_score; space separated list of matched PRMs in spec1; space separated list of matched PRMs in spec2
//
bool Save_results(const char *filename, vector<PairAlign> &aligns) {
	ofstream output(filename);
	if (!output) { cerr << "ERROR: cannot open " << filename << "\n";   aligns.resize(0);  return false; }
	unsigned int i,j,k;

	output << aligns.size() << endl;
	for(i=0; i<aligns.size() ; i++) {
		output << aligns[i].spec1 << " " << aligns[i].spec2 << " " << aligns[i].shifts.size() << endl;
		for(j=0; j<aligns[i].shifts.size() ; j++){
			output << aligns[i].shifts[j] << " " << aligns[i].results[j].score << ";";
			for (k=0; k<aligns[i].results[j].idxMatched1.size() ; k++) output << " " << aligns[i].results[j].idxMatched1[k]+1;  // +1 for matlab arrays starting at 1
			output << ";";
			for (k=0; k<aligns[i].results[j].idxMatched2.size() ; k++) output << " " << aligns[i].results[j].idxMatched2[k]+1;
			output << endl;
		}
	}

	output.close();    return true;
}

bool combineResultsOCC(const char* filename, const char* outfile, SpecSet& spectra, vector<int>& charges, vector<int>& matchedPeaks, vector<float>& minRatio) {
	list<Results_OCC> aligns;
	if (! Load_resultsOCC(filename, aligns)) {cerr << "ERROR: cannot open " << filename << "\n"; return false;}
	list<Results_OCC> outAligns;
	Results_OCC res;
	list<Results_OCC>::iterator alignIt;
	set<int> used;
	int index;
	for (int crit = 0; crit < charges.size(); crit++) {
		index = 0;
		int charge = charges[crit]; int mPeaks = matchedPeaks[crit]; float minR = minRatio[crit] - 0.00001;
		for (alignIt = aligns.begin(); alignIt != aligns.end(); alignIt++) {
			if (used.count(index) > 0) {index ++; continue;}
			res = *alignIt;
			if (spectra[res.spec1].parentCharge >= charge && res.matchedPeaks >= mPeaks && min(res.score1, res.score2) > minR) {
				outAligns.push_back(res); used.insert(index);
			}
			index ++;
		}
	}
	return Save_resultsOCC(outfile, outAligns);
}

bool Save_resultsOCC(const char* filename, list<Results_OCC> &aligns) {
	ofstream output(filename);
	list<Results_OCC>::iterator it;
	if (!output) { cerr << "ERROR: cannot open " << filename << "\n";   aligns.resize(0);  return false; }

	output << aligns.size() << endl;
	for(it = aligns.begin(); it != aligns.end() ; it++) (*it).output(output, ' ');

	output.close();    return true;
}

bool Load_resultsOCC(const char *filename, list<Results_OCC> &aligns) {
	BufferedLineReader read;
	char* next;
	char* result;
	Results_OCC res;
	vector<float> v(7);
	int x;
	if (! read.Load(filename)) return false;
	int end = getInt(read.getline(0));
	aligns.clear();
	for (unsigned int i = 1; i <= end; i++) {
		next = read.getline(i);
		result = strtok( next, " " );
		x = 0;
		while( result != NULL ) {
			const char* floatStr = result;
			v[x] = getFloat(floatStr);
			result = strtok( NULL, " " );
			x ++;
		}
		res.load(v);
		aligns.push_back(res);
	}
	read.reset();
	return true;
}

//
// Save_results - save the alignment results in text format
// Line 1: Number_of_pairs
// Per pair, first line:       spec1 spec2 number_of_shift_lines
// Per pair, following lines:  shift_value shift_score; space separated list of matched PRMs in spec1; space separated list of matched PRMs in spec2
//
bool Save_resultsCS(const char *filename, list<Results_CS> &aligns, const char sep) {
	ofstream output(filename);
	if (!output) { cerr << "ERROR: cannot open " << filename << "\n"; return false; }
	unsigned int i,j,k;
	list<Results_CS>::iterator iter = aligns.begin();

	output << aligns.size() << endl;
	for(iter = aligns.begin(); iter != aligns.end() ; iter++)
		output << iter->spec1 << sep << iter->spec2 << sep << iter->specC << sep << iter->shift1 << sep << iter->shift2 << endl;

	output.close();    return true;
}


// findConnectorSpectra - Connect contigs by finding single PRM spectra to
//     connect pairs of contigs (consensus spectra from multiple alignments)
//
// contigs - Input contigs
// specSet - Input PRM spectra (possible connectors)
// minPeakDist  - Minimum distance between matched peaks (in the merged version of the overlapping 3 spectra)
// mergedPeakCount - A peak is reported in the merged spectrum if it occurs in at least mergedPeakCount spectra (1,2 or 3).
// results - all significant connector triples
// ratios  - set of triples with the matched score percentages for contig1/contig2/connector spectrum
//
void findConnectorSpectra(Clusters &contigs, SpecSet &specSet,
							float peakTol, float minRatio, float minPeakDist, int mergedPeakCount,
							list<Results_CS> &results, list<vector<float> > &ratios) {
	unsigned int numContigs = contigs.size(), numContigPairs = (numContigs*(numContigs-1))/2;
	vector<TwoValues<int> > contigPairs(numContigPairs);            // Indices of paired contigs
	vector<TwoValues<vector<float> > > contigShifts(numContigPairs);  // Possible shifts between paired contigs (pos.0 contig2 unchanged, pos.1 contig2 reversed)
	vector<TwoValues<vector<float> > > connectorShifts(numContigs);   // Possible shifts between all contigs and a given connector (pos.0 contig unchanged, pos.1 contig reversed)
	vector<TwoValues<float> > minMaxScores_ct(contigs.size());        // Minimum scores to meet minRatioCt criterion (pos.0), total spectrum score (pos.1)
	Spectrum contig1, contig2;
	unsigned int pairIdx, cIdx1, cIdx2, csIdx, peakIdx;

cerr<<"WARNING: batch.cpp::findConnectorSpectra() was not revised to handle high precision data!\n";

	// Compute lists of possible shifts between all contigs
	pairIdx=0;
	for(cIdx1=0; cIdx1<numContigs; cIdx1++) {
		contigs.getSpecIfB(cIdx1,contig1);
		minMaxScores_ct[cIdx1][1]=0; for(peakIdx=0; peakIdx<contig1.size(); peakIdx++) minMaxScores_ct[cIdx1][1]+=contig1[peakIdx][1]; minMaxScores_ct[cIdx1][0]=minMaxScores_ct[cIdx1][1]*minRatio;
		for(cIdx2=cIdx1+1; cIdx2<numContigs; cIdx2++) {
			contigPairs[pairIdx].set(cIdx1,cIdx2);

			contigs.getSpecIfB(cIdx2,contig2);
			computeShifts2(contig1,contig2,contigShifts[pairIdx][0]);

			// TODO: Reverse contigs2 and invoke computeShifts2 again
			contigShifts[pairIdx][1].resize(0);

			pairIdx++;
		}
	}

	// Variables used below
	Spectrum mergedSpec;    // Projection of the 3 aligned spectra
	vector<vector<int> > peaksIndices;  // Indices of the merged peaks in each of the original 3 spectra
										  //   peaksIndices[i][j] - index of peak i in spectrum j. -1 if none.
	vector<vector<int> > peakNeighs;    // Peak neighbors in mergedSpec according to minPeakDist
	vector<int> peakMatches;            // Maximal sparse subset of peaks from mergedSpec (peak indices)
	vector<float> matchScores(3);         // Match scores per spectrum: contig1 (pos.0), contig2 (pos.1), connector spectrum (pos.2)
	Results_CS tmpResult;                 // Temporary object used to insert into results
	float bestScore, curScore;
	unsigned int bestShift12Idx, bestShift1csIdx;
	vector<float> bestMatchScores(3);     // Match scores for the best shifts

	// Iterate over all possible connectors
	for(csIdx=0; csIdx<specSet.size(); csIdx++) {
		TwoValues<float> minMaxMatchScoreCs(0,0); for(peakIdx=0; peakIdx<specSet[csIdx].size(); peakIdx++) minMaxMatchScoreCs[1]+=specSet[csIdx][peakIdx][1]; minMaxMatchScoreCs[0]=minMaxMatchScoreCs[1]*minRatio;
		// Compute shifts between connector and all contigs
		for(cIdx1=0; cIdx1<numContigs; cIdx1++) {
			contigs.getSpecIfB(cIdx1,contig1);
			computeShifts2(contig1,specSet[csIdx],connectorShifts[cIdx1][0]);
// TODO: Reverse the contig and invoke computeShifts2 again
			connectorShifts[cIdx1][1].resize(0);
		}

		// Try this connector spectrum for all possible pairs of contigs
		for(pairIdx=0; pairIdx<numContigPairs; pairIdx++) {
			cIdx1 = contigPairs[pairIdx][0];    cIdx2 = contigPairs[pairIdx][1];
			contigs.getSpecIfB(cIdx1,contig1);  contigs.getSpecIfB(cIdx2,contig2);
			vector<float> &shifts12 = contigShifts[pairIdx][0];
			vector<float> &shifts1cs = connectorShifts[cIdx1][0];
			vector<bool> shifts2cs_ok;
			int offset = Utils::list_to_vector_bool(connectorShifts[cIdx2][0], shifts2cs_ok, peakTol, 10);
//cerr << "Triple ("<<cIdx1<<","<<cIdx2<<","<<csIdx<<"):\n";

			// Find pairs of shifts for which there is at least a pair of matching peaks on every aligned contig/spectrum
			int eligible=0; bestScore=0;
			for(unsigned int shift12Idx=0; shift12Idx<shifts12.size(); shift12Idx++)
				for(unsigned int shifts1csIdx=0; shifts1csIdx<shifts1cs.size(); shifts1csIdx++) {
					float shift2cs = shifts1cs[shifts1csIdx] - shifts12[shift12Idx];
					int shift2csInt = offset+(int)round(10*shift2cs);
					if(!shifts2cs_ok[shift2csInt]) continue;

					MergeSpectra3(contig1,contig2,specSet[csIdx],shifts12[shift12Idx],shifts1cs[shifts1csIdx],peakTol,mergedPeakCount,mergedSpec,peaksIndices);
					if(mergedSpec.size()==0) continue;
					eligible++;
					getNeighborsL(mergedSpec, minPeakDist, peakTol, peakNeighs);
					curScore = denovo_LtoR(mergedSpec, peakNeighs, peakMatches);

					if(curScore>bestScore) {
						matchScores[0]=0;matchScores[1]=0;matchScores[2]=0;
						for(peakIdx=0; peakIdx<peakMatches.size(); peakIdx++) {
							if(peaksIndices[peakMatches[peakIdx]][0]>=0) matchScores[0]+=contig1[peaksIndices[peakMatches[peakIdx]][0]][1];
							if(peaksIndices[peakMatches[peakIdx]][1]>=0) matchScores[1]+=contig2[peaksIndices[peakMatches[peakIdx]][1]][1];
							if(peaksIndices[peakMatches[peakIdx]][2]>=0) matchScores[2]+=specSet[csIdx][peaksIndices[peakMatches[peakIdx]][2]][1];
						}

						if(matchScores[0]>=minMaxScores_ct[cIdx1][0] and matchScores[1]>=minMaxScores_ct[cIdx2][0] and matchScores[2]>=minMaxMatchScoreCs[0]) {
//cerr<<"mergedSpec\n"; mergedSpec.output(cerr);
//cerr << " -- eligible shifts ["<<shifts12[shift12Idx]<<","<<shifts1cs[shifts1csIdx]<<","<<shift2cs<<"]\n";
//cerr << "Matched peaks for score "<<curScore<<":\n";
//for(unsigned int i=0; i<peakMatches.size(); i++) cerr<<peakMatches[i]<<": "<<mergedSpec[peakMatches[i]][0]<<"\t"<<mergedSpec[peakMatches[i]][1]<<endl;
							bestShift12Idx=shift12Idx; bestShift1csIdx=shifts1csIdx; bestScore = curScore;
							bestMatchScores[0]=matchScores[0]; bestMatchScores[1]=matchScores[1]; bestMatchScores[2]=matchScores[2];
						}
					}

					// Compute match ratios and reject if any is below the given threshold
				}

			if(bestScore>0) {
//cerr << " -- best match at ["<<shifts12[bestShift12Idx]<<","<<shifts1cs[bestShift1csIdx]<<"], score "<<bestScore<<" (out of "<<eligible<<" eligible shifts computed).\n";
				bestMatchScores[0]/=minMaxScores_ct[cIdx1][1];  // Convert match scores to match ratios
				bestMatchScores[1]/=minMaxScores_ct[cIdx2][1];
				bestMatchScores[2]/=minMaxMatchScoreCs[1];
				ratios.push_back(bestMatchScores);

				tmpResult.spec1 = cIdx1;  tmpResult.spec2 = cIdx2;   tmpResult.specC = csIdx;
				tmpResult.shift1 = shifts12[bestShift12Idx];    tmpResult.score1 = bestMatchScores[0];
				tmpResult.shift2 = shifts1cs[bestShift1csIdx];  tmpResult.score2 = bestMatchScores[1];
				tmpResult.spec2rev = false;
				results.push_back(tmpResult);
			}
		}
	}
}
