#ifndef SPECTRAL_PAIRS_H
#define SPECTRAL_PAIRS_H

#include "spectrum.h"
#include "batch.h"
#include "label.h"
#include "msn.h"
#include "dekel_align.h"
#include "alignment_scoring.h"
#include <fstream>

/**
 *@deprecated use SplitSpectra/SplitPairs instead
 */
void SplitPairs(SpecSet &specSet, vector<Results_ASP> &aligns, float peakTol,
		int maxAAjump, float penalty_sameVert, float penalty_ptm,
		bool forceSymmetry, SpecSet &specSetNew,
		vector<Results_ASP> &alignsNew,
		vector<vector<TwoValues<int> > > &matches, vector<bool> &pairFlipped,
		vector<vector<float> > *dbg_matchScores = 0, ofstream *debug = 0);

/**
 * Decides on a consensus orientation for  a set of edges. Flips spectra
 * according to assigned orientation.
 *
 *@param specSet - star spectra
 *@param aligns  - pairwise alignments
 *@param peakTol - peak mass tolerance
 *@param pmTol   - parent mass tolerance
 *@param maxAAjump - maximum amino-acid mass jump. Consecutive matched peaks with mass
 *                     difference in this range must match an amino acid mass.
 *@param penalty_sameVert
 *@param penalty_ptm
 *@param matches     - position [i] contains set of peak indices for matched peaks in aligns[i].
 *@param specFlipped - indicates whether spectra were reversed when added to the ABruijn graph
 *@param modPos      - location of the mass difference when the spectrum was added to the ABruijn graph
 *@param minMatchedPeaks - Minimum number of matched peaks to keep matches (ABruijn glues) between two spectra
 *@param minEdgesToComponent - Minimum number of edges to spectra in the component for additional spectra to be added
 *@param forceSymmetry
 *@param alignStats  - position [i] %matched score in spec1/spec2 (pos 0/1) and #matched peaks (pos 2) for aligns[i]
 *@param labelsP
 *@return
 */
template<class T> void SplitPairs(SpecSet &specSet, vector<T> &aligns,
		float peakTol, float pmTol, int maxAAjump, float maxModMass, float penalty_sameVert,
		float penalty_ptm, vector<vector<TwoValues<int> > > &matches,
		vector<bool> &specFlipped, vector<float> &modPos,
		unsigned int minMatchedPeaks=0, unsigned int minEdgesToComponent=0,
		bool forceSymmetry = false, vector<vector<float> > *alignStats = NULL,
		vector<SpectrumPeakLabels> *labelsP = NULL) {
	vector<list<int> > alignsEntries(specSet.size());         // List of pairs in aligns that incide on the corresponding vertex
	vector<float> flipScores(specSet.size());                 // Accumulated don't_flip/flip score balance per vertex
	SpecSet specSetRev; specSetRev.resize(specSet.size());    // Holds reversed versions of the spectra in the component
	vector<int> numPairsPerSpectrum(specSet.size());		  // Number of pairs to from spectrum i to spectra included in the component
	list<int> toProcess;                                      // Indices of spectra left to process
	TwoValues<vector<vector<TwoValues<int> > > > tmpMatches;  // Temporary matches when a neighbor spectrum S1 in the i-th pair is matched directly/reversed (its pair S2 is assumed to already be in the component)
	vector<TwoValues<float> > tmpModPos(aligns.size());       // Temporary storage for the position of the modification
	vector<bool> matchComputed(aligns.size());                // Indicates whether the i-th entry in tmpAligns was already computed
	Results_PA curPair;                                       // Used to access aligns for template class T = Results_ASP or Results_PA
	bool isASP;
	unsigned int specsInComponent=0;					      // Number of spectra already included in the component

	// Initializations
	matches.resize(aligns.size());       modPos.resize(aligns.size());
	tmpMatches[0].resize(aligns.size()); tmpMatches[1].resize(aligns.size());
	for(unsigned int i=0;i<aligns.size();i++) {
		if(specSet[aligns[i].spec1].size()==0 or specSet[aligns[i].spec1].size()!=specSetRev[aligns[i].spec1].size())
		{ specSet[aligns[i].spec1].reverse(0, &specSetRev[aligns[i].spec1]); toProcess.push_back(aligns[i].spec1); }
		if(specSet[aligns[i].spec2].size()==0 or specSet[aligns[i].spec2].size()!=specSetRev[aligns[i].spec2].size())
		{ specSet[aligns[i].spec2].reverse(0, &specSetRev[aligns[i].spec2]); toProcess.push_back(aligns[i].spec2); }
		matchComputed[i] = false;
		alignsEntries[aligns[i].spec1].push_back(i);   alignsEntries[aligns[i].spec2].push_back(i);
		matches[i].resize(0);
	}

	list<int>::iterator pIter;
	for(pIter=toProcess.begin(); pIter!=toProcess.end(); pIter++) flipScores[*pIter]=0;

	// Use edge scores as percentages of spectrum score to select the edge/vertex to initialize the direction selection
	vector<float> specScores(specSet.size());
	float bestScore=0, curScore; int bestScoreIdx=-1, specToAdd=-1;  // specToAdd is index of the new spectrum being added to the component
	for(unsigned int i=0;i<specSet.size(); i++) {
		specScores[i]=0;   numPairsPerSpectrum[i]=0;
		for(unsigned int j=0; j<specSet[i].size(); j++) specScores[i]+=specSet[i][j][1];
	}
	// Select starting spectrum by participation in highest-scoring alignment
	for(unsigned int i=0;i<aligns.size(); i++) {
		curScore=min(aligns[i].score1/specScores[aligns[i].spec1],aligns[i].score2/specScores[aligns[i].spec2]);
		numPairsPerSpectrum[aligns[i].spec1]++;   numPairsPerSpectrum[aligns[i].spec2]++;
		if(bestScoreIdx<0 or bestScore<curScore) { bestScoreIdx=i; bestScore=curScore; specToAdd=aligns[bestScoreIdx].spec1; }
	}
	// Select starting spectrum by highest number of pairs
	pIter=toProcess.begin();
	specToAdd=*pIter;   bestScore=numPairsPerSpectrum[specToAdd];   pIter++;
	for(; pIter!=toProcess.end(); pIter++) {
		if(numPairsPerSpectrum[*pIter]>bestScore) { bestScore=numPairsPerSpectrum[*pIter]; specToAdd=*pIter; }
		numPairsPerSpectrum[*pIter]=0;
	}
//	cout<<"  - (SplitPairs) Component initialized with spectrum "<<specToAdd<<" with "<<bestScore<<" neighbors\n";

	flipScores[specToAdd]=-1;                  // Arbitrarily decide that this spectrum is not flipped
	pIter=toProcess.begin(); while(*pIter!=specToAdd) pIter++; toProcess.erase(pIter);
	while(specToAdd>=0) {  // Process the inclusion of specToAdd to the oriented component
		// Flip/don't-flip specToAdd
		int flipDir=0;  // Spectrum not flipped (default)
		if(flipScores[specToAdd]>0) {
			specSet[specToAdd]=specSetRev[specToAdd]; flipDir=1; specFlipped[specToAdd]=true;
			if(labelsP) (*labelsP)[specToAdd].reverse();
		}
		specsInComponent++;

		// Process matches between newly included specToAdd and its alignments to spectra not yet included in the oriented component
		for(pIter = alignsEntries[specToAdd].begin(); pIter!=alignsEntries[specToAdd].end(); pIter++) {
			curPair = aligns[*pIter];   isASP = (fabs(curPair.shift1)<=pmTol or fabs(curPair.shift2)<=pmTol);
			if(!matchComputed[*pIter]) {
				// otherSpec is not in the component yet. Compute how the addition of specToAdd changes otherSpec's orientation selection.
				int otherSpec, otherSpecPos;   // otherSpec=index of current paired spectrum, otherSpecPos=position of otherSpec index in aligns
				if(curPair.spec1==specToAdd) { otherSpec=curPair.spec2; otherSpecPos=1; }
				else { otherSpec=curPair.spec1; otherSpecPos=0; }

				TwoValues<float> matchScore1(0,0),matchScore2(0,0);
				vector<int> indices;
				vector<Spectrum> results(4);

				TwoValues<float> score1, score2;
				vector<int> matchA1, matchA2, matchB1, matchB2;
				bool gluesAccepted=false;  // A pair's glues are accepted if the pair matches the minimum number of peaks (depending on orientation of the paired spectra)
				int extraMatchPeaks=0;     // ASP pairs always match either/both PRMs at 0/18 or/and PM-19/PM-1 so these have higher numbers of peaks that need to match

//cerr<<"Spectrum "<<specToAdd<<":\n"; specSet[specToAdd].output(cerr);
//cerr<<"Spectrum "<<otherSpec<<":\n"; specSet[otherSpec].output(cerr);

				if(isASP) {
					if(fabs(curPair.shift1)<=peakTol) extraMatchPeaks++;  // matching 0/18 does not count towards achieving minMatchedPeaks
					if(fabs(curPair.shift2)<=peakTol) extraMatchPeaks++;  // matching PM-19/PM-1 does not count towards achieving minMatchedPeaks

					tmpModPos[*pIter][0] = dekel_align(&specSet[specToAdd],&specSet[otherSpec],peakTol,&results[0],&results[1],maxAAjump,penalty_sameVert,penalty_ptm,forceSymmetry,false);
/*if(specToAdd==9642) {
	cerr<<"  ==> ASP ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<results[0].size()<<" / "<<minMatchedPeaks+extraMatchPeaks<<"\n";
	if(otherSpec==9630)
		for(unsigned int i=0; i<results[0].size(); i++)
			cerr<<" ------ ("<<results[0][i][0]<<","<<results[0][i][1]<<")\n";
}*/
					if(results[0].size()<minMatchedPeaks+extraMatchPeaks) tmpMatches[0][*pIter].resize(0);
					else {
						for(unsigned int i=0; i<results[0].size(); i++) { matchScore1[0]+=results[0][i][1]; matchScore2[0]+=results[1][i][1]; }
						specSet[specToAdd].massesToIndices(results[0].peakList, indices, peakTol);
						tmpMatches[0][*pIter].resize(indices.size());
						for(unsigned int i=0; i<indices.size(); i++) tmpMatches[0][*pIter][i][1-otherSpecPos]=indices[i];
						specSet[otherSpec].massesToIndices(results[1].peakList, indices, peakTol);
						for(unsigned int i=0; i<indices.size(); i++) tmpMatches[0][*pIter][i][otherSpecPos]=indices[i];
						gluesAccepted = true;
					}
					//if(specToAdd==111 and otherSpec==246)
					//	{ cerr<<"Direct match:\n";
					//	for(unsigned int i=0;i<results[0].size();i++) cerr<<"  ("<<results[0][i][0]<<","<<results[0][i][1]<<")/("<<results[1][i][0]<<","<<results[1][i][1]<<")\n"; }

					tmpModPos[*pIter][1] = dekel_align(&specSet[specToAdd],&specSetRev[otherSpec],peakTol,&results[2],&results[3],maxAAjump,penalty_sameVert,penalty_ptm,forceSymmetry,false);
/*if(specToAdd==9642){
	cerr<<"  ==> ASP rev ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<results[2].size()<<" / "<<minMatchedPeaks+extraMatchPeaks<<"\n";
	if(otherSpec==9630)
		for(unsigned int i=0; i<results[0].size(); i++)
			cerr<<" ------ ("<<results[0][i][0]<<","<<results[0][i][1]<<")\n";
}*/
					if(results[2].size()<minMatchedPeaks+extraMatchPeaks) tmpMatches[1][*pIter].resize(0);
					else {
						for(unsigned int i=0; i<results[2].size(); i++) { matchScore1[1]+=results[2][i][1]; matchScore2[1]+=results[3][i][1]; }
						specSet[specToAdd].massesToIndices(results[2].peakList, indices, peakTol);
						tmpMatches[1][*pIter].resize(indices.size());
						for(unsigned int i=0; i<indices.size(); i++) tmpMatches[1][*pIter][i][1-otherSpecPos]=indices[i];
						specSetRev[otherSpec].massesToIndices(results[3].peakList, indices, peakTol);
						for(unsigned int i=0; i<indices.size(); i++) tmpMatches[1][*pIter][i][otherSpecPos]=indices[i];
						gluesAccepted = true;
					}
					//if(specToAdd==111 and otherSpec==246)
					//	{ cerr<<"Reverse match:\n";
					//	for(unsigned int i=0;i<results[2].size();i++) cerr<<"  ("<<results[2][i][0]<<","<<results[2][i][1]<<")/("<<results[3][i][0]<<","<<results[3][i][1]<<")\n"; }
				} else {
					// Select shift for direct orientation
					float shift = 0;
					if(specToAdd==curPair.spec1) shift=curPair.shift1; else shift=-curPair.shift1;
					ScoreOverlap6(specSet[specToAdd], specSet[otherSpec], shift, peakTol, matchA1, matchA2);
//if(specToAdd==9642)
//	cerr<<"  ==> PA A1 ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<matchA1.size()<<" / "<<minMatchedPeaks<<"\n";
					if(matchA1.size()<minMatchedPeaks) { matchA1.resize(0); matchA2.resize(0); } else gluesAccepted = true;
					score1.set(0,0); for(unsigned int i=0; i<matchA1.size();i++) { score1[0]+=specSet[specToAdd][matchA1[i]][1]; score1[1]+=specSet[otherSpec][matchA2[i]][1]; }
					if(specToAdd==curPair.spec1) shift=curPair.shift2; else shift=-curPair.shift2;
					ScoreOverlap6(specSet[specToAdd], specSet[otherSpec], shift, peakTol, matchB1, matchB2);
//if(specToAdd==9642)
//	cerr<<"  ==> PA B1 ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<matchB1.size()<<" / "<<minMatchedPeaks<<"\n";
					if(matchB1.size()<minMatchedPeaks) { matchB1.resize(0); matchB2.resize(0); } else gluesAccepted = true;
					score2.set(0,0); for(unsigned int i=0; i<matchB1.size();i++) { score2[0]+=specSet[specToAdd][matchB1[i]][1]; score2[1]+=specSet[otherSpec][matchB2[i]][1]; }

					if( (score1[0]/specScores[specToAdd]+score1[1]/specScores[otherSpec]) > (score2[0]/specScores[specToAdd]+score2[1]/specScores[otherSpec]) )
					{ matchScore1[0]=score1[0]; matchScore2[0]=score1[1]; tmpMatches[0][*pIter].resize(matchA1.size()); for(unsigned int i=0;i<matchA1.size(); i++) { tmpMatches[0][*pIter][i][1-otherSpecPos]=matchA1[i]; tmpMatches[0][*pIter][i][otherSpecPos]=matchA2[i];} }
					else { matchScore1[0]=score2[0]; matchScore2[0]=score2[1]; tmpMatches[0][*pIter].resize(matchB1.size()); for(unsigned int i=0;i<matchB1.size(); i++) { tmpMatches[0][*pIter][i][1-otherSpecPos]=matchB1[i]; tmpMatches[0][*pIter][i][otherSpecPos]=matchB2[i];} }

					// Select shift for reversed orientation
					if(specToAdd==curPair.spec1) shift=curPair.shift1; else shift=-curPair.shift1;
					ScoreOverlap6(specSet[specToAdd], specSetRev[otherSpec], shift, peakTol, matchA1, matchA2);
//if(specToAdd==9642)
//	cerr<<"  ==> PA A1 rev ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<matchA1.size()<<" / "<<minMatchedPeaks<<"\n";
					if(matchA1.size()<minMatchedPeaks) { matchA1.resize(0); matchA2.resize(0); } else gluesAccepted = true;
					score1.set(0,0); for(unsigned int i=0; i<matchA1.size();i++) { score1[0]+=specSet[specToAdd][matchA1[i]][1]; score1[1]+=specSetRev[otherSpec][matchA2[i]][1]; }
					if(specToAdd==curPair.spec1) shift=curPair.shift2; else shift=-curPair.shift2;
					ScoreOverlap6(specSet[specToAdd], specSetRev[otherSpec], shift, peakTol, matchB1, matchB2);
//if(specToAdd==9642)
//	cerr<<"  ==> PA B1 rev ("<<curPair.spec1<<","<<curPair.spec2<<") matches "<<matchB1.size()<<" / "<<minMatchedPeaks<<"\n";
					if(matchB1.size()<minMatchedPeaks) { matchB1.resize(0); matchB2.resize(0); } else gluesAccepted = true;
					score2.set(0,0); for(unsigned int i=0; i<matchB1.size();i++) { score2[0]+=specSet[specToAdd][matchB1[i]][1]; score2[1]+=specSetRev[otherSpec][matchB2[i]][1]; }

					if( (score1[0]/specScores[specToAdd]+score1[1]/specScores[otherSpec]) > (score2[0]/specScores[specToAdd]+score2[1]/specScores[otherSpec]) )
					{ matchScore1[1]=score1[0]; matchScore2[1]=score1[1]; tmpMatches[1][*pIter].resize(matchA1.size()); for(unsigned int i=0;i<matchA1.size(); i++) { tmpMatches[1][*pIter][i][1-otherSpecPos]=matchA1[i]; tmpMatches[1][*pIter][i][otherSpecPos]=matchA2[i];} }
					else { matchScore1[1]=score2[0]; matchScore2[1]=score2[1]; tmpMatches[1][*pIter].resize(matchB1.size()); for(unsigned int i=0;i<matchB1.size(); i++) { tmpMatches[1][*pIter][i][1-otherSpecPos]=matchB1[i]; tmpMatches[1][*pIter][i][otherSpecPos]=matchB2[i];} }
				}

//				flipScores[otherSpec] += (matchScore2[1]-matchScore2[0])/specScores[otherSpec];
				if(gluesAccepted) {
					flipScores[otherSpec] += (matchScore2[1]-matchScore2[0])/specScores[otherSpec] + (matchScore1[1]-matchScore1[0])/specScores[specToAdd];
					numPairsPerSpectrum[otherSpec]++;
				}
				matchComputed[*pIter] = true;
			} else {
				// Copy adequate entry (according to flipDir) from tmpMatch to matches
//cerr<<" * copying match "<<curPair.spec1<<" vs "<<curPair.spec2<<endl;
				matches[*pIter].resize(tmpMatches[flipDir][*pIter].size());
				for(unsigned int i=0;i<matches[*pIter].size();i++) matches[*pIter][i]=tmpMatches[flipDir][*pIter][i];
//				if(alignRatios and flipDir==1) { float tmp=(*alignRatios)[*pIter][0]; (*alignRatios)[*pIter][0]=(*alignRatios)[*pIter][1]; (*alignRatios)[*pIter][1]=tmp; }
				if(isASP) modPos[*pIter] = tmpModPos[*pIter][flipDir];
			}
		}

		if(toProcess.size()==0) specToAdd=-1; else {
//			list<int>::iterator bestNewSpec=toProcess.begin();  pIter=bestNewSpec;  pIter++;
//			while(pIter!=toProcess.end()) { if(fabs(flipScores[*pIter])>fabs(flipScores[*bestNewSpec])) bestNewSpec=pIter; pIter++; }
//			specToAdd = *bestNewSpec;   toProcess.erase(bestNewSpec);
			list<int>::iterator bestNewSpec=toProcess.begin();
			while(bestNewSpec!=toProcess.end()) {
				if(numPairsPerSpectrum[*bestNewSpec]>=min(minEdgesToComponent,specsInComponent)) break;
				bestNewSpec++;
			}
			if(bestNewSpec!=toProcess.end()) {
				pIter=bestNewSpec;  pIter++;
				while(pIter!=toProcess.end()) {
					if(numPairsPerSpectrum[*pIter]>=min(minEdgesToComponent,specsInComponent) and fabs(flipScores[*pIter])>fabs(flipScores[*bestNewSpec])) bestNewSpec=pIter;
					pIter++;
				}
				specToAdd = *bestNewSpec;   toProcess.erase(bestNewSpec);
			} else specToAdd = -1;
		}
	}
//	cout<<"  - (SplitPairs) Component includes "<<specsInComponent<<" spectra.\n";

	if(alignStats) {  // Cannot be computed correctly until all flipping decisions are made
		alignStats->resize(aligns.size());
		for(unsigned int pairIdx=0;pairIdx<aligns.size();pairIdx++) {
			float score1=0, score2=0;   int s1=aligns[pairIdx].spec1, s2=aligns[pairIdx].spec2;
			(*alignStats)[pairIdx].resize(3);
			for(unsigned int i=0; i<matches[pairIdx].size(); i++)
				{ score1+=specSet[s1][matches[pairIdx][i][0]][1];  score2+=specSet[s2][matches[pairIdx][i][1]][1]; }
			(*alignStats)[pairIdx][0] = score1/specScores[s1];
			(*alignStats)[pairIdx][1] = score2/specScores[s2];
			(*alignStats)[pairIdx][2] = matches[pairIdx].size();
		}
	}
}

/**
 * Decides on a consensus orientation for  a set of edges. Flips spectra
 * according to assigned orientation.
 *
 *@param specSet
 *@param aligns
 *@param peakTol
 *@param maxAAjump
 *@param penalty_sameVert
 *@param penalty_ptm
 *@param matches
 *@param specFlipped
 *@param modPos
 *@param forceSymmetry
 *@param labelsP
 */
void SplitPairs2(SpecSet &specSet, vector<Results_ASP> &aligns, float peakTol,
		int maxAAjump, float penalty_sameVert, float penalty_ptm, vector<
				vector<TwoValues<int> > > &matches, vector<bool> &specFlipped,
		vector<float> &modPos, bool forceSymmetry = false, vector<
				SpectrumPeakLabels> *labelsP = NULL);

/**
 * Like SplitPairs2 but also processes Partial Overlap alignments (alignsPA).
 *
 *@param specSet
 *@param aligns
 *@param alignsPA
 *@param peakTol
 *@param maxAAjump
 *@param penalty_sameVert
 *@param penalty_ptm
 *@param matches
 *@param matchesPA
 *@param specFlipped
 *@param modPos
 *@param forceSymmetry
 *@param labelsP
 *@param alignRatios
 *@param alignRatiosPA
 */
void SplitPairs3(SpecSet &specSet, vector<Results_ASP> &aligns, vector<
		Results_PA> &alignsPA, float peakTol, int maxAAjump,
		float penalty_sameVert, float penalty_ptm, vector<
				vector<TwoValues<int> > > &matches, vector<vector<
				TwoValues<int> > > &matchesPA, vector<bool> &specFlipped,
		vector<float> &modPos, bool forceSymmetry = false, vector<
				SpectrumPeakLabels> *labelsP = NULL,
		vector<TwoValues<float> > *alignRatios = NULL,
		vector<TwoValues<float> > *alignRatiosPA = NULL);

/**
 * TODO: add description
 *
 *@param specSet
 *@param specSetSplit
 */
void SplitSpectra(SpecSet &specSet, SpecSet &specSetSplit);

/**
 * TODO: add description
 *
 *@param specSetSplit
 *@param aligns
 *@param peakTol
 *@param maxAAjump
 *@param penalty_sameVert
 *@param penalty_ptm
 *@param forceSymmetry
 *@param alignsNew
 *@param matches
 *@param pairFlipped
 *@param dbg_matchScores
 *@param debug
 */
void SplitAligns(SpecSet &specSetSplit, vector<Results_ASP> &aligns,
		float peakTol, int maxAAjump, float penalty_sameVert,
		float penalty_ptm, bool forceSymmetry, vector<Results_ASP> &alignsNew,
		vector<vector<TwoValues<int> > > &matches, vector<bool> &pairFlipped,
		vector<vector<float> > *dbg_matchScores = 0, ofstream *debug = 0);

/**
 * TODO: add description
 *
 *@param labels
 *@param newLabels
 */
void SplitLabels(vector<SpectrumPeakLabels> &labels,
		vector<SpectrumPeakLabels> &newLabels);

/**
 * TODO: add description
 *
 *@param specSet
 *@param consensus
 *@param peakTol
 *@param resolution
 */
void ComputeSpectralStars(SpecSet &specSet, Spectrum &consensus, float peakTol,
		float resolution);

/*
 * projectSpectrum - "projects" all specSet spectra with indices in specsToProcess onto specBase:
 *   projection of A onto B: Let M=mass(B)-mass(A). The projection of A onto B increases the scores of
 *   peaks in B by the scores of peaks in A whose masses match a modified version of A (with mod mass M).
 *
 *@param specSet    Set of all spectra
 *@param specBase   Spectrum to project onto
 *@param specsToProcess Indices of spectra (in specSet) to project from
 *@param curDeltas  Locations of the deltas (modifications) for the current subset of processed spectra (from specsToProcess, used in recursion)
 *@param bestScore  Best projection score in specBase over all possible deltas on all neighbors in specsToProcess
 *@param bestDeltas The locations of the deltas resulting in the highest score (bestScore)
 *@param peakTol    Tolerance for mass errors (in Daltons)
 *@param finalProj  Projected version of specBase after projecting all specsToProcess with bestDeltas
 *@param idxMatched Indices of matched peaks in specBase/specsToProcess.front() (only when specsToProcess.size()==1)
 */
void ProjectSpectrum(SpecSet &specSet, const Spectrum &specBase,
		list<int> &specsToProcess, list<int> &curDeltas, float &bestScore,
		list<int> &bestDeltas, float peakTol, Spectrum *finalProj = NULL,
		vector<TwoValues<int> > *idxMatched = NULL, unsigned int minNumMatchedPeaks=0);

#endif
