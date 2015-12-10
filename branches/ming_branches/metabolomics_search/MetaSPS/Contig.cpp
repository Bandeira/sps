/*
 * Contig.cpp
 *
 *  Created on: Jun 23, 2011
 *      Author: aguthals
 */

#include "Contig.h"
#include "Logger.h"
#include "ExecFramework/ExecAssembly.h"
#include "prm_alignment.h"

namespace specnets {
Contig::Contig(void) :
	index(0), reversed(false), assembledStars(0x0), innerEdges(0x0), rootRef(
			0x0), childSpectra(0x0), Spectrum() {
	create();
}

Contig::Contig(int idx, SpecSet* contigs, abinfo_t* _assembledStars) :
	index(idx), reversed(false), assembledStars(0x0), innerEdges(0x0), rootRef(
			0x0), childSpectra(0x0), Spectrum() {
	create();
	initialize(idx, contigs, _assembledStars);
}

Contig::Contig(Contig* other) :
	index(other->index), reversed(false), assembledStars(0x0), innerEdges(0x0),
			rootRef(0x0), childSpectra(0x0), Spectrum() {
	create();
	copy(other);
}

Contig::~Contig(void) {
	delete assembledStars;
	delete innerEdges;
	delete rootRef;
	delete childSpectra;
}

Contig& Contig::operator=(const Contig &other) {
	if (this == &other) {
		return (*this);
	}
	copy(other);
	return (*this);
}

void Contig::initialize(int idx, SpecSet* contigs, abinfo_t* _assembledStars) {
	if (idx < 0 || idx >= contigs->size()) {
		ERROR_MSG("Invalid index " << idx
				<< " for meta-contig initialization from specset of size "
				<< contigs->size());
		return;
	}
	if (_assembledStars != 0 && _assembledStars->count(idx) == 0) {
		ERROR_MSG("Invalid index " << idx
				<< " not found in abinfo for meta-contig initialization");
		return;
	}
	((Spectrum &) (*this)) = (*contigs)[idx];
	index = idx;
	reversed = false;
	assembledStars->clear();
	if (_assembledStars != 0) {
		(*assembledStars)[0] = (*_assembledStars)[index];
	}
	innerEdges->clear();
	rootRef->clear();
	(*rootRef)[index] = pair<float, float> (0, 0);
	childSpectra->clear();
	(*childSpectra)[index] = pair<Spectrum, bool> ((*contigs)[idx], false);
	endGaps = pair<float, float> (0, 0);
}

void Contig::copy(const Contig& other) {
	((Spectrum &) (*this)) = (Spectrum) other;
	index = other.index;
	reversed = other.reversed;
	endGaps = other.endGaps;
	assembledStars->clear();
	assembledStars->insert(other.assembledStars->begin(),
			other.assembledStars->end());

	innerEdges->clear();
	innerEdges->insert(innerEdges->begin(), other.innerEdges->begin(),
			other.innerEdges->end());

	rootRef->clear();
	rootRef->insert(other.rootRef->begin(), other.rootRef->end());

	childSpectra->clear();
	childSpectra->insert(other.childSpectra->begin(), other.childSpectra->end());
}

void Contig::reverse(void) {
	Spectrum::reverse(0.0 - AAJumps::massH2O);
	reversed = !reversed;

	float temp = endGaps.first;
	endGaps.first = endGaps.second;
	endGaps.second = temp;

	bool rev1;
	for (map<int, pair<Spectrum, bool> >::iterator childIt =
			childSpectra->begin(); childIt != childSpectra->end(); childIt++) {
		childIt->second.first.reverse(0.0 - AAJumps::massH2O);
		rev1 = childIt->second.second;
		childIt->second.second = !rev1;
	}

	for (list<PRMSpecEdge>::iterator innerIt = innerEdges->begin(); innerIt
			!= innerEdges->end(); innerIt++) {
		innerIt->reverse(true, true);
	}

	for (map<int, pair<float, float> >::iterator refIt = rootRef->begin(); refIt
			!= rootRef->end(); refIt++) {
		temp = refIt->second.second;
		refIt->second.second = refIt->second.first;
		refIt->second.first = temp;
	}
}

void Contig::assembleConsensus(float pkTol, float pmTol) {
	SpecSet inputSpecs(childSpectra->size());
	int idxUse = 0;
	vector<int> idxRef(inputSpecs.size());
	for (map<int, pair<Spectrum, bool> >::iterator childIt =
			childSpectra->begin(); childIt != childSpectra->end(); childIt++) {
		inputSpecs[idxUse] = childIt->second.first;
		idxRef[idxUse] = childIt->first;
	}

	inputSpecs.setPeakTolerance(pkTol, false);

	SpectrumPairSet inputPairs(innerEdges->size());
	idxUse = 0;
	for (list<PRMSpecEdge>::iterator innerIt = innerEdges->begin(); innerIt
			!= innerEdges->end(); innerIt++) {
		inputPairs[idxUse] = (SpectrumPair &) (*innerIt);
		++idxUse;
	}

	ParameterList inputParams;
	inputParams.addIfDoesntExist("PENALTY_PTM", "-200");
	inputParams.addIfDoesntExist("PENALTY_SAME_VERTEX", "-1000000");
	inputParams.addIfDoesntExist("GRAPH_TYPE", "2");
	inputParams.addIfDoesntExist("MAX_AA_JUMP", "1");
	inputParams.addIfDoesntExist("MAX_MOD_MASS", "100.0");
	inputParams.addIfDoesntExist("TOLERANCE_PEAK", parseFloat(pkTol, 5));
	inputParams.addIfDoesntExist("TOLERANCE_PM", parseFloat(pmTol, 5));
	inputParams.addIfDoesntExist("MIN_MATCHED_PEAKS", "0");
	inputParams.addIfDoesntExist("MIN_EDGES_TO_COMPONENT", "0");
	inputParams.addIfDoesntExist("PATH_MIN_SPECS", "0");
	inputParams.addIfDoesntExist("PATH_MIN_PEAKS", "0");
	inputParams.addIfDoesntExist("SPEC_TYPE_MSMS", "0");
	inputParams.addIfDoesntExist("NO_SEQUENCING", "0");
	inputParams.addIfDoesntExist("ADD_ENDPOINTS", "0");
	inputParams.addIfDoesntExist("OUTPUT_COMPLETE_ABRUIJN", "");

	Clusters outputClusters;
	abinfo_t outputAbinfo;

	ExecAssembly assemblyObj(inputParams, &inputSpecs, &inputPairs,
			&outputClusters, &outputAbinfo);

	//Logger currentLogger(Logger::getDefaultLogger());
	//Logger::setDefaultLogger(Logger::getLogger(1));

	assemblyObj.invoke();

	//Logger::setDefaultLogger(currentLogger);

	if (outputClusters.consensus.size() == 0
			|| outputClusters.consensus[0].size() == 0) {
		WARN_MSG("\ExecAssembly RETURNED 0 COMPONENTS, MERGING FOR COMPONENT " << index
				<< " FAILED!\n");
		return;
	} else if (outputClusters.consensus.size() > 1) {
		WARN_MSG("\nMASAB SPLIT COMPONENT " << index << " INTO " << outputClusters.consensus.size()
				<< " CONTIGS, going with the first one\n");
	}

	(Spectrum&) (*this) = outputClusters.consensus[0];
	assembledStars->clear();
	(*assembledStars)[0] = outputAbinfo[outputAbinfo.begin()->first];

	parentMass = peakList.back()[0] + AAJumps::massMH;

	setPeakTolerance(pkTol);
	setParentMassTol(pmTol);

	float minLeftEdgeF = 0;
	float maxRightEdge = peakList.back()[0];
	for (map<int, pair<Spectrum, bool> >::iterator childIt =
			childSpectra->begin(); childIt != childSpectra->end(); childIt++) {

		float FFShift = (childIt->first == index) ? 0
				: (*rootRef)[childIt->first].first;
		minLeftEdgeF = min(minLeftEdgeF, FFShift);
		float rightEdge = FFShift
				+ childIt->second.first.back()->operator [](0);
		maxRightEdge = max(maxRightEdge, rightEdge);
		//if (debug) {
		//	cout << "\nContig " << *nodeIt << "(shift=" << FFShift << "):\n";
		//	oriented[*nodeIt].output(cout);
		//}
	}

	maxRightEdge -= minLeftEdgeF;

	pair<vector<int> , vector<double> > first_mass =
			(*assembledStars)[0].second.front();
	int locIdxF = first_mass.first.front();
	endGaps.first = first_mass.second.front() + (*rootRef)[locIdxF].first
			- minLeftEdgeF;
	endGaps.second = maxRightEdge - back()->operator [](0) - endGaps.first;

	/*
	 if (debug) {
	 cout << "\ngapF = " << gapF << " = " << first_mass.second.front()
	 << " + " << contigIndent[locIdxF][0] << "\n";
	 cout << "gapR = " << gapR << " = " << maxRightEdge << " - "
	 << putSpec.peakList.back()[0] << " - " << gapF << "\n";
	 cout << "First mass from " << locIdxF << "\n";
	 cout << "minLeftEdgeF " << minLeftEdgeF << "\n";
	 }
	 */
}

void Contig::merge(PRMSpecEdge* inEdge, Contig* inOther, float pktol, float pmTol) {

	Contig* other = inOther;
	Contig tempOther;

	PRMSpecEdge* edge = inEdge;
	PRMSpecEdge tempEdge;

	if (edge->spec2rev) {
		tempOther = *inOther;
		tempOther.reverse();
		other = &tempOther;

		tempEdge = *inEdge;
		tempEdge.reverse(other->index);
		edge = &tempEdge;
	}

	for (map<int, pair<float, float> >::iterator refIt =
			other->rootRef->begin(); refIt != other->rootRef->end(); refIt++) {
		(*rootRef)[refIt->first] = pair<float, float> (edge->getShift(index, other->index)
				+ refIt->second.first, edge->getReversedShift(index, other->index) + refIt->second.second);
	}

	PRMSpecEdge bestConn;
	int numMP = -1;
	pair<int, pair<float, float> > shiftScore;
	PRMAlignment alignmentObj;
	alignmentObj.spec2rev = false;

	for (map<int, pair<Spectrum, bool> >::iterator childIt =
			childSpectra->begin(); childIt != childSpectra->end(); childIt++) {

		alignmentObj.spec1 = childIt->first;
		alignmentObj.setSpec1(&childIt->second.first);

		for (map<int, pair<Spectrum, bool> >::iterator otherIt =
				other->childSpectra->begin(); otherIt
				!= other->childSpectra->end(); otherIt++) {

			alignmentObj.spec2 = otherIt->first;
			alignmentObj.setSpec2(&otherIt->second.first);
			alignmentObj.shift1 = (*rootRef)[otherIt->first].first
					- (*rootRef)[childIt->first].first;
			alignmentObj.shift2 = (*rootRef)[otherIt->first].second
					- (*rootRef)[childIt->first].second;

			shiftScore = alignmentObj.getShiftScore(alignmentObj.shift1, pktol,
					0);

			if (shiftScore.first > numMP) {
				numMP = shiftScore.first;
				bestConn = (SpectrumPair&) alignmentObj;
			}
		}
	}

	if (numMP < 0) {
		ERROR_MSG("Could not find overlapping contigs when merging meta-contigs " << index << " and " << other->index);
	}

	innerEdges->push_back(bestConn);

	childSpectra->insert(other->childSpectra->begin(),
			other->childSpectra->end());

	assembleConsensus(pktol, pmTol);

}

void Contig::create(void) {
	assembledStars = new abinfo_t;
	innerEdges = new list<PRMSpecEdge> ;
	rootRef = new map<int, pair<float, float> > ;
	childSpectra = new map<int, pair<Spectrum, bool> > ;
}
}
