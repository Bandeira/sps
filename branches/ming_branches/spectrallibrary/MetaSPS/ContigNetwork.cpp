/*
 * ContigNetwork.cpp
 *
 *  Created on: Jan 5, 2012
 *      Author: aguthals
 */

#include "ContigNetwork.h"

using namespace std;

namespace specnets {

ContigNetwork::ContigNetwork() :
	size(0), rootSpectra(0x0), edges(0x0), nodes(0x0), graph(0x0),
			assembledStars(0x0) {
}

ContigNetwork::ContigNetwork(const SpecSet& _rootSpectra,
		const SpectrumPairSet& alignments, const abinfo_t& _assembledStars) :
	size(0), rootSpectra(0x0), edges(0x0), nodes(0x0), graph(0x0),
			assembledStars(0x0) {
	initialize(_rootSpectra, alignments, _assembledStars);
}

ContigNetwork::~ContigNetwork() {
	if (rootSpectra != 0x0) {
		delete rootSpectra;
	}

	if (nodes != 0x0) {
		delete nodes;
	}

	if (edges != 0x0) {
		delete edges;
	}

	if (graph != 0x0) {
		delete graph;
	}

	if (assembledStars != 0x0) {
		delete assembledStars;
	}

	maxNodeIdx = 0;
	maxEdgeIdx = 0;
}

void ContigNetwork::initialize(const SpecSet& _rootSpectra,
		const SpectrumPairSet& alignments, const abinfo_t& _assembledStars) {

	size = _rootSpectra.size();
	maxNodeIdx = 0;
	maxEdgeIdx = 0;

	rootSpectra = new SpecSet(_rootSpectra.size());
	rootSpectra->operator =(_rootSpectra);

	nodes = new map<int, Contig> ;
	for (int i = 0; i < size; i++) {
		Contig newContig(i, rootSpectra);
		maxNodeIdx++;
		addNode(newContig);
	}

	graph = new map<int, map<int, int> > ;
	edges = new map<int, PRMSpecEdge> ;
	for (int i = 0; i < alignments.size(); i++) {
		PRMSpecEdge newEdge(i, alignments[i], *rootSpectra);
		maxEdgeIdx++;
		addEdge(newEdge);
	}

	assembledStars = new abinfo_t;
	assembledStars->operator =(_assembledStars);

	// Default values
	pkTol = 0.5;
	pmTol = 1.0;
}

int ContigNetwork::assembleIteratively(float minScore) {

	int numMerged = 0;

	for (map<int, Contig>::iterator nodeIt = nodes->begin(); nodeIt
			!= nodes->end(); nodeIt++) {
		nodeIt->second.setPeakTolerance(pkTol);
		nodeIt->second.setParentMassTol(pmTol);
	}

	PRMSpecEdge* edgeIt = getHighestScoringEdge();
	while (edgeIt != (PRMSpecEdge*) 0 && edgeIt->getMinScore() > minScore) {

		int edgeIdx = edgeIt->index;
		int node1 = edgeIt->spec1;

		edgeIt = (PRMSpecEdge*) 0;

		mergeEdge(edgeIdx);
		numMerged++;

		rescoreConnectingEdges(node1);

		edgeIt = getHighestScoringEdge();
	}

	return numMerged;
}

bool ContigNetwork::mergeEdge(int edgeIdx) {
	PRMSpecEdge* edge = getEdge(edgeIdx);

	if (edge == (PRMSpecEdge*) 0) {
		return false;
	}

	int node1 = edge->spec1;
	int node2 = edge->spec2;

	if (edge->spec2rev) {
		reverseNode(node2);
	}

	Contig* mergedContig = getNode(node1);
	Contig* contig2 = getNode(node2);

	mergedContig->merge(edge, contig2, pkTol, pmTol);
	mergedContig->setPeakTolerance(pkTol);
	mergedContig->setParentMassTol(pmTol);

	list<PRMSpecEdge> edgesToAdd;
	list<int> edgesToRemove;

	map<int, int>* node2Neighbors = &(*graph)[node2];
	for (map<int, int>::iterator neighIt = node2Neighbors->begin(); neighIt
			!= node2Neighbors->end(); neighIt++) {
		int node3 = neighIt->first;

		if (node3 == node1) {
			continue;
		}

		int otherEdge3Idx = neighIt->second;
		PRMSpecEdge otherEdge3(*getEdge(otherEdge3Idx));

		PRMSpecEdge::appendEdges(edge, &otherEdge3, &otherEdge3);

		if (!containsEdge(node1, node3)) {
			edgesToAdd.push_back(otherEdge3);
		} else {
			PRMSpecEdge* origEdge = getEdge(node1, node3);

			pair<float, float> scoreOrigRes = getEdgeScore(*origEdge);
			pair<float, float> scoreOther3Res = getEdgeScore(otherEdge3);
			float scoreOrig = min(scoreOrigRes.first, scoreOrigRes.second);
			float scoreOther3 =
					min(scoreOther3Res.first, scoreOther3Res.second);

			if (scoreOrig < scoreOther3) {
				edgesToAdd.push_back(otherEdge3);
				edgesToRemove.push_back(origEdge->index);
			}
		}
	}

	edge = (PRMSpecEdge*) 0;
	contig2 = (Contig*) 0;

	removeNode(node2);

	for (list<int>::iterator removeIt = edgesToRemove.begin(); removeIt
			!= edgesToRemove.end(); removeIt++) {
		removeEdge(*removeIt);
	}
	for (list<PRMSpecEdge>::iterator addIt = edgesToAdd.begin(); addIt
			!= edgesToAdd.end(); addIt++) {
		addEdge(*addIt);
	}

	return true;
}

bool ContigNetwork::rescoreConnectingEdges(int nodeIdx) {
	if (!containsNode(nodeIdx)) {
		ERROR_MSG("Node " << nodeIdx << " does not exist");
		return false;
	}

	map<int, int>* neighborRef = &(*graph)[nodeIdx];
	for (map<int, int>::iterator neighIt = neighborRef->begin(); neighIt
			!= neighborRef->end(); neighIt++) {
		PRMSpecEdge* edge = getEdge(neighIt->second);
		pair<float, float> scoreOrigRes = getEdgeScore(*edge);
		edge->score1 = scoreOrigRes.first;
		edge->score2 = scoreOrigRes.second;
	}

	return true;
}

list<pair<int, int> >* ContigNetwork::getNeighborList(int metaContigIdx) const {
	list<pair<int, int> >* neighbors = new list<pair<int, int> > ;

	if (!containsNode(metaContigIdx)) {
		return neighbors;
	}
	map<int, int>* neighborRef = &(*graph)[metaContigIdx];
	for (map<int, int>::iterator neighIt = neighborRef->begin(); neighIt
			!= neighborRef->end(); neighIt++) {
		neighbors->push_back(pair<int, int> (neighIt->first, neighIt->second));
	}
	return neighbors;
}

bool ContigNetwork::containsNode(int metaContigIdx) const {
	return graph->count(metaContigIdx) > 0;
}

bool ContigNetwork::containsEdge(int node1, int node2) const {
	if ((!containsNode(node1)) || (!containsNode(node2))
			|| (*graph)[node1].count(node2) == 0) {
		return false;
	}
	return true;
	/*

	 if (edges->count((*graph)[node1][node2]) == 0) {
	 ERROR_MSG("Graph is inconsistent");
	 return false;
	 }
	 return true;
	 */
}

bool ContigNetwork::containsEdge(int edgeIdx) const {
	if (edges->count(edgeIdx) == 0) {
		return false;
	}
	return true;
	/*
	 PRMSpecEdge* edge = (*edges)[edgeIdx];

	 if ((!containsNode(edge->spec1)) || (!containsNode(edge->spec2))
	 || (*graph)[edge->spec1].count(edge->spec2) == 0) {
	 ERROR_MSG("Graph is inconsistent");
	 return false;
	 }
	 return true;
	 */
}

Contig* ContigNetwork::getNode(int index) const {
	if (!containsNode(index)) {
		ERROR_MSG("Node " << index << " does not exist");
		return (Contig*) 0;
	} else {
		return &(*nodes)[index];
	}
}

PRMSpecEdge* ContigNetwork::getEdge(int index) const {
	if (!containsEdge(index)) {
		ERROR_MSG("Edge " << index << " does not exist");
		return (PRMSpecEdge*) 0;
	} else {
		return &(*edges)[index];
	}
}

PRMSpecEdge* ContigNetwork::getEdge(int node1, int node2) const {
	if (!containsEdge(node1, node2)) {
		ERROR_MSG("Edge between nodes " << node1 << " and " << node2 << " does not exist");
		return (PRMSpecEdge*) 0;
	} else {
		return &(*edges)[(*graph)[node1][node2]];
	}
}

bool ContigNetwork::addNode(const Contig& newNode) {
	if (containsNode(newNode.index)) {
		ERROR_MSG("Node " << newNode.index << " already exists");
		return false;
	}

	(*nodes)[newNode.index] = newNode;

	map<int, int> subGraph;
	(*graph)[newNode.index] = subGraph;

	return true;
}

bool ContigNetwork::addEdge(const PRMSpecEdge& newEdge) {
	if (containsEdge(newEdge.spec1, newEdge.spec2)) {
		ERROR_MSG("Edge between nodes " << newEdge.spec1 << " and " << newEdge.spec2 << " already exists");
		return false;
	} else if (containsEdge(newEdge.index)) {
		ERROR_MSG("Edge at index " << newEdge.index << " already exists");
		return false;
	}

	addEdge(newEdge.spec1, newEdge.spec2, newEdge.index);

	(*edges)[newEdge.index] = newEdge;

	return true;
}

bool ContigNetwork::removeNode(int nodeIdx) {

	if (!containsNode(nodeIdx)) {
		ERROR_MSG("Node " << nodeIdx << " does not exist");
		return false;
	}

	map<int, int>* neighborRef = &(*graph)[nodeIdx];
	for (map<int, int>::iterator neighIt = neighborRef->begin(); neighIt
			!= neighborRef->end(); neighIt++) {
		removeEdge(nodeIdx, neighIt->second);
	}

	graph->erase(nodeIdx);

	return true;
}

bool ContigNetwork::removeEdge(int edgeIdx) {
	PRMSpecEdge* edgeToRemove = getEdge(edgeIdx);

	if (edgeToRemove == 0) {
		return false;
	}

	/*
	 if (!containsEdge(edgeToRemove->spec1, edgeToRemove->spec2)) {
	 return false;
	 }
	 */

	int node1 = edgeToRemove->spec1;
	int node2 = edgeToRemove->spec2;

	edgeToRemove = (PRMSpecEdge*) 0;

	(*graph)[node1].erase(node2);
	(*graph)[node2].erase(node1);
	edges->erase(edgeIdx);

	return true;
}

bool ContigNetwork::removeEdge(int node1, int node2) {
	PRMSpecEdge* edgeToRemove = getEdge(node1, node2);

	if (edgeToRemove == 0) {
		return false;
	} else {
		int edgeIdx = edgeToRemove->index;
		edgeToRemove = (PRMSpecEdge*) 0;
		return removeEdge(edgeIdx);
	}
}

bool ContigNetwork::reverseNode(int metaContigIdx) {
	if (!containsNode(metaContigIdx)) {
		return false;
	}
	(*nodes)[metaContigIdx].reverse();

	map<int, int>* neighborRef = &(*graph)[metaContigIdx];
	for (map<int, int>::iterator neighIt = neighborRef->begin(); neighIt
			!= neighborRef->end(); neighIt++) {
		int edgeIdx = neighIt->second;
		PRMSpecEdge* connectingEdge = &(*edges)[edgeIdx];
		connectingEdge->reverse(metaContigIdx);
	}

	return true;
}

float ContigNetwork::getConsensusShift(PRMSpecEdge& edge, int nodeFrom,
		int nodeTo) const {

	Contig* node1 = getNode(nodeFrom);
	Contig* node2 = getNode(nodeTo);
	PRMSpecEdge* edgeUse = &edge;

	Contig node2Copy;
	PRMSpecEdge edgeCopy;

	if (edge.spec2rev) {
		node2Copy = *node2;
		node2Copy.reverse();
		node2 = &node2Copy;

		edgeCopy = edge;
		edgeCopy.reverse(nodeTo);
		edgeUse = &edgeCopy;
	}

	float totalScoreF, totalScoreR;
	float leftEdgeF = 0, rightEdgeF = 0;
	for (map<int, pair<float, float> >::iterator refIt =
			node1->rootRef->begin(); refIt != node1->rootRef->end(); refIt++) {
		if (refIt->first == edgeUse->spec1) {
			continue;
		}
		leftEdgeF = min(leftEdgeF, refIt->second.first);
	}
	for (map<int, pair<float, float> >::iterator refIt =
			node2->rootRef->begin(); refIt != node2->rootRef->end(); refIt++) {
		if (refIt->first == edgeUse->spec2) {
			continue;
		}
		rightEdgeF = min(rightEdgeF, refIt->second.first);
	}
	/*
	 if (debug) {
	 cout << "Root shift between " << idx1 << " and " << idx2 << " = "
	 << fShift << ", " << rShift << "\n";
	 cout << "consensusFShift = " << fShift << " - " << leftEdgeF << " + "
	 << rightEdgeF << " = " << fShift - leftEdgeF + rightEdgeF
	 << "\n";
	 cout << "consensusRShift = " << rShift << " - " << leftEdgeF << " + "
	 << rightEdgeR << " = " << rShift - leftEdgeF + rightEdgeR
	 << "\n";
	 }*/
	return edgeUse->getShift(nodeFrom, nodeTo) - leftEdgeF + rightEdgeF
			- node1->endGaps.first + node2->endGaps.first;
}

pair<float, float> ContigNetwork::getEdgeScore(PRMSpecEdge& edge) {

	float shift = getConsensusShift(edge, edge.spec1, edge.spec2);

	Contig* node1 = getNode(edge.spec1);
	Contig* node2 = getNode(edge.spec2);

	Contig node2Copy;

	if (edge.spec2rev) {
		node2Copy = *node2;
		node2Copy.reverse();
		node2 = &node2Copy;
	}

	alignmentObj.setSpec1((Spectrum*) node1);
	alignmentObj.setSpec2((Spectrum*) node2);

	pair<int, pair<float, float> > alignRes = alignmentObj.getShiftScore(shift,
			pmTol, 0);

	return pair<float, float> (alignRes.second.first, alignRes.second.second);
}

PRMSpecEdge* ContigNetwork::getHighestScoringEdge() {

	if (edges->size() == 0) {
		return (PRMSpecEdge*) 0;
	}

	map<int, PRMSpecEdge>::iterator edgeIt = edges->begin();

	PRMSpecEdge* bestEdge = &edgeIt->second;

	edgeIt++;

	for (; edgeIt != edges->end(); edgeIt++) {
		bestEdge
				= (edgeIt->second.getMinScore() > bestEdge->getMinScore()) ? &edgeIt->second
						: bestEdge;
	}
	return bestEdge;
}

void ContigNetwork::addEdge(int spec1, int spec2, int edgeIdx) {

	map<int, int> emptySubGraph;
	if (graph->count(spec1) == 0) {
		(*graph)[spec1] = emptySubGraph;
		(*graph)[spec1][spec2] = edgeIdx;
	} else {
		if ((*graph)[spec1].count(spec2) == 0) {
			WARN_MSG("Over-writting edge between " << spec1 << " and " << spec2 << "(edge index from " << (*graph)[spec1][spec2] << " to " << edgeIdx << ")");
		}
		(*graph)[spec1][spec2] = edgeIdx;
	}

	if (graph->count(spec2) == 0) {
		(*graph)[spec2] = emptySubGraph;
		(*graph)[spec2][spec1] = edgeIdx;
	} else {
		if ((*graph)[spec2].count(spec1) == 0) {
			WARN_MSG("Over-writting edge between " << spec2 << " and " << spec1 << "(edge index from " << (*graph)[spec2][spec1] << " to " << edgeIdx << ")");
		}
		(*graph)[spec2][spec1] = edgeIdx;
	}
}

int ContigNetwork::getNewNodeIndex() {
	int idx = maxNodeIdx;
	maxNodeIdx++;
	return idx;
}

int ContigNetwork::getNewEdgeIndex() {
	int idx = maxEdgeIdx;
	maxEdgeIdx++;
	return idx;
}

}
