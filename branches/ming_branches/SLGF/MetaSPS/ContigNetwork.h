/*
 * ContigNetwork.h
 *
 *  Created on: Jan 5, 2012
 *      Author: aguthals
 */

#ifndef CONTIGNETWORK_H_
#define CONTIGNETWORK_H_

using namespace std;

#include "Contig.h"
#include "prm_alignment.h"
#include "PRMSpecEdge.h"

namespace specnets {

/**
 * This class was designed to merge contig PRM spectra into meta-contigs,
 * but it can be used to merge any group of aligned PRM spectra into contigs
 */
class ContigNetwork {
public:

	// Should not be needed but ExecAssembly needs them
	float pkTol;
	float pmTol;

	/**
	 * Default constructor. Must call initialize() after calling this
	 */
	ContigNetwork(void);

	/**
	 * Primary constructor, calls initialize with the same parameters.
	 * @param _rootSpectra PRM spectra that will be merged into contigs. Each of these is
	 *   initialized to a Contig containing one spectrum.
	 * @param alignments alignments of _rootSpectra, each of which will impose an alignment
	 *   edge between contigs.
	 * @param _assembledStars If merging congis into meta-contigs (where _rootSpectra are
	 *   the contig PRM spectra), specify this so the abinfo of output meta-contigs can be
	 *   computed.
	 */
	ContigNetwork(const SpecSet& _rootSpectra,
			const SpectrumPairSet& alignments, const abinfo_t& _assembledStars);

	/**
	 * Deconstructor. De-allocates all class variables stored on the heap
	 */
	virtual ~ContigNetwork(void);

	/**
	 * Initializes the graph to contain a node for evey contig and an edge for every
	 *   alignment. Also allocates memory for class data structures.
	 * @param _rootSpectra PRM spectra that will be merged into contigs. Each of these is
	 *   initialized to a Contig containing one spectrum.
	 * @param alignments alignments of _rootSpectra, each of which will impose an
	 *   alignment edge between contigs.
	 * @param _assembledStars If merging congis into meta-contigs (where _rootSpectra are
	 *   the contig PRM spectra), specify this so the abinfo of output meta-contigs can be
	 *   computed.
	 */
	virtual void initialize(const SpecSet& _rootSpectra,
			const SpectrumPairSet& alignments, const abinfo_t& _assembledStars);

	/**
	 * Iteratively (1) merges the highest scoring edge and (2) updates the scores of
	 *   neighboring edges until the highest scoring edges has score < minScore
	 * @param minScore minimum allowable edge/alignment score
	 * @return the number of edges that were merged
	 */
	virtual int assembleIteratively(float minScore);

	/**
	 * Merges the edge at edgeIdx. This will reverse a node connected to the edge if required
	 *   by the alignment (by calling reverseNode()). Then it will merge the two Contigs (by
	 *   calling Contig::merge()), update neighboring edges (without rescoring them), and
	 *   remove one of the two merged Contigs (the other one becoming the merged Contig)
	 * @param edgeIdx Index of edge to be merged. This edge will no longer exist if the
	 *   function is successful.
	 * @return true if an edge at edgeIdx exists and it was merged, false if not
	 */
	virtual bool mergeEdge(int edgeIdx);

	/**
	 * Re-assigns the alignment scores of all edges connected to the Contig at nodeIdx.
	 *   Alignment scores are computed by computeEdgeScore()
	 * @param nodeIdx index of Contig to re-score edges to.
	 * @return true if the Contig at nodeIdx exists and edges were re-scored, false if not
	 */
	virtual bool rescoreConnectingEdges(int nodeIdx);

	/**
	 * Generates a list of all neighbors to a given Contig node along with the indices of all
	 *   connecting edges. The returned data structure must be de-allocated by the caller.
	 * @param metaContigIdx index of Contig node
	 * @return pointer to a list of all neighbors to metaContigIdx. Each element contains
	 *   first=index of neighboring node and second=index of connecting edge. This must be
	 *   de-allocated by the caller via "delete"
	 */
	virtual list<pair<int, int> >* getNeighborList(int metaContigIdx) const;

	/**
	 * Returns true if a node with index metaContigIdx exists in the graph, false if not
	 */
	virtual bool containsNode(int metaContigIdx) const;

	/**
	 * Returns true if an edge between node1 and node2 exists in the graph, false if not
	 */
	virtual bool containsEdge(int node1, int node2) const;

	/**
	 * Returns true if an edge with index edgeIdx exists in the graph, false if not
	 */
	virtual bool containsEdge(int edgeIdx) const;

	/**
	 * Returns a pointer to the Contig at index in the graph.
	 * @param index
	 * @return Pointer to Contig at that index. If Contig does not exist, An error message
	 *   is printed and 0 is returned.
	 */
	virtual Contig* getNode(int index) const;

	/**
	 * Returns a pointer to the edge at index in the graph.
	 * @param index
	 * @return Pointer to edge at that index. If edge does not exist, An error message
	 *   is printed and 0 is returned.
	 */
	virtual PRMSpecEdge* getEdge(int index) const;

	/**
	 * Returns a pointer to the edge between node1 and node2.
	 * @param node1
	 * @param node2
	 * @return Pointer to the edge. If the edge does not exist, An error message
	 *   is printed and 0 is returned.
	 */
	virtual PRMSpecEdge* getEdge(int node1, int node2) const;

	/**
	 * Adds a Contig to the graph as a new node. This fails if a node at newNode.index
	 *   already exists.
	 * @param newNode Contig that will be copied into the graph
	 * @return true if node was added, false if not
	 */
	virtual bool addNode(const Contig& newNode);

	/**
	 * Adds an edge to the graph. This fails if an edge at newEdge.index
	 *   already exists ir an edge between nodes newEdge.spec1 and newEdge.spec2 already exists
	 * @param newEdge Edge that will be copied into the graph
	 * @return true if edge was added, false if not
	 */
	virtual bool addEdge(const PRMSpecEdge& newEdge);

	/**
	 * Removes the node at nodeIdx from the graph
	 * @param nodeIdx
	 * @return true if node existed and was removed, false if not
	 */
	virtual bool removeNode(int nodeIdx);

	/**
	 * Removes the edge at edgeIdx from the graph
	 * @param edgeIdx
	 * @return true if the edge existed and was removed, false if not
	 */
	virtual bool removeEdge(int edgeIdx);

	/**
	 * Removes the edge between node1 and node2 from the graph
	 * @param node1
	 * @param node2
	 * @return true if the edge existed and was removed, false if not
	 */
	virtual bool removeEdge(int node1, int node2);

	/**
	 * Reverses the Contig at index metaContigIdx. Also updates all connecting edges
	 * @param metaContigIdx
	 * @return true if node existed and was reversed, false if not
	 */
	virtual bool reverseNode(int metaContigIdx);

	/**
	 * Computes the shift between the consensus PRM spectra of nodes connected by a given edge
	 * @param edge edge between 2 Contig nodes
	 * @param nodeFrom either edge->spec1 or edge->spec2
	 * @param nodeTo either edge->spec2 or edge->spec1
	 * @return shift of Contig at nodeTo wrt Contig at nodeFrom
	 */
	virtual float
	getConsensusShift(PRMSpecEdge& edge, int nodeFrom, int nodeTo) const;

	/**
	 * Computes the score of the alignment between the consensus PRM spectra of Contig nodes
	 *   connected by an input edge.
	 * @param edge pointer to edge between 2 nodes
	 * @return pair.first=score of overlap with edge.spec1, pair-second=score of overlap with edge.spec2
	 */
	virtual pair<float, float> getEdgeScore(PRMSpecEdge& edge);

	/**
	 * Returns a pointer to the highest scoring edge in the graph
	 */
	virtual PRMSpecEdge* getHighestScoringEdge();

	/**
	 * Gets a new unique node index
	 */
	virtual int getNewNodeIndex();

	/**
	 * Gets a new unique edge index
	 */
	virtual int getNewEdgeIndex();

protected:

	virtual void addEdge(int spec1, int spec2, int edgeIdx);

	// Keep track of maximum node and edge indices so none overlap
	int maxNodeIdx;
	int maxEdgeIdx;

	/**
	 * Number of meta-contigs
	 */
	int size;

	/**
	 * Alignment object for re-scoring edges
	 */
	PRMAlignment alignmentObj;

	/**
	 * SPS contigs
	 */
	SpecSet* rootSpectra;

	/**
	 * abinfo for SPS contigs
	 */
	abinfo_t* assembledStars;

	/**
	 * Holds edge references for easy lookup
	 */
	map<int, PRMSpecEdge>* edges;

	/**
	 * Assembled meta-contigs
	 */
	map<int, Contig>* nodes;

	/**
	 * key = c1 meta-contig index
	 *   value->key = c2 meta-contig index
	 *     value->value = index of edge connecting c1 and c2
	 */
	map<int, map<int, int> >* graph;
};

}

#endif /* CONTIGNETWORK_H_ */
