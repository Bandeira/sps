/*
 * AbruijnGraph.h
 *
 *  Created on: May 1, 2012
 *      Author: aguthals
 */

#ifndef ABRUIJNGRAPH_H_
#define ABRUIJNGRAPH_H_

#include "BaseGraph.h"
#include "aminoacid.h"
#include "AbruijnNode.h"
#include "AbruijnEdge.h"
#include "SpectrumAlignment.h"
#include "SpectrumAlignmentSet.h"
#include "SpecSet.h"
#include "spectrum.h"
#include "JumpEdge.h"

#include <vector>
#include <string>
#include <stdlib.h>

using namespace std;
using namespace specnets;

namespace specnets
{
  class Node;
  class Edge;
  class BaseGraph;
}

namespace abruijn
{
  class AbruijnNode;
  class AbruijnEdge;
  class SpectrumAlignment;
  class SpectrumAlignmentSet;

  typedef list<pair<AbruijnNode*, AbruijnNode*> > AbruijnAlignment;

  class AbruijnGraph : public specnets::BaseGraph,
                       public specnets::Node
  {
  private:

    class ParallelPath
    {

    public:

      ParallelPath();

      ParallelPath(const AbruijnGraph& abG, AbruijnNode& start);

      void initialize(const AbruijnGraph& abG, AbruijnNode& start);

      void advance(const AbruijnGraph& abG, const JumpEdge& advanceEdge);

      bool addNode(const unsigned int& pathIdx, const long& node);

      void splice(const ParallelPath& other);

      vector<JumpEdge> m_path;
      vector<vector<long> > m_nodes;
      vector<pair<JumpEdge, unsigned int> > m_outgoingEdges;

      double m_weight;
    };

  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    static const unsigned int MAX_NUM_JUMPS;

    static const unsigned int MAX_NUM_MODS_PER_JUMP;

    static const double MAX_JUMP_MASS;

    static const string STARTPT_SPEC_ID;

    static const string ENDPT_SPEC_ID;

    static unsigned int PATH_EXPAND_LIMIT;

    static bool DEBUG_EXPANSION;

    static bool REVERSE_STARS;

    static bool ENFORCE_B_ENDPTS;

    static double GAP_PENALTY;

    static AbruijnGraph* castGraphPtr(specnets::BaseGraph* graphPtr);

    static AbruijnGraph* castNodePtr(specnets::Node* nodePtr);

    static string nodeSetToString(const set<AbruijnNode*>& inputNodeSet);

    AbruijnGraph();

    AbruijnGraph(const AbruijnGraph& other);

    AbruijnGraph(const pair<pair<vector<int> , vector<int> > , vector<pair<
                     vector<int> , vector<double> > > >& inAbinfo,
                 const SpecSet& assembledSpecs,
                 const float& peakTol,
                 pair<pair<vector<int> , vector<int> > , vector<pair<
                     vector<int> , vector<double> > > >* outAbinfo = 0);

    virtual ~AbruijnGraph(void)
    {
    }

    void initialize();

    virtual void compress(vector<long>* outputNewNodeIdxs = 0,
                          vector<long>* outputNewEdgeIdxs = 0);

    virtual void clear(void)
    {
      initialize();
    }

    void
    reSequencePPaths(const pair<pair<vector<int> , vector<int> > , vector<pair<
                         vector<int> , vector<double> > > >& abinfo,
                     const SpecSet& assembledSpecs,
                     const float& peakTol,
                     pair<pair<vector<int> , vector<int> > , vector<pair<
                         vector<int> , vector<double> > > >* outAbinfo = 0);

    virtual specnets::Node* cloneNode(specnets::Node& copyNode) const;

    virtual specnets::Edge* cloneEdge(specnets::Edge& copyEdge) const;

    virtual specnets::Node* createNode(void) const;

    virtual specnets::Edge* createEdge(void) const;

    inline AbruijnNode* getNode(unsigned long nodeIndex) const
    {
      return AbruijnNode::castNodePtr(BaseGraph::getNode(nodeIndex));
    }

    inline AbruijnEdge* getEdge(unsigned long edgeIndex) const
    {
      return AbruijnEdge::castEdgePtr(BaseGraph::getEdge(edgeIndex));
    }

    AbruijnGraph &operator=(const AbruijnGraph &other);

    virtual void copy(specnets::BaseGraph& otherGraph);

    virtual void
    addBinaryVersionInfo(map<string, unsigned short>& versions) const;

    void outputConsensusSpectrum(Spectrum& outputSpec) const;

    AbruijnNode* getMultiJumpDestination(AbruijnEdge* abEdge,
                                         bool goReverse,
                                         double& foundJumpMass,
                                         double& foundWeight,
                                         string& foundPath,
                                         string& foundPathFor,
                                         string& foundPathRev) const;

    //void appendGraph(const AbruijnGraph& otherGraph);

    //void removeEmptyNodes(void);
    /*
     void loadAssembly(const specnets::SpecSet& prmSpectra,
     const SpectrumAlignmentSet& prmAligns,
     bool addEndpointGlues = true);

     //string getConsensusSeq(const bool reportMods = false) const;
     * */

  protected:

    //SpecSet m_assembledSpecs;
    //SpectrumAlignmentSet m_alignments;

    Path m_consensusPath;
    vector<AbruijnNode*> m_consensusPathVerts;
    tr1::unordered_map<string, map<MZRange,
        pair<AbruijnNode*, TwoValues<bool> > > > m_nodeLookup;

    set<AbruijnNode*> m_nodesToExpand;
    map<string, double> m_startPtAssembly;

    specnets::AAJumps& m_globalJumps;

    bool m_reversed;

    double m_consensusPathScore;

    void
    m_reSequencePPaths(const pair<pair<vector<int> , vector<int> > , vector<
                           pair<vector<int> , vector<double> > > >& abinfo,
                       const SpecSet& assembledSpecs,
                       const float& peakTol,
                       const bool reverseStars,
                       const bool useBEndpts,
                       pair<pair<vector<int> , vector<int> > , vector<pair<
                           vector<int> , vector<double> > > >* outAbinfo = 0);

    void m_internalCopy(const AbruijnGraph &other);

    void addSpectrum(const specnets::Spectrum& prmSpec, bool allEndPts = false);

    inline void addSpectra(const specnets::SpecSet& prmSpectra)
    {
      for (unsigned int i = 0; i < prmSpectra.size(); i++)
      {
        addSpectrum(prmSpectra[i]);
      }
    }

    void addGlues(const SpectrumAlignment& prmAlign);

    inline void addGlues(const SpectrumAlignmentSet& prmAligns)
    {
      for (unsigned int i = 0; i < prmAligns.size(); i++)
      {
        addGlues(prmAligns[i]);
      }
    }

    AbruijnNode* addAbruijnNode(AbruijnNode& copyNode);

    AbruijnEdge* addLabelFreeAbEdge(AbruijnNode* from,
                                    AbruijnNode* to,
                                    const double& mass,
                                    const double& weight,
                                    const double& rWeight);

    void recoverSourceNodeScores();

    void reverseEdgeWeights();

    void addEndpointEdges(const set<string>& assembledSpecIDs,
                          const SpecSet& assmebledSpecs);

    void addAbruijnEdges(const Spectrum& prmSpec,
                         const unsigned int& peakIdxFrom,
                         const unsigned int& peakIdxTo,
                         const double& avgSpecScore,
                         list<AbruijnEdge*>* addedEdges = 0);

    void expandForwardReversePaths(const float& peakTol);

    void expandPaths(AbruijnNode* source, const bool& goReverse);

    void mergePaths(AbruijnNode* start, const bool& goReverse);

    void processAllPaths(const bool& goReverse, const bool& expandYes);

    /**
     * Combines edges with the same label between the same two nodes (adds together weight)
     */
    void mergeParallelEdges(AbruijnNode* from, AbruijnNode* to, bool mergeMods =
        false);

    /*
     * Calls mergeParallelEdges on all connected node pairs
     */
    void mergeAllParallelEdges(bool mergeMods = false);

    /*
     * Prunes any edge between two nodes that is not the maximal scoring edge over all others
     *   with the same jump mass
     */
    void pruneParallelEdgesByMass(AbruijnNode* from,
                                  AbruijnNode* to,
                                  float peakTol);

    /*
     * Calls pruneParallelEdgesByMass on all connected node pairs
     */
    void pruneAllParallelEdgesByMass(float peakTol);

    /**
     * Combines edges with no label but the same mass between the same two nodes (adds together weight)
     */
    void mergeParallelLabelFreeEdgesByMass(AbruijnNode* from,
                                           AbruijnNode* to,
                                           float peakTol);

    /**
     * Calls mergeParallelLabelFreeEdgesByMass on all connected node pairs
     */
    void mergeAllParallelLabelFreeEdgesByMass(float peakTol);

    bool computeHeaviestPathDAG();

    //void assembleConsensus(const bool& addEndpointGlues);

    void removeExpandedPath(AbruijnEdge** pathPartEdge);

    void injectExpandedPath(AbruijnEdge* templateEdge,
                            const unsigned int& step,
                            const bool& goReverse);
    /*
     void getEndPointSpectra(map<string, specnets::MZRange>& leftAlignedShifts,
     specnets::Spectrum& outputStartPtSpec,
     specnets::Spectrum& outputEndPtSpec,
     SpectrumAlignmentSet& outputStartPtAligns,
     SpectrumAlignmentSet& outputEndPtAligns);
     */

    void removeSymmetricYNodes();

    void removeSourceGreenNodes();

    void glueNodes(AbruijnNode* node1, AbruijnNode** node2);

    AbruijnNode* mergeNodes(const list<AbruijnEdge*>& edgeList,
                            AbruijnNode* source,
                            const string& fLabel,
                            const string& rLabel,
                            const string& pathLabel,
                            set<AbruijnEdge*>& edgesToRemove,
                            list<AbruijnEdge*>& edgesToRemoveNext,
                            const bool& lookReverse);

  };
}

#endif /* ABRUIJNGRAPH_H_ */
