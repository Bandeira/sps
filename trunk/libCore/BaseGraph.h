/*
 * BaseGraph.h
 *
 *  Created on: Feb 8, 2012
 *      Author: aguthals
 */

#ifndef BASEGRAPH_H_
#define BASEGRAPH_H_

#include <set>
#include <map>
#include <vector>
#include <stdlib.h>

#include "Node.h"
#include "Edge.h"
#include "Logger.h"
#include "tuple.h"

using namespace std;

namespace specnets
{
  class Node;
  class Edge;
  class EMExcep;
  class ENFExcep;
  class ENFExcep;

  typedef vector<Edge*> Path;
  typedef map<Node*, pair<double, Edge*> > Tree;

  typedef vector<Node*> NodeSet;
  typedef vector<Edge*> EdgeSet;
  typedef vector<long> IndexVector;
  typedef list<unsigned long> FreeList;
  typedef vector<pair<IndexVector, FreeList> > AdjacList;

  class BaseGraph
  {
  public:

    static const unsigned short BIN_VERSION;

    static const unsigned short BIN_SUBVERSION;

    static const string BIN_VERSION_ID;

    static const string BIN_SUBVERSION_ID;

    BaseGraph(unsigned long numNodes = 0, unsigned long numEdges = 0);

    BaseGraph(const BaseGraph& other);

    virtual ~BaseGraph(void);

    virtual string toString(void) const;

    virtual Node* cloneNode(Node& copyNode) const;

    virtual Edge* cloneEdge(Edge& copyEdge) const;

    virtual Node* createNode(void) const;

    virtual Edge* createEdge(void) const;

    virtual void validateGraph(void) const;

    virtual void validateEdge(Edge* edge) const;

    BaseGraph &operator=(const BaseGraph &other);

    virtual void copy(BaseGraph& other);

    virtual void appendGraph(const BaseGraph& otherGraph);

    virtual bool saveGraphviz(const string& filename) const;

    virtual bool saveBinaryFile(const string& filename);

    virtual void
    addBinaryVersionInfo(map<string, unsigned short>& versions) const;

    virtual bool loadBinaryFile(const string& filename);

    virtual bool saveGraphToBinaryStream(FILE* fp);

    virtual bool
    loadGraphFromBinaryStream(FILE* fp, map<string, unsigned short>& versions);

    virtual void modifiedGraph();

    void initialize(unsigned long numNodes = 0, unsigned long numEdges = 0);

    Node* addNode(Node& copyNode);

    Node* addNode(void);

    Edge* addEdge(const unsigned long& originIdx, const unsigned long& destIdx);

    Edge* addEdge(const unsigned long& originIdx,
                  const unsigned long& destIdx,
                  Edge& copyEdge);

    Edge* moveEdge(const unsigned long& edgeIdx,
                   const unsigned long& originIdx,
                   const unsigned long& destIdx);

    inline Edge* addEdge(Node* origin, Node* dest)
    {
      return addEdge(origin->getIndex(), dest->getIndex());
    }

    inline Edge* addEdge(Node* origin, Node* dest, Edge& copyEdge)
    {
      return addEdge(origin->getIndex(), dest->getIndex(), copyEdge);
    }

    void removeNode(unsigned long nodeIdx);

    inline void removeNode(Node** node)
    {
      removeNode((*node)->getIndex());
      *node = 0;
    }

    void removeEdge(unsigned long edgeIdx);

    inline void removeEdge(Edge** edge)
    {
      removeEdge((*edge)->getIndex());
      *edge = 0;
    }

    inline Node* getNode(unsigned long nodeIndex) const
    {
      if (nodeIndex >= m_nodes.size() || m_nodes[nodeIndex] == 0)
      {
        ERROR_MSG("Cannot find node " << nodeIndex);
        abort();
      }
      return m_nodes[nodeIndex];
    }

    inline Edge* getEdge(unsigned long edgeIndex) const
    {
      if (edgeIndex >= m_edges.size() || m_edges[edgeIndex] == 0)
      {
        ERROR_MSG("Cannot find edge " << edgeIndex);
        abort();
      }
      return m_edges[edgeIndex];
    }

    inline Node* lookupNode(unsigned long nodeIndex) const
    {
      return m_nodes[nodeIndex];
    }

    inline Edge* lookupEdge(unsigned long edgeIndex) const
    {
      return m_edges[edgeIndex];
    }

    inline bool containsNode(Node* node) const
    {
      return (node != 0) && (node->getIndex() < m_nodes.size())
          && (m_nodes[node->getIndex()] == node);
    }

    inline bool containsEdge(Edge* edge) const
    {
      return (edge != 0) && (edge->getIndex() < m_edges.size())
          && (m_edges[edge->getIndex()] == edge);
    }

    inline const IndexVector& getOutEdges(unsigned long nodeIndex) const
    {
      if (nodeIndex >= m_nodes.size() || m_nodes[nodeIndex] == 0)
      {
        ERROR_MSG("Cannot find node " << nodeIndex);
        abort();
      }
      return m_outEdges[nodeIndex].first;
    }

    inline const IndexVector& getOutEdges(const Node* node) const
    {
      return getOutEdges(node->m_index);
    }

    inline const IndexVector& getInEdges(unsigned long nodeIndex) const
    {
      if (nodeIndex >= m_nodes.size() || m_nodes[nodeIndex] == 0)
      {
        ERROR_MSG("Cannot find node " << nodeIndex);
        abort();
      }
      return m_inEdges[nodeIndex].first;
    }

    inline const IndexVector& getInEdges(const Node* node) const
    {
      return getInEdges(node->m_index);
    }

    inline unsigned long getNumOutEdges(const Node* node) const
    {
      return getOutEdges(node).size() - m_outEdges[node->m_index].second.size();
    }

    inline unsigned long getNumInEdges(const Node* node) const
    {
      return getInEdges(node).size() - m_inEdges[node->m_index].second.size();
    }

    inline unsigned long numEdges(void) const
    {
      return m_edges.size() - m_freeEdges.size();
    }

    inline unsigned long numNodes(void) const
    {
      return m_nodes.size() - m_freeNodes.size();
    }

    inline unsigned long maxEdgeIdx(void) const
    {
      return m_edges.size() - 1;
    }

    inline unsigned long maxNodeIdx(void) const
    {
      return m_nodes.size() - 1;
    }

    inline string getLabel(void) const
    {
      return m_label;
    }

    virtual void compress(vector<long>* outputNewNodeIdxs = 0,
                          vector<long>* outputNewEdgeIdxs = 0);

    void getLightestPaths(Node* source, Tree& outputPaths);

    pair<bool, double> getHeaviestPathDAG(Node* source,
                                          Node* sink,
                                          Path& outputPath);

    bool getTopologicalOrderingDAG(list<Node*>& outputOrder);

    Node* beginBFS(Node* source);

    Node* beginDFS(Node* source);

    Node* nextBFSDFS(void);

  protected:
    string m_label;

  private:

    NodeSet m_nodes;
    EdgeSet m_edges;

    FreeList m_freeNodes;
    FreeList m_freeEdges;

    AdjacList m_outEdges;
    AdjacList m_inEdges;

    list<Edge*> m_queueStack;
    set<Node*> m_visited;
    bool m_QueueIsStack;
    bool m_okTraverse;

    void m_beginTraverse(Node* source);

    void m_clearNodes(void);

    unsigned long m_getFreeNodeIdx(void);

    unsigned long m_getFreeEdgeIdx(void);

    unsigned long m_getFreeInAdacIdx(unsigned long nodeIdx);

    unsigned long m_getFreeOutAdacIdx(unsigned long nodeIdx);

  };
}

#endif /* BASEGRAPH_H_ */
