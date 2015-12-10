/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 10/12/09
 */

#include "CombineContigs.h"

using namespace std;

const char* DEBUG_OUTPUT = "debug_combinations.csv";
const float minContigOverlap = 400.0;
const float minSecondaryRatio = 0;
const float minParentMassSep = 10.0;

struct SortEdges : public std::binary_function<TwoValues<float> , TwoValues<
    float> , bool>
{
  bool operator()(TwoValues<float> left, TwoValues<float> right) const
  {
    return left[1] > right[1];
  }
  ;
};

void copy(CombineContigsParams& from, CombineContigsParams& to)
{
  to.contigs = from.contigs;
  to.star_spectra = from.star_spectra;
  //to.consensus = from.consensus;
  to.in_abinfo = from.in_abinfo;
  //to.out_abinfo = from.out_abinfo;
  //to.out_reversed = from.out_reversed;
  to.contig_alignments = from.contig_alignments;

  to.peak_tol = from.peak_tol;
  to.parent_mass_tol = from.parent_mass_tol;
  to.contig_overlap_score = from.contig_overlap_score;
  to.min_matched_peaks_contig_align = from.min_matched_peaks_contig_align;
  to.resolution = from.resolution;

  to.matepairs = from.matepairs;
  to.matepair_alignments = from.matepair_alignments;

  to.min_matepair_charge = from.min_matepair_charge;
  to.min_ratio_matepair_align = from.min_ratio_matepair_align;
  to.min_ratio_matepair_contig_align = from.min_ratio_matepair_contig_align;
  to.min_matched_peaks_matepair_align = from.min_matched_peaks_matepair_align;
  to.min_num_matepairs = from.min_num_matepairs;
  to.min_edges_to_component = from.min_edges_to_component;
  to.min_component_size = from.min_component_size;

  to.contig_prot_idx = from.contig_prot_idx;
  to.contig_peak_idx = from.contig_peak_idx;
  to.proteins = from.proteins;
  to.protein_idx_ignore = from.protein_idx_ignore;
}

/*
 struct CombineContigsParams {
 //required contig-contig merging parameters
 SpecSet* contigs;
 SpecSet consensus;
 abinfo_t* in_abinfo;
 abinfo_t out_abinfo;
 vector<int> out_reversed;
 list<Results_OCC>* contig_alignments;

 float peak_tol;
 float parent_mass_tol;
 float resolution;
 float contig_overlap_score;
 int min_matched_peaks_contig_align;

 //required contig-matepair merging parameters
 SpecSet* matepairs;
 list<Results_OCC>* matepair_alignments;

 short min_matepair_charge;
 float min_ratio_matepair_align;
 float min_ratio_matepair_contig_align;
 int min_matched_peaks_matepair_align;
 int min_num_matepairs;
 int min_edges_to_component;
 int min_component_size;

 //required matchma parameters
 vector<vector<int>* contig_prot_idx;
 SpecSet* contig_peak_idx;
 DB_fasta* proteins;
 set<int>* protein_idx_ignore;
 };
 */

CombineContigs::CombineContigs()
{
}
CombineContigs::~CombineContigs()
{
}

CombineContigs::CombineContigs(CombineContigsParams* params)
{
  construct(params);
}

void CombineContigs::construct(CombineContigsParams* params)
{
  merging_params = params;

  if (params->contig_prot_idx != NULL && params->contig_peak_idx != NULL
      && params->proteins != NULL && params->protein_idx_ignore != NULL)
  {
    overlaps = *(params->contig_peak_idx);
    prot_match = *(params->contig_prot_idx);
    fasta = *(params->proteins);
    idxNotUsed = *(params->protein_idx_ignore);
    haveRes = true;
  }
  else
  {
    haveRes = false;
  }

  contigs = *(merging_params->contigs);
  contig_abinfo = *(params->in_abinfo);

  if (params->contig_alignments != NULL)
    contig_alignments = *(params->contig_alignments);

  if (params->matepairs != NULL && params->matepair_alignments != NULL)
  {
    connectors = *(params->matepairs);
    connector_alignments = *(params->matepair_alignments);
    haveMatePairs = true;
  }
  else
  {
    haveMatePairs = false;
  }

  for (unsigned int i = 0; i < contigs.size(); i++)
  {
    if (contigs[i].size() == 0)
      continue;
    contigs[i].parentMass = contigs[i].peakList.back()[0] + AAJumps::massMH;
    //contigs[i].annotation.resize(0);
    //contigs[i].ionTypes.resize(0);
  }
  peakTol = params->peak_tol;
  parentMassTol = params->parent_mass_tol;
  intPeakTol = floatToInt(peakTol / merging_params->resolution);
  intParentMassTol = floatToInt(parentMassTol / merging_params->resolution);
  init();
}

// mergeType = 0 (combine w/ contig shifts), 1 (combine w/ connector shifts)
void CombineContigs::combineEdges(short mergeType)
{

  init();

  short precCharge = merging_params->min_matepair_charge;
  float minRatioConn = merging_params->min_ratio_matepair_align;
  float minRatioConnContig = merging_params->min_ratio_matepair_contig_align;
  int minNumMatchedPeaksConn = merging_params->min_matched_peaks_matepair_align;
  int minNumConnectors = merging_params->min_num_matepairs;
  int minEdgesToComponent = merging_params->min_edges_to_component;
  float minCombScore = merging_params->contig_overlap_score;
  int minComp = merging_params->min_component_size;
  int minContigMP = merging_params->min_matched_peaks_contig_align;

  float contigLen = 0, numContigs = 0, maxContigL = 0, contigA = 0;
  int maxI = 0, maxL = 0, maxLI = 0;
  for (int i = 0; i < contigs.size(); i++)
  {
    if (contigs[i].size() == 0)
      continue;
    contigLen += contigs[i].parentMass;
    numContigs += 1.0;
    if (contigs[i].parentMass > maxContigL)
    {
      maxContigL = contigs[i].parentMass;
      maxI = i;
    }
    if (haveRes && prot_match[i][0] >= 0)
    {
      float cover = overlaps[i][overlaps[i].size() - 1][1] - overlaps[i][0][1]
          + 1.0;
      contigA += cover;
      if ((int)cover > maxL)
      {
        maxL = (int)cover;
        maxLI = i;
      }
    }
  }
  printf("\nMaximum contig length (Da) : %.1f (index %d)\n", maxContigL, maxI);
  printf("Average contig length (Da) : %.1f (%.0f contigs found of possible %d)\n",
         contigLen / numContigs,
         numContigs,
         contigs.size());
  if (haveRes)
  {
    printf("Maximum contig AA length : %d (index %d)\n", maxL, maxLI);
    printf("Average contig AA length : %.1f\n", contigA / numContigs);
  }
  root_alignments.clear();
  root_alignments.resize(contigs.size());
  list<TwoValues<int> > comp_alignments;
  TwoValues<int> node_pair;

  graph.clear();
  map<int, map<int, int> > graphBefMerg;
  map<int, int>::iterator mapIt;
  edges.resize(0);
  vector<vector<float> > edgesBefMerg;
  list<TwoValues<float> > scoredEdges;
  list<TwoValues<float> >::iterator scoredEdgeIt;
  components.clear();
  components.resize(contigs.size());
  vector<bool> contigsMerged(contigs.size());
  for (int i = 0; i < contigsMerged.size(); i++)
    contigsMerged[i] = false;
  vector<set<int> > componentsBefMerg;
  set<int> erasedEdges;
  set<int>::iterator nodeIt;
  set<int>::iterator nodeIt2;
  TwoValues<float> res, res1, res2;
  vector<float> edgeTo(6);
  vector<float> edgeFrom(6);
  reversed.clear();
  reversed.resize(contigs.size());
  for (int i = 0; i < reversed.size(); i++)
    reversed[i] = false;
  vector<bool> reversedContigMerg;
  vector<bool> reversedBefMerg;
  SpecSet myContigsContigMerg;
  SpecSet myContigsBefMerg;
  for (int i = 0; i < contigs.size(); i++)
    components[i].insert(i);

  int maxContigEdgeIdx = -1;

  if (mergeType == 0)
  {
    getCondensedContigAlignmentGraph(graph,
                                     edges,
                                     scoredEdges,
                                     minCombScore,
                                     minContigMP);
  }
  vector<vector<float> > originalEdges = edges;

  if (mergeType == 1)
  {
    scoredEdge[0] = -1.0;
    scoredEdge[1] = -1.0;
    scoredEdges.push_back(scoredEdge);
  }

  cout << "Initialized with " << edges.size() << " edges\n";

  graphBefMerg = graph;
  edgesBefMerg = edges;
  componentsBefMerg = components;
  reversedBefMerg = reversed;
  myContigsBefMerg = oriented;
  //outputGraph(cout, graph, edges);

  vector<int> edgeRef(edges.size());
  for (int i = 0; i < edgeRef.size(); i++)
    edgeRef[i] = -1;

  float FFShift, FRShift, RFShift, RRShift, oldFFShift, oldFRShift, oldRFShift,
      oldRRShift, conFFShift, conFRShift, conRFShift, conRRShift, edgeScore;
  float goodEdges = 0, totalEdges = 0, edgesAdded = 0, goodEdgesAdded = 0,
      zero = 0;

  FILE* output = fopen(DEBUG_OUTPUT, "w");

  fprintf(output, "Adding Contig Edges\n");
  cout << "\nadding contig edges...\n";

  if (haveRes)
  {
    fprintf(output,
            "Protein Idx1%sProtein Idx2%sContig Idx1%sContig Idx2%sExperimental Shift FF%sExperimental Shift RR%sExpected Shift FF%sExpected Shift RR%sEqual?%sRemoved?%sContig2 Reversed?%sContig1 Reversed Global?%sContig2 Reversed Global?%sForward Score%sReverse Score%sComponent Overlap (Da)\n",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
  }
  else
  {
    fprintf(output,
            "Contig Idx1%sContig Idx2%sExperimental Shift FF%sExperimental Shift RR%sRemoved?%sContig2 Reversed?%sForward Score%sReverse Score%sComponent Overlap (Da)\n",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
  }

  bool checkOverlap = true;
  bool mergingMatePairs = false;
  int countmatep = 0, countcont = 0, countmerg = 0;

  scoredEdgeIt = scoredEdges.begin();
  while (scoredEdgeIt != scoredEdges.end())
  {
    int edgeIdx = floatToInt((*scoredEdgeIt)[0]);
    int origIdx = edgeIdx;
    edgeScore = (*scoredEdgeIt)[1];
    int idx1, idx2;

    if (edgeIdx >= 0)
    {
      while (edgeRef[edgeIdx] != -1)
        edgeIdx = edgeRef[edgeIdx];
      idx1 = floatToInt(edges[edgeIdx][0]);
      idx2 = floatToInt(edges[edgeIdx][1]);

      if (erasedEdges.count(edgeIdx) > 0 || components[idx1].count(idx2) > 0
          || components[idx2].count(idx1) > 0 || graph.count(idx1) == 0
          || graph[idx1].count(idx2) == 0)
      {
        scoredEdgeIt++;
        continue;
      }

    }
    else
    {
      mergingMatePairs = true;
      cout << "\n\nFinding Connector Edges\n";
      //outputComponents(output, graph, edges, components, oriented, reversed, false, false);

      fprintf(output, "\n\nAdding Connector Edges\n");

      if (haveRes)
      {
        fprintf(output,
                "Protein Idx1%sProtein Idx2%sContig Idx1%sContig Idx2%sExperimental Shift FF%sExperimental Shift RR%sExpected Shift FF%sExpected Shift RR%sEqual?%sRemoved?%sContig2 Reversed?%sContig1 Reversed Global?%sContig2 Reversed Global?%sForward Score%sReverse Score\n",
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP);
      }
      else
      {
        fprintf(output,
                "Contig Idx1%sContig Idx2%sExperimental Shift FF%sExperimental Shift RR%sRemoved?%sContig2 Reversed?%sForward Score%sReverse Score\n",
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP,
                CSV_SEP);
      }

      getCondensedAlignmentGraph(components,
                                 graph,
                                 edges,
                                 reversed,
                                 scoredEdges,
                                 root_alignments,
                                 precCharge,
                                 minRatioConnContig,
                                 minRatioConn,
                                 minNumMatchedPeaksConn,
                                 minNumConnectors,
                                 minEdgesToComponent);
      int oldSize = edgeRef.size();
      edgeRef.resize(edges.size());
      for (int i = oldSize; i < edgeRef.size(); i++)
        edgeRef[i] = -1;
      cout << "\nadding connector edges...\n";
      //printf("\nCombined %s %.2f of possible edges (%.0f/%.0f)\n", "%", (100.0*edgesAdded)/totalEdges, edgesAdded, totalEdges);
      scoredEdgeIt = scoredEdges.begin();
      checkOverlap = true;
      //minCombScore = 2.0;
      edgeIdx = floatToInt((*scoredEdgeIt)[0]);
      edgeScore = (*scoredEdgeIt)[1];
      idx1 = floatToInt(edges[edgeIdx][0]);
      idx2 = floatToInt(edges[edgeIdx][1]);
    }

    oldFFShift = edges[edgeIdx][2];
    oldFRShift = edges[edgeIdx][3];
    oldRFShift = edges[edgeIdx][4];
    oldRRShift = edges[edgeIdx][5];

    /**
     * Test merge here
     */

    /**
     *
     */

    float fShift = edges[edgeIdx][2], rShift = edges[edgeIdx][3];

    bool debug = false;//(idx1 == 1029 && idx2 == 137);

    if (debug)
    {
      int origIdx1 = floatToInt(originalEdges[edgeIdx][0]);
      int origIdx2 = floatToInt(originalEdges[edgeIdx][1]);
      float origFFShift = originalEdges[edgeIdx][2];
      float origFRShift = originalEdges[edgeIdx][3];
      cout << "\nMerging " << idx1 << " (" << origIdx1 << ") and " << idx2
          << " (" << origIdx2 << "), shift = " << origFFShift << ", "
          << origFRShift << "\n";
      cout << "Contig " << origIdx1 << " (" << reversed[origIdx1] << "):\n";
      contigs[origIdx1].output(cout);
      /*
       cout << "\nContig " << origIdx1 << " (reversed):\n";
       Spectrum tempContig1 = contigs[origIdx1];
       tempContig1.reverse(0.0 - AAJumps::massH2O, 0);
       tempContig1.output(cout);
       */
      cout << "\nContig " << origIdx2 << " (" << reversed[origIdx2] << "):\n";
      contigs[origIdx2].output(cout);
      cout << "\nContig " << origIdx2 << " (reversed):\n";
      Spectrum tempContig = contigs[origIdx2];
      tempContig.reverse(0.0 - AAJumps::massH2O, 0);
      tempContig.output(cout);
    }

    res = getConsensusOverlap(idx1,
                              idx2,
                              fShift,
                              rShift,
                              graph,
                              edges,
                              components,
                              root_alignments,
                              false);

    float avgF = res[0], avgR = res[1];

    bool rev = avgR > avgF;
    //rev = (prot_match[idx1][2] != prot_match[idx2][2]) && (reversed[idx1] == reversed[idx2]);
    bool remove = max(avgF, avgR) <= 0 || (checkOverlap && (max(avgF, avgR)
        < minCombScore));
    //remove = false;

    const char* revStr = (rev) ? "1" : "0";

    float shift_use = (rev) ? rShift : fShift;
    float overlap_area = getComponentOverlapArea(idx1,
                                                 idx2,
                                                 shift_use,
                                                 rev,
                                                 graph,
                                                 edges,
                                                 components);

    if (!remove)
    {
      //if (rev)
      //  reverseNode(idx2, graph, edges, components, reversed);
      remove = !tryMergeContigs(idx1,
                                idx2,
                                edgeIdx,
                                graph,
                                edges,
                                components,
                                root_alignments,
                                rev,
                                debug);
      //if (rev)
      //  reverseNode(idx2, graph, edges, components, reversed);
    }

    const char* remStr = (remove) ? "1" : "0";

    conFFShift = (rev) ? oldFRShift : oldFFShift;
    conFRShift = (rev) ? oldFFShift : oldFRShift;
    conRFShift = (rev) ? oldRRShift : oldRFShift;
    conRRShift = (rev) ? oldRFShift : oldRRShift;
    bool correct;
    if (haveRes && prot_match[idx1][0] >= 0 && prot_match[idx2][0] >= 0)
    {
      totalEdges += 1.0;
      bool reverse1 = prot_match[idx1][2] == 1;
      bool reverse2 = prot_match[idx2][2] == 1;
      float shift = conFFShift;
      float shiftR = conRRShift;

      bool sameProt = prot_match[idx1][0] == prot_match[idx2][0];
      float protShift1 = getContigShift(idx1, reverse1);
      float protShift2 = getContigShift(idx2, reverse2);
      float protFFShift = protShift2 - protShift1;
      float protRRShift = (protShift1 + contigs[idx1].parentMass) - (protShift2
          + contigs[idx2].parentMass);
      float temF = protFFShift, temR = protRRShift;
      protFFShift = (reverse1 == reversed[idx1]) ? temF : temR;
      protRRShift = (reverse1 == reversed[idx1]) ? temR : temF;
      correct = isEqual(shift, protFFShift, 2.0) && prot_match[idx1][0]
          == prot_match[idx2][0];
      if (!remove && !sameProt)
      {
        cout << "INCORRECTLY MERGING CONTIGS FROM DIFFERENT PROTEINS [contigs "
            << idx1 << " (prot " << prot_match[idx1][0] << ") and " << idx2
            << " (prot " << prot_match[idx2][0] << "), avgOverlap = "
            << max(avgF, avgR) << "]\n";
      }
      else if (!remove && !correct)
      {
        cout << "CONTIG-CONTIG SHIFT DIFFERS FROM EXPECTED BY "
            << abs(protFFShift - shift) << " [contigs " << idx1 << " and "
            << idx2 << " (prot " << prot_match[idx1][0] << "), avgOverlap = "
            << max(avgF, avgR) << "]\n";
        //cout << fShift << " or " << rShift << ": chose " << shift << "\n";
        //cout << prot_match[idx1][2] << ", " << reversed[idx1] << ": " << prot_match[idx2][2] << ", " << reversed[idx2] << "\n";
      }
      protFFShift = (sameProt) ? protFFShift : 0 / zero;
      protRRShift = (sameProt) ? protRRShift : 0 / zero;
      const char* corrStr = (sameProt && correct) ? "1" : "0";
      fprintf(output,
              "%d%s%d%s%d%s%d%s%.2f%s%.2f%s%.2f%s%.2f%s%s%s%s%s%s%s%d%s%d%s%.2f%s%.2f%s%.2f\n",
              prot_match[idx1][0],
              CSV_SEP,
              prot_match[idx2][0],
              CSV_SEP,
              idx1,
              CSV_SEP,
              idx2,
              CSV_SEP,
              shift,
              CSV_SEP,
              shiftR,
              CSV_SEP,
              protFFShift,
              CSV_SEP,
              protRRShift,
              CSV_SEP,
              corrStr,
              CSV_SEP,
              remStr,
              CSV_SEP,
              revStr,
              CSV_SEP,
              prot_match[idx1][2],
              CSV_SEP,
              prot_match[idx2][2],
              CSV_SEP,
              avgF,
              CSV_SEP,
              avgR,
              CSV_SEP,
              overlap_area);
    }
    else if (haveRes)
    {
      fprintf(output,
              "%d%s%d%s%d%s%d%s%.2f%s%.2f%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%.2f%s%.2f%s%.2f\n",
              prot_match[idx1][0],
              CSV_SEP,
              prot_match[idx2][0],
              CSV_SEP,
              idx1,
              CSV_SEP,
              idx2,
              CSV_SEP,
              conFFShift,
              CSV_SEP,
              conRRShift,
              CSV_SEP,
              "",
              CSV_SEP,
              "",
              CSV_SEP,
              "",
              CSV_SEP,
              remStr,
              CSV_SEP,
              revStr,
              CSV_SEP,
              "",
              CSV_SEP,
              "",
              CSV_SEP,
              avgF,
              CSV_SEP,
              avgR,
              CSV_SEP,
              overlap_area);
    }
    else
    {
      const char* remStr = (remove) ? "1" : "0";
      const char* revStr = (rev) ? "1" : "0";
      float shift = conFFShift;
      float shiftR = conRRShift;
      fprintf(output,
              "%d%s%d%s%.2f%s%.2f%s%s%s%s%s%.2f%s%.2f%s%.2f\n",
              idx1,
              CSV_SEP,
              idx2,
              CSV_SEP,
              shift,
              CSV_SEP,
              shiftR,
              CSV_SEP,
              remStr,
              CSV_SEP,
              revStr,
              CSV_SEP,
              avgF,
              CSV_SEP,
              avgR,
              CSV_SEP,
              overlap_area);
    }
    if (haveRes && correct)
      goodEdges += 1.0;
    if (!remove)
      edgesAdded += 1.0;
    if (haveRes && correct && !remove)
      goodEdgesAdded += 1.0;

    if (remove)
    {
      erasedEdges.insert(edgeIdx);
      if (graph[idx1][idx2] == edgeIdx)
      {
        erasedEdges.insert(graph[idx2][idx1]);
        graph[idx1].erase(idx2);
        graph[idx2].erase(idx1);
      }
      scoredEdgeIt++;
      if (debug)
        break;
      continue;
    }

    if (mergingMatePairs)
    {
      countmatep += edgeNum[edgeIdx];
      countmerg += 1;
    }
    else
    {
      countcont += 1;
    }

    edgeTo[0] = (float)idx1;
    edgeTo[1] = (float)idx2;
    edgeTo[2] = conFFShift;
    edgeTo[3] = conFRShift;
    edgeTo[4] = conRFShift;
    edgeTo[5] = conRRShift;
    edgeFrom[0] = (float)idx2;
    edgeFrom[1] = (float)idx1;
    edgeFrom[2] = 0.0 - conFFShift;
    edgeFrom[3] = 0.0 - conRFShift;
    edgeFrom[4] = 0.0 - conFRShift;
    edgeFrom[5] = 0.0 - conRRShift;
    edges[edgeIdx] = edgeTo;
    edges[graph[idx2][idx1]] = edgeFrom;

    erasedEdges.insert(graph[idx2][idx1]);
    graph[idx2].erase(idx1);

    if (rev)
      reverseNode(idx2, graph, edges, components, reversed);

    int most_matched_peaks = 0;
    TwoValues<int> best_node_pair;
    best_node_pair[0] = idx1;
    best_node_pair[1] = idx2;
    for (nodeIt = components[idx1].begin(); nodeIt != components[idx1].end(); nodeIt++)
    {
      float shift1 = (*nodeIt == idx1) ? 0 : edges[graph[idx1][*nodeIt]][2];
      for (nodeIt2 = components[idx2].begin(); nodeIt2
          != components[idx2].end(); nodeIt2++)
      {
        if (*nodeIt2 == *nodeIt)
          continue;
        float shift2 = (*nodeIt2 == idx2) ? 0 : edges[graph[idx2][*nodeIt2]][2];
        float locFShift = shift_use - shift1 + shift2;
        TwoValues<int> mp = getContigOverlapPeaks(oriented[*nodeIt],
                                                  oriented[*nodeIt2],
                                                  locFShift,
                                                  0,
                                                  parentMassTol,
                                                  minContigOverlap,
                                                  debug);
        if (mp[0] > most_matched_peaks)
        {
          best_node_pair[0] = *nodeIt;
          best_node_pair[1] = *nodeIt2;
          most_matched_peaks = mp[0];
        }
      }
    }
    if (most_matched_peaks == 0)
    {
      cout << "\nCOULD NOT FIND OVERLAPPING CONTIGS BETWEEN COMONENTS " << idx1
          << " AND " << idx2 << "\n";
    }
    if (debug)
    {
      cout << "\nbest node pair: " << best_node_pair[0] << " and "
          << best_node_pair[1] << "(" << most_matched_peaks
          << " matching peaks)\n";
    }
    //list<TwoValues<int> > idx2_aligns = root_alignments[idx2];
    //idx2_aligns.push_back(best_node_pair);
    root_alignments[idx1].push_back(best_node_pair);
    for (list<TwoValues<int> >::iterator alIt = root_alignments[idx2].begin(); alIt
        != root_alignments[idx2].end(); alIt++)
    {
      root_alignments[idx1].push_back(*alIt);
    }

    /*if (idx1 == 3247 || idx1 == 1029 || idx1 == 137 || idx1 == 1557 || idx1 == 1293) {
     cout << "idx1 = " << idx1 << "\n";
     cout << "idx2 = " << idx2 << "\n";
     cout << "components[idx1].size() = " << components[idx1].size() << "\n";
     cout << "components[idx2].size() = " << components[idx2].size() << "\n";
     cout << "root_alignments[idx1].size() = " << root_alignments[idx1].size() << "\n";
     cout << "root_alignments[idx2].size() = " << root_alignments[idx2].size() << "\n";
     }*/

    components[idx1].insert(idx2);

    for (nodeIt = components[idx2].begin(); nodeIt != components[idx2].end(); nodeIt++)
    {
      int node2 = *nodeIt;

      if (node2 == idx2)
        continue;

      int edge2 = graph[idx2][node2];

      FFShift = conFFShift + edges[edge2][2];
      FRShift = conFFShift + edges[edge2][3];
      RFShift = conRRShift + edges[edge2][4];
      RRShift = conRRShift + edges[edge2][5];

      edgeTo[0] = (float)idx1;
      edgeTo[1] = (float)node2;
      edgeTo[2] = FFShift;
      edgeTo[3] = FRShift;
      edgeTo[4] = RFShift;
      edgeTo[5] = RRShift;

      components[idx1].insert(node2);
      edges[edge2] = edgeTo;
      graph[idx1][node2] = edge2;
    }

    Spectrum consensusIdx1;
    pair<pair<vector<int> , vector<int> > , vector<pair<vector<int> , vector<
        double> > > > abIdx1;
    TwoValues<float> extraShift1 = mergeContigs(idx1,
                                                graph,
                                                edges,
                                                components,
                                                root_alignments,
                                                consensusIdx1,
                                                abIdx1,
                                                debug);


    consensus[idx1] = consensusIdx1;
    extraShifts[idx1] = extraShift1;
    parent_abinfo[idx1] = abIdx1;
    contigsMerged[idx1] = true;
    contigsMerged[idx2] = true;

    map<int, map<int, int> > graphSep = graph;
    vector<vector<float> > edgesSep = edges;

    for (mapIt = graph[idx2].begin(); mapIt != graph[idx2].end(); mapIt++)
    {
      int node2 = mapIt->first;
      int edge2 = mapIt->second;
      if (components[idx1].count(node2) > 0)
        continue;

      FFShift = conFFShift + edges[edge2][2];
      FRShift = conFFShift + edges[edge2][3];
      RFShift = conRRShift + edges[edge2][4];
      RRShift = conRRShift + edges[edge2][5];

      edgeTo[0] = (float)idx1;
      edgeTo[1] = (float)node2;
      edgeTo[2] = FFShift;
      edgeTo[3] = FRShift;
      edgeTo[4] = RFShift;
      edgeTo[5] = RRShift;
      edgeFrom[0] = (float)node2;
      edgeFrom[1] = (float)idx1;
      edgeFrom[2] = 0.0 - FFShift;
      edgeFrom[3] = 0.0 - RFShift;
      edgeFrom[4] = 0.0 - FRShift;
      edgeFrom[5] = 0.0 - RRShift;

      if (graph[idx1].count(node2) > 0)
      {
        bool edge1OK = graph[idx1][node2] <= maxContigEdgeIdx;
        bool edge2OK = graph[idx2][node2] <= maxContigEdgeIdx;

        float prevShiftF = edges[graph[idx1][node2]][2];
        float prevShiftR = edges[graph[idx1][node2]][3];
        float prevShiftRF = edges[graph[idx1][node2]][4];
        float prevShiftRR = edges[graph[idx1][node2]][5];
        if (!isEqual(prevShiftF, FFShift, 2.0) && !isEqual(prevShiftR,
                                                           FRShift,
                                                           2.0)
            && !isEqual(prevShiftRF, RFShift, 2.0) && !isEqual(prevShiftRR,
                                                               RRShift,
                                                               2.0))
        {
          bool deb = false;//(idx1 == 253 && idx2 == 523);
          //if (deb) {cout << "idx1 = " << idx1 << ", idx2 = " << idx2 << ", node2 = " << node2 << "\n"; cout.flush();}
          res1 = getConsensusOverlap(idx1,
                                     node2,
                                     prevShiftF,
                                     prevShiftR,
                                     graphSep,
                                     edgesSep,
                                     components,
                                     root_alignments,
                                     deb);
          res2 = getConsensusOverlap(idx1,
                                     node2,
                                     FFShift,
                                     FRShift,
                                     graphSep,
                                     edgesSep,
                                     components,
                                     root_alignments,
                                     deb);
          float score1 = max(res1[0], res1[1]);
          float score2 = max(res2[0], res2[1]);

          //cout << "Conflicting edges between " << node2 << " and component (" << idx1 << "," << idx2 << ")\n"; cout << idx1 << " -- " << node2 << ": FFShift = " << prevShiftF << ", FRShift = " << prevShiftR << ", RFShift = " << prevShiftRF << ", RRShift = " << prevShiftRR << "\n" << idx1 << " -- " << idx2 << " -- " << node2 << ": FFShift = " << FFShift << ", FRShift = " << FRShift << ", RFShift = " << RFShift << ", RRShift = " << RRShift << "\n";
          /*
           bool r1 = prot_match[idx1][2] == 1;
           bool r2 = prot_match[node2][2] == 1;
           float ps1 = getContigShift(idx1, r1);
           float ps2f = getContigShift(node2, r2);
           cout << "correct shift = " << ps2f - ps1 << "\n";
           */
          if ((edge1OK && !edge2OK) || (score1 >= score2))
          {
            //cout << "Removing edge from " << idx2 << " to " << node2 << "\n";
            erasedEdges.insert(graph[idx2][node2]);
            erasedEdges.insert(graph[node2][idx2]);
            edgeRef[graph[idx2][node2]] = graph[idx1][node2];
            edgeRef[graph[node2][idx2]] = graph[node2][idx1];
            graph[node2].erase(idx2);
          }
          else
          {//score2 > minSecondaryRatio &&
            //cout << "Removing edge from " << idx1 << " to " << node2 << "\n";
            erasedEdges.insert(graph[idx1][node2]);
            erasedEdges.insert(graph[node2][idx1]);
            edgeRef[graph[idx1][node2]] = graph[idx2][node2];
            edgeRef[graph[node2][idx1]] = graph[node2][idx2];
            //if (components[idx2].count(node2) > 0) idx1sToAdd.insert(node2);//components[idx1].insert(node2);
            int edgeToNode = graph[idx2][node2];
            int edgeFromNode = graph[node2][idx2];
            edges[edgeToNode] = edgeTo;
            edges[edgeFromNode] = edgeFrom;
            graph[node2].erase(idx2);
            graph[node2][idx1] = edgeFromNode;
            graph[idx1][node2] = edgeToNode;
          }/* else {
           //cout << "Removing edge from " << idx2 << " to " << node2 << "\n";
           //cout << "Removing edge from " << idx1 << " to " << node2 << "\n";
           erasedEdges.insert(graph[idx2][node2]);
           erasedEdges.insert(graph[node2][idx2]);
           graph[node2].erase(idx2);
           erasedEdges.insert(graph[idx1][node2]);
           erasedEdges.insert(graph[node2][idx1]);
           graph[idx1].erase(node2);
           graph[node2].erase(idx1);
           }*/
        }
        else
        {
          //cout << "Agreeable edges between " << node2 << " and component (" << idx1 << "," << idx2 << ")\n";
          /*if (components[idx2].count(node2) > 0) {
           //cout << idx2 << " already in component.\n";
           edges[graph[idx2][node2]] = edgeTo;
           idx1sToAdd.insert(node2);//components[idx1].insert(node2);
           }
           else if (components[idx1].count(node2) > 0) {
           //cout << idx1 << " already in component.\n";
           erasedEdges.insert(graph[idx2][node2]);
           erasedEdges.insert(graph[node2][idx2]);
           graph[node2].erase(idx2);
           }*/
          if (edge2OK && !edge1OK)
          {
            //cout << "choosing edge " << edge2 << " 2( " << maxContigEdgeIdx << " ) ( " << edgeTo[2] << " )\n";
            erasedEdges.insert(graph[idx1][node2]);
            erasedEdges.insert(graph[node2][idx1]);
            edgeRef[graph[idx1][node2]] = graph[idx2][node2];
            edgeRef[graph[node2][idx1]] = graph[node2][idx2];
            //if (components[idx2].count(node2) > 0) idx1sToAdd.insert(node2);//components[idx1].insert(node2);
            int edgeToNode = graph[idx2][node2];
            int edgeFromNode = graph[node2][idx2];
            edges[edgeToNode] = edgeTo;
            edges[edgeFromNode] = edgeFrom;
            graph[node2].erase(idx2);
            graph[node2][idx1] = edgeFromNode;
            graph[idx1][node2] = edgeToNode;
          }
          else
          {
            //cout << "choosing edge " << graph[idx1][node2] << " 1( " << maxContigEdgeIdx << " ) ( " << edges[graph[idx1][node2]][2] << " )\n";
            erasedEdges.insert(graph[idx2][node2]);
            erasedEdges.insert(graph[node2][idx2]);
            edgeRef[graph[idx2][node2]] = graph[idx1][node2];
            edgeRef[graph[node2][idx2]] = graph[node2][idx1];
            graph[node2].erase(idx2);
          }
        }

      }
      else
      {
        edges[edge2] = edgeTo;
        graph[idx1][node2] = edge2;

        int node2Edge = graph[node2][idx2];
        edges[node2Edge] = edgeFrom;
        graph[node2].erase(idx2);
        graph[node2][idx1] = node2Edge;
      }
    }
    //cout << "erasing\n";
    //cout.flush();
    graph.erase(idx2);
    components[idx2].clear();
    consensus[idx2].resize(0);
    parent_abinfo.erase(idx2);

    //outputGraph(cout, graph, edges, &erasedEdges, 0);
    //break;
    //if (idx1 == 90 && idx2 == 153) {cout << "\n" << idx1 << " and " << idx2 << " after: " << edgeIdx << "\n\n"; outputGraph(cout, graph, edges, &erasedEdges, &components);}
    if (scoredEdgeIt != scoredEdges.end())
      scoredEdgeIt++;
    scoredEdges.erase(scoredEdges.begin(), scoredEdgeIt);

    recomputeEdges(idx1,
                   graph,
                   edges,
                   components,
                   erasedEdges,
                   edgeRef,
                   root_alignments,
                   scoredEdges);

    scoredEdgeIt = scoredEdges.begin();

    if (debug)
      break;
  }

  set<int> countNodes;
  vector<int> compSize(components.size());
  for (int i = 0; i < components.size(); i++)
  {
    compSize[i] = components[i].size();
    if (components[i].size() == 0)
    {
      consensus[i].resize(0);
      continue;
    }
    if (components[i].size() < 1)
      continue;
    set<int> seen;
    for (nodeIt = components[i].begin(); nodeIt != components[i].end(); nodeIt++)
    {
      seen. insert(*nodeIt);
    }
    //c  ompSize[i] = seen.size();
    if (mergeType == 1 && seen.size() < minComp)
    {
      components[i].clear();
      continue;
    }
    countNodes.insert(i);
    for (nodeIt = components[i].begin(); nodeIt != components[i].end(); nodeIt++)
    {
      countNodes.insert(*nodeIt);
    }
  }
  //printf("\nCombined %s %.2f good edges (%.0f/%.0f)\n", "%", (100.0*goodEdges)/edgesAdded, goodEdges, edgesAdded);
  outputComponents(output,
                   graph,
                   edges,
                   components,
                   oriented,
                   reversed,
                   false,
                   false,
                   0,
                   true);
  fclose(output);

  if (mergeType == 1)
    printf("\nStatistics for merged contigs after contig/mate-pair shifts used:\n");
  else
    printf("\nStatistics for merged contigs after contig/contig shifts used:\n");

  outputComponents(0,
                   graph,
                   edges,
                   components,
                   oriented,
                   reversed,
                   true,
                   false,
                   &countNodes);

  printf("\nStatistics for contigs before any shifts used:\n");
  outputComponents(0,
                   graphBefMerg,
                   edgesBefMerg,
                   componentsBefMerg,
                   myContigsBefMerg,
                   reversedBefMerg,
                   true,
                   false,
                   &countNodes,
                   true);

  printf("\nUsed %d contig/contig alignments and %d contig/mate-pair alignments (%d contigs merged with mate-pairs)\n\n",
         countcont,
         countmatep,
         countmerg);

  cout << "merging_params->min_component_size: "
      << merging_params->min_component_size << "\n";
  for (int i = 0; i < consensus.size(); i++)
  {
    if (compSize[i] < merging_params->min_component_size)
    {
      consensus[i].resize(0);
    }
  }

  merging_params->consensus.resize(consensus.size());
  for (int i = 0; i < consensus.size(); i++)
  {
    merging_params->consensus[i] = consensus[i];
  }

  SpecSet stars = *(merging_params->star_spectra);
  stars.addZPMpeaks(peakTol, 0, false);

  /*
   unsigned int idxcheck = 11;
   cout << "\n" <<idxcheck << ":\n";
   for (int i = 0; i < parent_abinfo[idxcheck].second.size(); i++)
   {
   cout << " -- " << i << " -> ";
   for (int j = 0; j < parent_abinfo[idxcheck].second[i].first.size(); j++)
   {
   cout << "(" << parent_abinfo[idxcheck].second[i].first[j] << ", " << reversed[parent_abinfo[idxcheck].second[i].first[j]] << ", " << parent_abinfo[idxcheck].second[i].second[j] << "); ";
   }
   cout << endl;
   }
   */

  Combine_abinfo_v1_0(contig_abinfo,
                      reversed,
                      stars,
                      contigsMerged,
                      contigs,
                      parent_abinfo,
                      consensus_abinfo);

  merging_params->out_reversed.resize(reversed.size());

  for (int i = 0; i < reversed.size(); i++)
    merging_params->out_reversed[i] = (reversed[i]) ? 1 : 0;

  merging_params->out_abinfo.clear();
  for (unsigned int i = 0; i < consensus.size(); i++)
  {
    merging_params->out_abinfo[i] = consensus_abinfo[i];
  }

}

bool CombineContigs::saveComponents(const char* outcomponents)
{
  FILE* compsOut = fopen(outcomponents, "w");
  if (compsOut == NULL)
  {
    cerr << "ERROR: Cannot write to file " << outcomponents << endl;
    return false;
  }

  int count = 0;
  fprintf(compsOut,
          "Component Idx%sContig Idx2%sContig Root Idx%sContig2 Shift\n",
          CSV_SEP,
          CSV_SEP,
          CSV_SEP);
  for (int i = 0; i < components.size(); i++)
  {
    if (components[i].size() < 2)
      continue;
    set<int> seen;
    map<int, set<int> > seenNodes;
    set<int> snode;
    for (set<int>::iterator nodeIt = components[i].begin(); nodeIt
        != components[i].end(); nodeIt++)
    {
      if (seen.count(*nodeIt) > 0)
        continue;
      seen.insert(*nodeIt);
      float FFShift = (i == *nodeIt) ? 0.0 : edges[graph[i][*nodeIt]][2];
      fprintf(compsOut,
              "%d%s%d%s%d%s%.1f\n",
              count,
              CSV_SEP,
              *nodeIt,
              CSV_SEP,
              i,
              CSV_SEP,
              FFShift);
    }
    count++;
  }
  fclose(compsOut);
  return true;
}

bool CombineContigs::getContigShiftStats(const char* outfile,
                                         int startMinNumMatchedPeaks,
                                         int endMinNumMatchedPeaks,
                                         int stepMinNumMatchedPeaks,
                                         float startMinRatio,
                                         float endMinRatio,
                                         float stepMinRatio,
                                         float startScore,
                                         float endScore,
                                         float stepScore)
{
  if (!haveRes)
  {
    cerr
        << "ERROR: Cannot output statistics if matchma results are not loaded\n";
    return false;
  }
  list<Results_OCC>::iterator pivot_res;
  int tp, fn, fp, t;
  FILE* output = fopen(outfile, "w");

  if (output == NULL)
  {
    cerr << "ERROR: Cannot write to file " << outfile << endl;
    return false;
  }

  fprintf(output,
          "Min Num Matched Peaks * Min Ratio%sTP%sFN%sFP%sNum Above Threshhold%sTP/(TP+FP)%sTP/(TP+FN)\n",
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP);
  for (float score = startScore; score < (endScore + stepScore); score
      += stepScore)
  {
    tp = 0;
    fn = 0;
    fp = 0;
    t = 0;
    for (pivot_res = contig_alignments.begin(); pivot_res
        != contig_alignments.end(); pivot_res++)
    {
      int idx1 = (*pivot_res).spec1, idx2 = (*pivot_res).spec2, mPeaks =
          (*pivot_res).matchedPeaks;
      if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0)
        continue;
      float shift = (*pivot_res).shift1, mRatio = min((*pivot_res).score1,
                                                      (*pivot_res).score2);
      float combScore = ((float)mPeaks) * mRatio;
      bool reverse2 = isEqual((*pivot_res).shift2, 1.0, 0.0001);
      bool thresh = combScore >= score;
      bool correct = validShift(shift, idx1, idx2, reverse2);

      if (thresh)
        t++;
      if (correct && thresh)
        tp++;
      else if (correct && !thresh)
        fn++;
      else if (!correct && thresh)
        fp++;

      /*if (minMatchedPeaks == endMinNumMatchedPeaks-1 && minRatio == startMinRatio && thresh && ! correct) {
       cout << "false-positive alignment between " << idx1 << " (protein " << prot_match[idx1][0] << ") and " << idx2 << " (protein " << prot_match[idx2][0] << ") with " << mPeaks << " matched peaks and " << minRatio << " min ratio\n";
       }*/
    }
    fprintf(output,
            "%.3f%s%d%s%d%s%d%s%d%s%.3f%s%.3f\n",
            score,
            CSV_SEP,
            tp,
            CSV_SEP,
            fn,
            CSV_SEP,
            fp,
            CSV_SEP,
            t,
            CSV_SEP,
            ((float)tp) / ((float)tp + (float)fp),
            CSV_SEP,
            ((float)tp) / ((float)tp + (float)fn));
  }

  fprintf(output,
          "\nMin Num Matched Peaks%sMin Ratio%sTP%sFN%sFP%sNum Above Threshhold%sTP/(TP+FP)%sTP/(TP+FN)\n",
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP);
  for (int minMatchedPeaks = startMinNumMatchedPeaks; minMatchedPeaks
      <= endMinNumMatchedPeaks; minMatchedPeaks += stepMinNumMatchedPeaks)
  {
    for (float minRatio = startMinRatio; minRatio
        < (endMinRatio + stepMinRatio); minRatio += stepMinRatio)
    {
      tp = 0;
      fn = 0;
      fp = 0;
      t = 0;
      for (pivot_res = contig_alignments.begin(); pivot_res
          != contig_alignments.end(); pivot_res++)
      {
        int idx1 = (*pivot_res).spec1, idx2 = (*pivot_res).spec2, mPeaks =
            (*pivot_res).matchedPeaks;
        if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0)
          continue;
        float shift = (*pivot_res).shift1, mRatio = min((*pivot_res).score1,
                                                        (*pivot_res).score2);
        bool reverse2 = isEqual((*pivot_res).shift2, 1.0, 0.0001);
        bool thresh = mRatio >= minRatio && mPeaks >= minMatchedPeaks;
        bool correct = validShift(shift, idx1, idx2, reverse2);

        if (thresh)
          t++;
        if (correct && thresh)
          tp++;
        else if (correct && !thresh)
          fn++;
        else if (!correct && thresh)
          fp++;

        /*if (minMatchedPeaks == endMinNumMatchedPeaks-1 && minRatio == startMinRatio && thresh && ! correct) {
         cout << "false-positive alignment between " << idx1 << " (protein " << prot_match[idx1][0] << ") and " << idx2 << " (protein " << prot_match[idx2][0] << ") with " << mPeaks << " matched peaks and " << minRatio << " min ratio\n";
         }*/
      }
      fprintf(output,
              "%d%s%.3f%s%d%s%d%s%d%s%d%s%.3f%s%.3f\n",
              minMatchedPeaks,
              CSV_SEP,
              minRatio,
              CSV_SEP,
              tp,
              CSV_SEP,
              fn,
              CSV_SEP,
              fp,
              CSV_SEP,
              t,
              CSV_SEP,
              ((float)tp) / ((float)tp + (float)fp),
              CSV_SEP,
              ((float)tp) / ((float)tp + (float)fn));
    }
  }

  TwoValues<int> mp, mp2;
  TwoValues<float> sc;
  fprintf(output,
          "\nMin Consecutive Matched Peaks%sMin Ratio%sTP%sFN%sFP%sNum Above Threshhold%sTP/(TP+FP)%sTP/(TP+FN)\n",
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP);
  for (int minMatchedPeaks = startMinNumMatchedPeaks; minMatchedPeaks
      <= endMinNumMatchedPeaks; minMatchedPeaks += stepMinNumMatchedPeaks)
  {
    for (float minRatio = startMinRatio; minRatio
        < (endMinRatio + stepMinRatio); minRatio += stepMinRatio)
    {
      tp = 0;
      fn = 0;
      fp = 0;
      t = 0;
      for (pivot_res = contig_alignments.begin(); pivot_res
          != contig_alignments.end(); pivot_res++)
      {
        int idx1 = (*pivot_res).spec1, idx2 = (*pivot_res).spec2, mPeaks =
            (*pivot_res).matchedPeaks;
        if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0
            || idxNotUsed.count(prot_match[idx2][0]) > 0
            || idxNotUsed.count(prot_match[idx1][0]) > 0)
          continue;
        float shift = (*pivot_res).shift1, mRatio = min((*pivot_res).score1,
                                                        (*pivot_res).score2);
        bool reverse2 = isEqual((*pivot_res).shift2, 1.0, 0.0001);

        Spectrum contig2 = contigs[idx2];
        if (reverse2)
        {
          contig2.reverse(0.0 - AAJumps::massH2O, 0);
        }
        bool reverse1 = prot_match[idx1][2] == 1;
        bool reverse2h = prot_match[idx2][2] == 1;
        float protShift1 = getContigShift(idx1, reverse1);
        float protShift2 = getContigShift(idx2, reverse2h);
        float protFFShift = protShift2 - protShift1;
        //float protRRShift = (protShift1+contigs[idx1].parentMass) - (protShift2+contigs[idx2].parentMass);
        //protFFShift = (reverse1) ? protRRShift : protFFShift;
        //if (protFFShift < 0 - contigs[idx2].parentMass || protFFShift > contigs[idx1].parentMass) continue;
        Spectrum contig1 = contigs[idx1];
        if (reverse1)
        {
          contig1.reverse(0.0 - AAJumps::massH2O, 0);
        }
        Spectrum contig2h = contigs[idx2];
        if (reverse2h)
        {
          contig2h.reverse(0.0 - AAJumps::massH2O, 0);
        }
        TwoValues<int> mp = getContigOverlapPeaks(contig1,
                                                  contig2h,
                                                  protFFShift,
                                                  contigs[idx1].parentMass,
                                                  0.5,
                                                  0,
                                                  false);
        //if (! validShift(protFFShift, idx1, idx2, reverse2)) {cout << "e";}
        //if (seenNodes.count(idx1) > 0 && seenNodes[idx1].count(idx2) > 0 && mp[0] < 4) cout << "\n" << idx1 << " (rev " << prot_match[idx1][2] << " prot " << prot_match[idx1][0] << ")to " << idx2 << " (rev " << prot_match[idx2][2] << " prot " << prot_match[idx2][0] << ")\n";
        if (mp[0] < 7)
          continue;

        mp = countConsecutiveMP(contigs[idx1],
                                contig2,
                                shift,
                                contigs[idx1].parentMass,
                                1.0,
                                minContigOverlap,
                                false);
        mPeaks = mp[0];
        bool thresh = mRatio >= minRatio && mPeaks >= minMatchedPeaks;
        bool correct = validShift(shift, idx1, idx2, reverse2);

        if (thresh)
          t++;
        if (correct && thresh)
          tp++;
        else if (correct && !thresh)
          fn++;
        else if (!correct && thresh)
          fp++;

        //if (minMatchedPeaks == endMinNumMatchedPeaks && minRatio == startMinRatio && thresh && ! correct) {
        //	cout << "false-positive alignment between " << idx1 << " (protein " << prot_match[idx1][0] << ") and " << idx2 << " (protein " << prot_match[idx2][0] << ") with " << mPeaks << " consecutive matched peaks and " << minRatio << " min ratio\n";
        //}
      }
      fprintf(output,
              "%d%s%.3f%s%d%s%d%s%d%s%d%s%.3f%s%.3f\n",
              minMatchedPeaks,
              CSV_SEP,
              minRatio,
              CSV_SEP,
              tp,
              CSV_SEP,
              fn,
              CSV_SEP,
              fp,
              CSV_SEP,
              t,
              CSV_SEP,
              ((float)tp) / ((float)tp + (float)fp),
              CSV_SEP,
              ((float)tp) / ((float)tp + (float)fn));
    }
  }
  /*
   fprintf(output, "\ nContig Index 1%sContig Index 2%sCorrect%sExpected Shift%sObserved Shift%sDifference%sContig 2 Reversed WRT Shift%sContig 1 Reversed Gloabl%sContig 2 Reversed Global%sMatched Peaks%sConsecutive Matching Peaks%sOverlapping Ratio%sProtein 1%sProtein 2\n", CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP);
   for (pivot_res = contig_alignments. begin(); pivot_res != contig_alignments.end(); pivot_res ++)
   {
   int idx1 = (*pivot_res).spec1, idx2 = (*pivot_res).spec2, mPeaks = (*pivot_res).matchedPeaks;
   if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0) continue;
   float shift = (*pivot_res).shift1, mRatio = min((*pivot_res).score1, (*pivot_res).score2);
   bool isreverse2 = isEqual((*pivot_res).shift2, 1.0, 0.0001);
   Spectrum contig2 = contigs[idx2];
   if (isreverse2)
   { contig2.reverse(0.0 - AAJumps::massH2O, 0);}
   bool debug = false;//(idx1 == 1 && idx2 == 634);
   mp = countConsecutiveMP(contigs[idx1], contig2, shift, contigs[idx1].parentMass, parentMassTol, minContigOverlap, false);
   mp2 = getContigOverlapPeaks(contigs[idx1], contig2, shift, contigs[idx1].parentMass, parentMassTol, minContigOverlap, false);
   sc = getContigOverlapScores(contigs[idx1], contig2, shift, contigs[idx1].parentMass, parentMassTol, minContigOverlap, debug);
   bool reverse1 = prot_match[idx1][2] == 1;
   bool reverse2 = prot_match[idx2][2] == 1;
   float protShift1 = getContigShift(idx1, reverse1);
   float protShi ft2 = getContigShift(idx2, reverse2);
   float protFFShift = protShift2 - protShift1;
   float protRRShift = (protShift1+contigs[idx1].parentMass) - (protShift2+contigs[idx2].parentMass);
   protFFShift = (reverse2 == isreverse2) ? protFFShift : protRRShift;
   bool correct = isEqual(shift, protFFShift, 2.0) && prot_match[idx1][0] == prot_match[idx2][0];
   int cor = (correct) ? 1 : 0;
   int isrev = (isreverse2) ? 1 : 0;
   int rev1 = (reverse1) ? 1 : 0;
   int rev2 = (reverse2) ? 1 : 0;
   float dist = abs(protFFShift - shift);
   fprintf(output, "%d%s%d%s%d%s%.2f%s%.2f%s%.2f%s%d%s%d%s%d%s%d%s%d%s%.2f%s%d%s%d\n", idx1, CSV_SEP, idx2, CSV_SEP, cor, CSV_SEP, protFFShift, CSV_SEP, shift, CSV_SEP, dist, CSV_SEP, isrev, CSV_SEP, rev1, CSV_SEP, rev2, CSV_SEP, mPeaks, CSV_SEP, mp[ 0], CSV_SEP, mRatio, CSV_SEP, prot_match[idx1][0], CSV_SEP, prot_match[idx2][0]);
   }*/
  fclose(output);

  return true;
}

void CombineContigs::init()
{
  merging_params->consensus = contigs;
  merging_params->out_abinfo = contig_abinfo;
  (merging_params->out_reversed).resize(contigs.size());
  for (unsigned int i = 0; i < (merging_params->out_reversed).size(); i++)
  {
    (merging_params->out_reversed)[i] = 0;
  }
  consensus = contigs;
  oriented = contigs;
  consensus_abinfo = contig_abinfo;
  extraShifts.resize(contigs.size());
  for (unsigned int i = 0; i < extraShifts.size(); i++)
  {
    extraShifts[i][0] = 0;
    extraShifts[i][1] = 0;
  }
  numCorrectPairs = 0;
}

void CombineContigs::outputConnectorOverlap(int connIdx,
                                            int contig1,
                                            int contig2,
                                            Results_OCC& res1,
                                            Results_OCC& res2)
{
  TwoValues<float> res;
  cout << "Overlap of connector " << connIdx << " and contig " << contig1
      << ":\n";
  res = getContigOverlapScores(connectors[connIdx],
                               contigs[contig1],
                               res1.shift1,
                               res1.shift2,
                               parentMassTol,
                               0,
                               true);

  cout << "\nOverlap of connector " << connIdx << " and contig " << contig2
      << ":\n";
  res = getContigOverlapScores(connectors[connIdx],
                               contigs[contig2],
                               res2.shift1,
                               res2.shift2,
                               parentMassTol,
                               0,
                               true);
}

void CombineContigs::getConnectorShiftStats(const char* outfile,
                                            int precCharge,
                                            int startMinNumMatchedPeaks,
                                            int endMinNumMatchedPeaks,
                                            int stepMinNumMatchedPeaks,
                                            float startMinRatio,
                                            float endMinRatio,
                                            float stepMinRatio)
{
  if (!haveRes)
    return;

  map<unsigned int, map<unsigned int, Results_OCC> > connContRes;
  map<unsigned int, map<unsigned int, list<unsigned int> > > graphCont;
  map<unsigned int, list<unsigned int> >::iterator graphContIt;
  map<unsigned int, set<unsigned int> > used;
  list<unsigned int>::iterator nodeIt;
  map<int, float> shiftScore;
  map<int, float> shiftScoreBT;
  map<int, float>::iterator shiftScoreIt;
  TwoValues<float> res;
  TwoValues<int> mp;

  getAlignmentGraph(&graphCont, &connContRes, precCharge, 2, 0.0, 0.6);
  int tp, fn, fp, t;
  FILE* output = fopen(outfile, "w");
  fprintf(output,
          "Min Num MatchedPeaks%smin Ratio%sTP%sFN%sFP%sNum Above Threshhold%sTP/(TP+FP)%sTP/(TP+FN)\n",
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP);
  for (int minMatchedPeaks = startMinNumMatchedPeaks; minMatchedPeaks
      <= endMinNumMatchedPeaks; minMatchedPeaks += stepMinNumMatchedPeaks)
  {
    cout << "Computing stats for " << minMatchedPeaks
        << " min matched peaks ...\n";
    cout.flush();
    for (float minRatio = startMinRatio; minRatio
        < (endMinRatio + stepMinRatio); minRatio += stepMinRatio)
    {
      tp = 0;
      fn = 0;
      fp = 0;
      t = 0;
      used.clear();
      for (int idx1 = 0; idx1 < contigs.size(); idx1++)
      {
        if (graphCont.count(idx1) == 0)
          continue;
        for (graphContIt = graphCont[idx1].begin(); graphContIt
            != graphCont[idx1].end(); graphContIt++)
        {
          int idx2 = graphContIt->first;
          if (used.count(idx1) > 0 && used[idx1].count(idx2) > 0)
            continue;
          if (used.count(idx2) > 0 && used[idx2].count(idx1) > 0)
            continue;
          if (used.count(idx1) == 0)
          {
            set<unsigned int> myS;
            myS.insert(idx2);
            used[idx1] = myS;
          }
          else
            used[idx1].insert(idx2);

          if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0
              || idxNotUsed.count(prot_match[idx1][0]) > 0
              || idxNotUsed.count(prot_match[idx2][0]) > 0)
          {
            continue;
          }

          bool reverse1 = prot_match[idx1][2] == 1;
          bool reverse2 = prot_match[idx2][2] == 1;
          float protShift1 = getContigShift(idx1, reverse1);
          float protShift2 = getContigShift(idx2, reverse2);
          float protFFShift = protShift2 - protShift1;
          //float protRRShift = (protShift1+contigs[idx1].parentMass) - (protShift2+contigs[idx2].parentMass);
          //protFFShift = (reverse1) ? protRRShift : protFFShift;
          //if (protFFShift < 0 - contigs[idx2].parentMass || protFFShift > contigs[idx1].parentMass) continue;
          /*
           Spectrum contig1 = contigs[idx1];
           if (reverse1) {contig1.reverse(0.0 - AAJumps::massH2O, 0);}
           Spectrum contig2 = contigs[idx2];
           if (reverse2) {contig2.reverse(0.0 - AAJumps::massH2O, 0);}
           TwoValues<int> mp = getContigOverlapPeaks(contig1, contig2, protFFShift, contigs[idx1].parentMass, 0.5, 0, false);
           if (mp[0] < 7) continue;
           */
          //if (! validShift(protFFShift, idx1, idx2, reverse2)) {cout << "e";}
          //if (seenNodes.count(idx1) > 0 && seenNodes[idx1].count(idx2) > 0 && mp[0] < 4) cout << "\n" << idx1 << " (rev " << prot_match[idx1][2] << " prot " << prot_match[idx1][0] << ") to " << idx2 << " (rev " << prot_match[idx2][2] << " prot " << prot_match[idx2][0] << ")\n";


          shiftScore.clear();
          shiftScoreBT.clear();
          float bestShift, bestShiftBT;
          float score;
          int numShifts;
          for (nodeIt = (graphContIt->second).begin(); nodeIt
              != (graphContIt->second).end(); nodeIt++)
          {
            Results_OCC res1 = connContRes[*nodeIt][idx1];
            Results_OCC res2 = connContRes[*nodeIt][idx2];
            int connIdx = *nodeIt;
            float minScore = min(res1.score1, res1.score2);
            minScore = min(minScore, res2.score1);
            minScore = min(minScore, res2.score2);
            int mPeaks = min(res1.matchedPeaks, res2.matchedPeaks);
            float fShift = res2.shift1 - res1.shift1;
            bestShift = fShift;
            score = minScore;
            int intShift = floatToInt(fShift / InputParams::Resolution);
            for (int s = intShift - intParentMassTol; s <= intShift
                + intParentMassTol; s++)
            {
              if (minScore < minRatio || mPeaks < minMatchedPeaks)
              {
                if (shiftScoreBT.count(s) == 0)
                {
                  shiftScoreBT[s] = minScore - (0.0001 * (float)abs(s
                      - intShift));
                }
                else
                {
                  shiftScoreBT[s] += minScore - (0.0001 * (float)abs(s
                      - intShift));
                }
              }
              else
              {
                if (shiftScore.count(s) == 0)
                {
                  shiftScore[s] = minScore
                      - (0.0001 * (float)abs(s - intShift));
                }
                else
                {
                  shiftScore[s] += minScore - (0.0001
                      * (float)abs(s - intShift));
                }
              }
            }
          }

          score = 0;
          for (shiftScoreIt = shiftScore.begin(); shiftScoreIt
              != shiftScore.end(); shiftScoreIt++)
          {
            float shiftScore = shiftScoreIt->second;
            if (shiftScore > score)
            {
              score = shiftScore;
              bestShift = ((float)shiftScoreIt->first)
                  * InputParams::Resolution;
            }
          }
          bool first = true;

          float FFShift, FRShift, dist;
          for (nodeIt = (graphContIt->second).begin(); nodeIt
              != (graphContIt->second).end(); nodeIt++)
          {
            Results_OCC res1 = connContRes[*nodeIt][idx1];
            Results_OCC res2 = connContRes[*nodeIt][idx2];
            float minScore = min(res1.score1, res1.score2);
            minScore = min(minScore, res2.score1);
            minScore = min(minScore, res2.score2);
            int mPeaks = min(res1.matchedPeaks, res2.matchedPeaks);
            if (minScore < minRatio || mPeaks < minMatchedPeaks)
              continue;
            float fShift = res2.shift1 - res1.shift1;
            if (first || abs(fShift - bestShift) < dist)
            {
              FFShift = fShift;
              FRShift = res2.shift2 - res1.shift1;
              dist = abs(fShift - bestShift);
              first = false;
            }
          }
          float obshift = (prot_match[idx1][2] == prot_match[idx2][2])
              ? FFShift : FRShift;

          score = 0;
          for (shiftScoreIt = shiftScoreBT.begin(); shiftScoreIt
              != shiftScoreBT.end(); shiftScoreIt++)
          {
            float shiftScore = shiftScoreIt->second;
            if (shiftScore > score)
            {
              score = shiftScore;
              bestShift = ((float)shiftScoreIt->first)
                  * InputParams::Resolution;
            }
          }
          first = true;
          for (nodeIt = (graphContIt->second).begin(); nodeIt
              != (graphContIt->second).end(); nodeIt++)
          {
            Results_OCC res1 = connContRes[*nodeIt][idx1];
            Results_OCC res2 = connContRes[*nodeIt][idx2];
            float minScore = min(res1.score1, res1.score2);
            minScore = min(minScore, res2.score1);
            minScore = min(minScore, res2.score2);
            int mPeaks = min(res1.matchedPeaks, res2.matchedPeaks);
            if (minScore >= minRatio && mPeaks >= minMatchedPeaks)
              continue;
            float fShift = res2.shift1 - res1.shift1;
            if (first || abs(fShift - bestShift) < dist)
            {
              FFShift = fShift;
              FRShift = res2.shift2 - res1.shift1;
              dist = abs(fShift - bestShift);
              first = false;
            }
          }
          float obshiftBT = (prot_match[idx1][2] == prot_match[idx2][2])
              ? FFShift : FRShift;

          bool isReversed = prot_match[idx1][2] != prot_match[idx2][2];

          bool isC = validShiftMod(obshift, idx1, idx2, isReversed);
          bool isCBT = validShiftMod(obshiftBT, idx1, idx2, isReversed);
          if (shiftScore.size() > 0)
          {
            if (isC)
            {
              tp++;
            }
            else
            {
              fp++;
            }
          }
          if (shiftScoreBT.size() > 0 && isCBT && !isC)
          {
            fn++;
          }
          if (shiftScore.size() > 0)
          {
            t++;
          }
          /*
           Spectrum contig1 = contigs[idx1];
           if (reverse1) {contig1.reverse(0.0 - AAJumps::massH2O, 0);}
           Spectrum contig2 = contigs[idx2];
           if (reverse2) {contig2.reverse(0.0 - AAJumps::massH2O, 0);}
           TwoValues<int> mp = getContigOverlapPeaks(contig1, contig2, protFFShift, contigs[idx1].parentMass, 0.5, 0, false);
           if (minMatchedPeaks == 6 && isEqual(minRatio, 0.2, 0.001) && mp[0] < 4 && shiftScore.size() > 0 && isC) {
           if (possiblePairs.count(idx1) > 0) possiblePairs[idx1].insert(idx2);
           else {
           set<int> poss; poss.insert(idx2);
           possiblePairs[idx1] = poss;
           }


           cout << "True positive contig/contig connection w/ < 4 matching peaks between " << idx1 << " (protein " << prot_match[idx1][0] << ") and " << idx2 << " (protein " << prot_match[idx2][0] << ") with:\n";
           for (nodeIt = (graphContIt->second).begin(); nodeIt != (graphContIt->second).end(); nodeIt ++) {
           Results_OCC res1 = connContRes[*nodeIt][idx1];
           Results_OCC res2 = connContRes[*nodeIt][idx2];
           float minScore = min(res1.score1, res1.score2); minScore = min(minScore, res2.score1); minScore = min(minScore, res2.score2);
           int mPeaks = min(res1.matchedPeaks, res2.matchedPeaks);
           if (minScore < minRatio || mPeaks < minMatchedPeaks) continue;
           cout << "mate-pair " << *nodeIt << ", " << mPeaks << " matched peaks, " << minScore << " min ratio, " << res1.shift1 << " Da shift of contig " << idx1 << ", and " << res2.shift1 << " Da shift of contig " << idx2 << "\n";
           }
           cout << "\n";
           cout.flush();

           }
           */
        }
      }
      fprintf(output,
              "%d%s%.3f%s%d%s%d%s%d%s%d%s%.3f%s%.3f\n",
              minMatchedPeaks,
              CSV_SEP,
              minRatio,
              CSV_SEP,
              tp,
              CSV_SEP,
              fn,
              CSV_SEP,
              fp,
              CSV_SEP,
              t,
              CSV_SEP,
              ((float)tp) / ((float)tp + (float)fp),
              CSV_SEP,
              ((float)tp) / ((float)tp + (float)fn));
    }
  }
  fprintf(output,
          "\nMatched b-ions%sMatched y-ions%s(1=true positive 0=false positive)\n",
          CSV_SEP,
          CSV_SEP);
  graphCont.clear();
  connContRes.clear();
  used.clear();
  getAlignmentGraph(&graphCont, &connContRes, precCharge, 6, 0.2, 0.6);
  for (int idx1 = 0; idx1 < contigs.size(); idx1++)
  {
    if (graphCont.count(idx1) == 0)
      continue;
    for (graphContIt = graphCont[idx1].begin(); graphContIt
        != graphCont[idx1].end(); graphContIt++)
    {
      int idx2 = graphContIt->first;
      if (used.count(idx1) > 0 && used[idx1].count(idx2) > 0)
        continue;
      if (used.count(idx2) > 0 && used[idx2].count(idx1) > 0)
        continue;
      if (used.count(idx1) == 0)
      {
        set<unsigned int> myS;
        myS.insert(idx2);
        used[idx1] = myS;
      }
      else
        used[idx1].insert(idx2);

      if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0)
      {
        continue;
      }

      for (nodeIt = (graphContIt->second).begin(); nodeIt
          != (graphContIt->second).end(); nodeIt++)
      {
        Results_OCC res1 = connContRes[*nodeIt][idx1];
        Results_OCC res2 = connContRes[*nodeIt][idx2];
        int connIdx = *nodeIt;
        //if ((idx1 == 207 && (idx2 == 22 || idx2 == 182)) || ((idx1 == 22 || idx1 == 182) && idx2 == 207)) cout << idx1 << "->" << idx2 << ": " << *nodeIt << "\n";
        float minScore = min(res1.score1, res1.score2);
        minScore = min(minScore, res2.score1);
        minScore = min(minScore, res2.score2);
        int mPeaks = min(res1.matchedPeaks, res2.matchedPeaks);
        float fShift = res2.shift1 - res1.shift1, rShift = res2.shift2
            - res1.shift1;
        float obshift = (prot_match[idx1][2] == prot_match[idx2][2]) ? fShift
            : rShift;
        bool isReversed = (prot_match[idx1][2] == prot_match[idx2][2]) ? false
            : true;
        bool reverse1 = prot_match[idx1][2] == 1;
        bool reverse2 = prot_match[idx2][2] == 1;
        float protShift1 = getContigShift(idx1, reverse1);
        float protShift2 = getContigShift(idx2, reverse2);
        float protFFShift = protShift2 - protShift1;
        float protRRShift = (protShift1 + contigs[idx1].parentMass)
            - (protShift2 + contigs[idx2].parentMass);
        protFFShift = (reverse2 == isReversed) ? protFFShift : protRRShift;
        bool correct = isEqual(obshift, protFFShift, 58.0)
            && prot_match[idx1][0] == prot_match[idx2][0];
        int cor = (correct) ? 1 : 0;

        TwoValues<int> symPeaks = getContigOverlapPeaks(connectors[connIdx],
                                                        contigs[idx1],
                                                        res1.shift1,
                                                        res1.shift2,
                                                        0.05,
                                                        0,
                                                        false);
        fprintf(output,
                "%d%s%d%s%d\n",
                symPeaks[0],
                CSV_SEP,
                symPeaks[1],
                CSV_SEP,
                cor);

        symPeaks = getContigOverlapPeaks(connectors[connIdx],
                                         contigs[idx2],
                                         res2.shift1,
                                         res2.shift2,
                                         0.05,
                                         0,
                                         false);
        fprintf(output,
                "%d%s%d%s%d\n",
                symPeaks[0],
                CSV_SEP,
                symPeaks[1],
                CSV_SEP,
                cor);
      }
    }
  }
  /*
   fprintf(output, "\nAll contig connections induced by mate-pair alignments w/ %d matching peaks and %.2f min ratio:\nContig Index 1%sContig Index 2%sConnector Index%sCorrect%sExpected Shift%sObserved Shift%sDifference%sContig 1 Shift%sContig 2 Shift%sMatched Peaks%sOverlapping Ratio%sProtein 1%sProtein 2\n", 10, 0.2, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP, CSV_SEP);
   graphCont.clear(); connContRes.clear(); used.clear();
   getAlignmentGraph(&graphCont, &connContRes, precCharge, 10, 0.2, 0.2);

   for (int idx1 = 0; idx1 < contigs.size(); idx1 ++) {
   if ( graphCont.count(idx1) == 0 ) continue;
   for (graphContIt = graphCont[idx1].begin(); graphContIt != graphCont[idx1].end(); graphContIt ++) {
   int idx2 = graphContIt->first;
   if (used.count(idx1) > 0 && used[idx1].count(idx2) > 0) continue;
   if (used.count(idx2) > 0 && used[idx2].count(idx1) > 0) continue;
   if (used.count(idx1) == 0) {
   set<unsigned int> myS;
   myS.insert(idx2);
   used[idx1] = myS;
   } else used[idx1].insert(idx2);

   if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0) {continue;}

   for (nodeIt = (graphContIt->sec ond).begin(); nodeIt != (graphContIt->second).end(); nodeIt ++)
   {
   Results_OCC res1 = connContRes[*nodeIt][idx1];
   Results_OCC res2 = connContRes[*nodeIt][idx2];
   int connIdx = *nodeIt;
   //if ((idx1 == 207 && (idx2 == 22 || idx2 == 182)) || ((idx1 == 22 || idx1 == 182) && idx2 == 207)) cout << idx1 << "->" << idx2 << ": " << *nodeIt << "\n";
   float minScore = min(res1.score1, res1.score2); minScore = min(minScore, res2.score1); minScore = min(minScore, res2.score2);
   int mPeaks = min(res1.matchedPeaks, res2.matchedPeaks);
   float obshift = (prot_match[idx1][2] == prot_match[idx2][2]) ? res2.shift1 - res1.shift1 : res2.shift2 - res1.shift1;
   bool isReversed = (prot_match[idx1][2] == prot_match[idx2][2]) ? false : true;
   bool reverse1 = prot_match[idx1][2] == 1;
   bool reverse2 = prot_match[idx2][2] == 1;
   float protShift1 = getContigShift(idx1, reverse1);
   float protShift2 = getContigShift(idx2, reverse2);
   float protFFShift = protShift2 - protShift1;
   float protRRShift = (protShift1+contigs[idx1].parentMass) - (protShift2+contigs[idx2].parentMass);
   protFFShift = (reverse2 == isReversed) ? protFFShift : protRRShift;
   bool correct = isEqual(obshift, protFFShift, 58.0) && prot_match[idx1][0] == prot_match[idx2][0];
   float dist = abs(protFFShift - obshift);
   int cor = (correct) ? 1 : 0;
   fprintf(output, "%d%s%d%s%d%s%d%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%d%s%.2f%s%d%s%d\n", idx1, CSV_SEP, idx2, CSV_SEP, connIdx, CSV_SEP, cor, CSV_SEP, protFFShift, CSV_SEP, obshift, CSV_SEP, dist, CSV_SEP, res1.shift1, CSV_SEP, res2.shift1, CSV_SEP, mPeaks, CSV_SEP, minScore, CSV_SEP, prot_match[idx1][0], CSV_SEP, prot_match[idx2][0]);
   }
   }
   }*/
  fclose(output);
}

void CombineContigs::getAlignmentGraph(map<unsigned int, map<unsigned int,
                                           list<unsigned int> > >* graphCont,
                                       map<unsigned int, map<unsigned int,
                                           Results_OCC> >* connContRes,
                                       short precCharge,
                                       unsigned int minNumMatchedPeaks,
                                       float minRatio,
                                       float minContigOverlapIntensity)
{
  list<Results_OCC>::iterator pivot_res;

  //idxNotUsed.insert(4); idxNotUsed.insert(5); idxNotUsed.insert(7); idxNotUsed.insert(8); idxNotUsed.insert(9); idxNotUsed.insert(10); idxNotUsed.insert(11); idxNotUsed.insert(15); idxNotUsed.insert(16); idxNotUsed.insert(17); idxNotUsed.insert(18); idxNotUsed.insert(19); idxNotUsed.insert(20); idxNotUsed.insert(21); idxNotUsed.insert(22);

  //for (set<int>::iterator setit = idxNotUsed.begin(); setit != idxNotUsed.end(); setit++) {cout << *setit << " - ";}
  map<unsigned int, list<unsigned int> > graphConnCont;
  map<unsigned int, list<unsigned int> >::iterator graphConnContIt;
  list<unsigned int>::iterator nodeIt;
  list<unsigned int>::iterator nodeIt2;
  list<unsigned int> newL;
  map<unsigned int, list<unsigned int> > newM;
  map<unsigned int, Results_OCC> newMRes;
  if (graphCont != 0)
    (*graphCont).clear();
  if (connContRes != 0)
    (*connContRes).clear();

  int total = 0;
  int used = 0;

  for (pivot_res = connector_alignments.begin(); pivot_res
      != connector_alignments.end(); pivot_res++)
  {
    Results_OCC result = *pivot_res;
    if (haveRes && (idxNotUsed.count(prot_match[result.spec2][0]) > 0))
      continue;

    total++;
    if (connectors[result.spec1].parentCharge < precCharge)
      continue;
    if (result.matchedPeaks < minNumMatchedPeaks)
      continue;
    if (result.score1 < minRatio || result.score2 < minContigOverlapIntensity)
      continue;
    used++;

    if (graphConnCont.count(result.spec1) == 0)
    {
      newL.clear();
      newL.push_back(result.spec2);
      graphConnCont[result.spec1] = newL;
    }
    else
      graphConnCont[result.spec1].push_back(result.spec2);
    if (connContRes != 0)
    {
      if ((*connContRes).count(result.spec1) == 0)
      {
        newMRes.clear();
        newMRes[result.spec2] = *pivot_res;
        (*connContRes)[result.spec1] = newMRes;
      }
      else
        (*connContRes)[result.spec1][result.spec2] = *pivot_res;
    }
  }

  if (graphCont != 0)
  {
    for (graphConnContIt = graphConnCont.begin(); graphConnContIt
        != graphConnCont.end(); graphConnContIt++)
    {
      for (nodeIt = (graphConnContIt->second).begin(); nodeIt
          != (graphConnContIt->second).end(); nodeIt++)
      {
        for (nodeIt2 = (graphConnContIt->second).begin(); nodeIt2
            != (graphConnContIt->second).end(); nodeIt2++)
        {
          if (*nodeIt2 == *nodeIt)
            continue;
          if ((*graphCont).count(*nodeIt) == 0)
          {
            newM.clear();
            newL.clear();
            newL.push_back(graphConnContIt->first);
            newM[*nodeIt2] = newL;
            (*graphCont)[*nodeIt] = newM;
          }
          else if ((*graphCont)[*nodeIt].count(*nodeIt2) == 0)
          {
            newL.clear();
            newL.push_back(graphConnContIt->first);
            (*graphCont)[*nodeIt][*nodeIt2] = newL;
          }
          else
            (*graphCont)[*nodeIt][*nodeIt2].push_back(graphConnContIt->first);
        }
      }
    }
  }

  cout << "\nUsed " << used << " of " << total
      << " possible connector-contig shifts.\n";
}

void CombineContigs::outputGraph(ostream& output,
                                 map<int, map<int, int> >& _graph,
                                 vector<vector<float> >& _edges,
                                 set<int>* erasedEdges,
                                 vector<set<int> >* _components)
{
  map<int, int>::iterator mapIt;
  set<int>::iterator nodeIt;
  output << "Nodes (edge):\n";
  for (int i = 0; i < contigs.size(); i++)
  {
    if (_graph.count(i) == 0)
      continue;
    output << i << " -> ";
    for (mapIt = _graph[i].begin(); mapIt != _graph[i].end(); mapIt++)
    {
      output << mapIt->first << "(" << mapIt->second << "), ";
    }
    output << "\n";
  }
  output << "\nEdges:\n";
  for (int i = 0; i < _edges.size(); i++)
  {
    if (erasedEdges != 0 && (*erasedEdges).count(i) > 0)
      continue;
    output << i << ": ";
    for (int j = 0; j < _edges[i].size(); j++)
    {
      output << _edges[i][j] << ", ";
    }
    output << "\n";
  }
  if (_components != 0)
  {
    output << "\nComponents:\n";
    for (int i = 0; i < (*_components).size(); i++)
    {
      if ((*_components)[i].size() < 2)
        continue;
      output << i << ": ";
      for (nodeIt = (*_components)[i].begin(); nodeIt
          != (*_components)[i].end(); nodeIt++)
      {
        output << *nodeIt;
        if (*nodeIt == i)
          cout << "(0), ";
        else
          cout << "(" << _edges[_graph[i][*nodeIt]][2] << ","
              << _edges[_graph[i][*nodeIt]][5] << "), ";
      }
      output << "\n";
    }
  }
}

float CombineContigs::getComponentOverlapArea(int idx1,
                                              int idx2,
                                              float shift_use,
                                              bool rev,
                                              map<int, map<int, int> >& _graph,
                                              vector<vector<float> >& _edges,
                                              vector<set<int> >& _components)
{

  float idx1_left_bound = 0;
  float idx1_right_bound = oriented[idx1].parentMass;

  float idx2_left_bound = shift_use;
  float idx2_right_bound = shift_use + oriented[idx2].parentMass;

  for (set<int>::iterator nodeIt = _components[idx1].begin(); nodeIt
      != _components[idx1].end(); nodeIt++)
  {
    if (*nodeIt == idx1)
      continue;
    float ffshift = _edges[_graph[idx1][*nodeIt]][2];

    idx1_left_bound = min(idx1_left_bound, ffshift);
    idx1_right_bound = max(idx1_right_bound, ffshift
        + oriented[*nodeIt].parentMass);
  }
  for (set<int>::iterator nodeIt = _components[idx2].begin(); nodeIt
      != _components[idx2].end(); nodeIt++)
  {
    if (*nodeIt == idx2)
      continue;
    float ffshift = (rev) ? _edges[_graph[idx2][*nodeIt]][5]
        : _edges[_graph[idx2][*nodeIt]][2];

    idx2_left_bound = min(idx2_left_bound, ffshift + shift_use);
    idx2_right_bound = max(idx2_right_bound, ffshift + shift_use
        + oriented[*nodeIt].parentMass);
  }

  float left_bound = max(idx1_left_bound, idx2_left_bound);
  float right_bound = min(idx1_right_bound, idx2_right_bound);

  return right_bound - left_bound;
}

void CombineContigs::outputComponents(FILE* output,
                                      map<int, map<int, int> >& _graph,
                                      vector<vector<float> >& _edges,
                                      vector<set<int> >& _components,
                                      SpecSet& orientedContigs,
                                      vector<bool>& _reversed,
                                      bool outputStats,
                                      bool printComps,
                                      set<int>* countNodes,
                                      bool justContigs)
{

  bool checkPossible = false;
  if (outputStats && haveRes && numCorrectPairs == 0)
  {
    checkPossible = true;
    int possible = 0;
    map<int, set<int> > seen;
    set<int> sit;
    for (int i = 0; i < contigs.size(); i++)
    {
      for (int j = 0; j < contigs.size(); j++)
      {
        int idx1 = i, idx2 = j;
        if (idx1 == idx2)
          continue;
        if (idxNotUsed.count(prot_match[idx1][0]) > 0
            || idxNotUsed.count(prot_match[idx2][0]) > 0)
          continue;
        if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0
            || prot_match[idx1][0] != prot_match[idx2][0])
          continue;
        if (seen.count(i) > 0 && seen[i].count(j) > 0)
          continue;
        if (seen.count(j) > 0 && seen[j].count(i) > 0)
          continue;
        if (seen.count(i) > 0)
          seen[i].insert(j);
        else
        {
          sit.clear();
          sit.insert(j);
          seen[i] = sit;
        }
        bool reverse1 = prot_match[idx1][2] == 1;
        bool reverse2 = prot_match[idx2][2] == 1;
        float protShift1 = getContigShift(idx1, reverse1);
        float protShift2 = getContigShift(idx2, reverse2);
        float protFFShift = protShift2 - protShift1;
        //float protRRShift = (protShift1+contigs[idx1].parentMass) - (protShift2+contigs[idx2].parentMass);
        //protFFShift = (reverse1) ? protRRShift : protFFShift;
        //if (protFFShift < 0 - contigs[idx2].parentMass || protFFShift > contigs[idx1].parentMass) continue;
        Spectrum contig1 = contigs[idx1];
        if (reverse1)
        {
          contig1.reverse(0.0 - AAJumps::massH2O, 0);
        }
        Spectrum contig2 = contigs[idx2];
        if (reverse2)
        {
          contig2.reverse(0.0 - AAJumps::massH2O, 0);
        }
        TwoValues<int> mp = getContigOverlapPeaks(contig1,
                                                  contig2,
                                                  protFFShift,
                                                  contigs[idx1].parentMass,
                                                  0.5,
                                                  0,
                                                  false);

        if (mp[0] < 6)
          continue;

        if (correctPairs.count(i) > 0)
          correctPairs[i].insert(j);
        else
        {
          sit.clear();
          sit.insert(j);
          correctPairs[i] = sit;
        }

        possible++;
      }
    }
    numCorrectPairs = possible;
  }
  bool addTo = seenBef.size() == 0;
  float totalLen = 0, numLen = 0, maxLen = 0, totalCover = 0, maxCover = 0;

  set<int>::iterator nodeIt;
  set<int>::iterator nodeIt2;
  if (haveRes)
  {
    if (output != 0)
      fprintf(output,
              "\nComponent Idx%sProtein Idx%sContig Root Idx%sContig Idx2%sExperimental FF Shift%sExperimental RR Shift%sExpected FF Shift%sExpected RR Shift%sHas Good Edge%sComponent Length (Da)%sComponent Length (AA)\n",
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP);
  }
  else
  {
    if (output != 0)
      fprintf(output,
              "\nComponent Idx%sContig Root Idx%sContig Idx2%sExperimental FF Shift%sExperimental RR Shift%sComponent Length (Da)\n",
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP,
              CSV_SEP);
  }
  map<int, int> componentRef;
  for (int i = 0; i < contigs.size(); i++)
  {
    componentRef[i] = -1;
  }
  map<int, set<int> > seenNodes;
  set<int> snode;
  float goodEdges = 0, edgesAdded = 0;
  int numCompsWSize = 0;
  int compIdx = 0;
  int minCompSize = 5;
  int good_aligned = 0;
  if (haveRes && printComps)
  {
    printf("\nComponent Index - Length (AA) - Protein Index - Protein Header\n");
  }
  for (int i = 0; i < _components.size(); i++)
  {
    //if (_components[i].size() == 0) continue;
    if (countNodes != 0 && (*countNodes).count(i) == 0)
      continue;

    if (_components[i].size() <= 1)
      continue;

    if (_components[i].size() >= minCompSize)
    {
      numCompsWSize++;
    }

    float len = 0, AAlen = 0, coverAA = 0;
    vector<int> protids(fasta.size());
    for (int k = 0; k < protids.size(); k++)
      protids[k] = 0;

    if (haveRes && prot_match[i][0] >= 0 && idxNotUsed.count(prot_match[i][0])
        == 0)
    {
      coverAA = abs(overlaps[i][overlaps[i].size() - 1][1] - overlaps[i][0][1])
          + 1;
      AAlen = max(AAlen, coverAA);
      len = contigs[i].parentMass;
    }
    if (!haveRes)
    {
      len = contigs[i].parentMass;
    }
    for (nodeIt = _components[i].begin(); nodeIt != _components[i].end(); nodeIt++)
    {
      componentRef[*nodeIt] = compIdx;
      float FFShift = (i == *nodeIt) ? 0.0 : _edges[_graph[i][*nodeIt]][2];
      float RRShift = (i == *nodeIt) ? 0.0 : _edges[_graph[i][*nodeIt]][5];
      if (haveRes && prot_match[*nodeIt][0] >= 0)
      {
        protids[prot_match[*nodeIt][0]] += 1;
      }
      if (haveRes)
      {
        for (nodeIt2 = _components[i].begin(); nodeIt2 != _components[i].end(); nodeIt2++)
        {
          int idx1 = *nodeIt, idx2 = *nodeIt2;
          if (idx1 == idx2)
            continue;
          if (idxNotUsed.count(prot_match[idx1][0]) > 0
              || idxNotUsed.count(prot_match[idx2][0]) > 0)
          {
            continue;
          }
          if ((seenNodes.count(idx1) > 0 && seenNodes[idx1].count(idx2) > 0)
              || (seenNodes.count(idx2) > 0 && seenNodes[idx2].count(idx1) > 0))
          {
            continue;
          }
          if (seenNodes.count(idx1) == 0)
          {
            snode.clear();
            snode.insert(idx2);
            seenNodes[idx1] = snode;
          }
          else
          {
            seenNodes[idx1].insert(idx2);
          }

          float shift = abs(_edges[_graph[i][idx1]][2]
              - _edges[_graph[i][idx2]][2]);
          if (_edges[_graph[i][idx1]][2] + contigs[idx1].parentMass
              > _edges[_graph[i][idx2]][2] + contigs[idx2].parentMass)
          {
            shift += contigs[idx1].parentMass;
          }
          else
          {
            shift += contigs[idx2].parentMass;
          }
          len = max(len, shift);

          if (prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0)
            continue;

          float FFShift2 = (i == idx2) ? 0.0 : _edges[_graph[i][idx2]][2];
          float RRShift2 = (i == idx2) ? 0.0 : _edges[_graph[i][idx2]][5];

          coverAA = overlaps[idx2][overlaps[idx2].size() - 1][1]
              - overlaps[idx2][0][1] + 1.0;
          coverAA = max(coverAA, overlaps[idx1][overlaps[idx1].size() - 1][1]
              - overlaps[idx1][0][1] + (float)1.0);

          if (prot_match[idx1][0] == prot_match[idx2][0]
              && (validShiftMod(FFShift2 - FFShift, idx1, idx2, _reversed[idx2])
                  || validShiftMod(RRShift2 - RRShift,
                                   idx1,
                                   idx2,
                                   !_reversed[idx2])))
          {
            goodEdges += 1.0;

            float maxPos = max(overlaps[idx2][overlaps[idx2].size() - 1][1],
                               overlaps[idx1][overlaps[idx1].size() - 1][1]);
            float minPos = min(overlaps[idx2][0][1], overlaps[idx1][0][1]);
            coverAA = max(coverAA, maxPos - minPos + (float)1.0);
            TwoValues<int> mp = getContigOverlapPeaks(orientedContigs[idx1],
                                                      orientedContigs[idx2],
                                                      FFShift2 - FFShift,
                                                      contigs[idx1].parentMass,
                                                      1.0,
                                                      0,
                                                      false);
            if (mp[0] >= 6 && ((correctPairs.count(idx1) > 0
                && correctPairs[idx1].count(idx2) > 0)
                || (correctPairs.count(idx2) > 0
                    && correctPairs[idx2].count(idx1) > 0)))
            {
              if (addTo)
              {
                if (seenBef.count(idx1) > 0)
                  seenBef[idx1].insert(idx2);
                else
                {
                  snode.clear();
                  snode.insert(idx2);
                  seenBef[idx1] = snode;
                }
              }
              else
              {
                if ((seenBef.count(idx1) == 0 || seenBef[idx1].count(idx2) == 0)
                    && (seenBef.count(idx2) == 0 || seenBef[idx2].count(idx1)
                        == 0))
                {
                  cout << "\nFOUND IT: " << idx1 << " " << idx2 << "\n";
                }
              }

              good_aligned++;
            }
          }

          AAlen = max(AAlen, coverAA);
          edgesAdded += 1.0;
        }
      }
      else
      {
        for (nodeIt2 = _components[i].begin(); nodeIt2 != _components[i].end(); nodeIt2++)
        {
          int idx1 = *nodeIt, idx2 = *nodeIt2;
          if (idx1 == idx2)
            continue;
          float shift = abs(_edges[_graph[i][idx1]][2]
              - _edges[_graph[i][idx2]][2]);
          if (_edges[_graph[i][idx1]][2] + contigs[idx1].parentMass
              > _edges[_graph[i][idx2]][2] + contigs[idx2].parentMass)
          {
            shift += contigs[idx1].parentMass;
          }
          else
          {
            shift += contigs[idx2].parentMass;
          }
          len = max(len, shift);
        }
      }
    }
    if (_components[i].size() >= minCompSize)
    {
      totalLen += len;
      numLen += 1.0;
      maxLen = max(maxLen, len);
    }

    if (haveRes)
    {
      int bestprot = -1, maxCount = 0;
      for (int k = 0; k < protids.size(); k++)
      {
        if (protids[k] > maxCount)
        {
          bestprot = k;
          maxCount = protids[k];
        }
      }
      if (_components[i].size() >= minCompSize)
      {
        maxCover = max(maxCover, AAlen);
        totalCover += AAlen;
      }

      const char* protident = "none";
      if (bestprot >= 0)
        protident = fasta.getID(bestprot);
      if (printComps)
        printf("%d - %.0f - %d - %s\n", compIdx, AAlen, bestprot, protident);
    }
    for (nodeIt = _components[i].begin(); nodeIt != _components[i].end(); nodeIt++)
    {
      float FFShift = (i == *nodeIt) ? 0.0 : _edges[_graph[i][*nodeIt]][2];
      float RRShift = (i == *nodeIt) ? 0.0 : _edges[_graph[i][*nodeIt]][5];
      if (haveRes)
      {
        string protFFStr = "";
        string protRRStr = "";
        const char* eqStr = "";
        if (prot_match[*nodeIt][0] >= 0 && prot_match[i][0] >= 0)
        {
          bool reverse1 = prot_match[i][2] == 1;
          bool reverse2 = prot_match[*nodeIt][2] == 1;

          bool sameProt = prot_match[i][0] == prot_match[*nodeIt][0];
          float protShift1 = getContigShift(i, reverse1);
          float protShift2 = getContigShift(*nodeIt, reverse2);
          float protFFShift = protShift2 - protShift1;
          float protRRShift = (protShift1 + contigs[i].parentMass)
              - (protShift2 + contigs[*nodeIt].parentMass);
          float temF = protFFShift, temR = protRRShift;
          protFFShift = (reverse1 == _reversed[i]) ? temF : temR;
          protRRShift = (reverse1 == _reversed[*nodeIt]) ? temR : temF;
          if (sameProt)
          {
            protFFStr = parseFloat(protFFShift, 3);
            protRRStr = parseFloat(protRRShift, 3);
            eqStr = (isEqual(protFFShift, FFShift, 58.0)) ? (char*)"1"
                : (char*)"0";
          }
          else
          {
            protFFStr = "nan";
            protRRStr = "nan";
            eqStr = "0";
          }
        }
        if (output != 0)
          fprintf(output,
                  "%d%s%d%s%d%s%d%s%.1f%s%.1f%s%s%s%s%s%s%s%.1f%s%.0f\n",
                  compIdx,
                  CSV_SEP,
                  prot_match[*nodeIt][0],
                  CSV_SEP,
                  i,
                  CSV_SEP,
                  *nodeIt,
                  CSV_SEP,
                  FFShift,
                  CSV_SEP,
                  RRShift,
                  CSV_SEP,
                  protFFStr.c_str(),
                  CSV_SEP,
                  protRRStr.c_str(),
                  CSV_SEP,
                  eqStr,
                  CSV_SEP,
                  len,
                  CSV_SEP,
                  AAlen);
      }
      else
      {
        if (output != 0)
          fprintf(output,
                  "%d%s%d%s%d%s%.1f%s%.1f%s%.1f\n",
                  compIdx,
                  CSV_SEP,
                  i,
                  CSV_SEP,
                  *nodeIt,
                  CSV_SEP,
                  FFShift,
                  CSV_SEP,
                  RRShift,
                  CSV_SEP,
                  len);
      }
    }
    compIdx++;
  }

  if (addTo)
  {
    for (map<int, set<int> >::iterator possIt = correctPairs.begin(); possIt
        != correctPairs.end(); possIt++)
    {
      int idx1 = possIt->first;
      set<int> idx2s = possIt->second;
      for (set<int>::iterator secIt = idx2s.begin(); secIt != idx2s.end(); secIt++)
      {
        int idx2 = *secIt;
        if (componentRef[idx1] != componentRef[idx2] && componentRef[idx1] >= 0
            && componentRef[idx2] >= 0)
        {// && (seenNodes.count(idx1) == 0 || seenNodes[idx1].count(idx2) == 0) && (seenNodes.count(idx2) == 0 || seenNodes[idx2].count(idx1) == 0)) {
          //cout << idx1 << "-" << componentRef[idx1] << ", " << idx2 << "-" << componentRef[idx2] << ": Missed true positive contig/contig connection w/ >= 6 matching peaks between " << idx1 << " (protein " << prot_match[idx1][0] << ") and " << idx2 << " (protein " << prot_match[idx2][0] << ")\n";
        }
      }
    }
  }
  /*
   if (checkPossible) {
   for (map<int, set<int> >::iterator possIt = possiblePairs.begin(); possIt != possiblePairs.end(); possIt ++) {
   int idx1 = possIt->first;
   set<int> idx2s = possIt->second;
   for (set<int>::iterator secIt = idx2s.begin(); secIt != idx2s.end(); secIt ++) {
   int idx2 = *secIt;
   if (componentRef[idx1] != componentRef[idx2] && componentRef[idx1] >= 0 && componentRef[idx2] >= 0 && (seenNodes.count(idx1) == 0 || seenNodes[idx1].count(idx2) == 0) && (seenNodes.count(idx2) == 0 || seenNodes[idx2].count(idx1) == 0)) {
   cout << idx1 << "-" << componentRef[idx1] << ", " << idx2 << "-" << componentRef[idx2] << ": True positive contig/contig connection w/ < 4 matching peaks between " << idx1 << " (protein " << prot_match[idx1][0] << ") and " << idx2 << " (protein " << prot_match[idx2][0] << ")\n";
   }
   }
   }
   }
   */
  if (outputStats)
  {
    printf("\nFound %d components (%d with %d contigs or more)\n",
           compIdx,
           numCompsWSize,
           minCompSize);
    if (haveRes)
    {
      printf("Combined %s %.2f good edges (%.0f/%.0f)\n", "%", (100.0
          * goodEdges) / edgesAdded, goodEdges, edgesAdded);
    }
    printf("Maximum component length (Da) : %.1f\n", maxLen);
    printf("Average component length (Da) : %.1f\n", totalLen / numLen);
    if (haveRes)
    {
      printf("Maximum component length (AA) : %.0f\n", maxCover);
      printf("Average component length (AA) : %.1f\n", totalCover / numLen);
      printf("Found %s %.2f (%d/%d) of possible contig pairs with 6 matching peaks\n",
             "%",
             (100.0 * (float)good_aligned) / ((float)numCorrectPairs),
             good_aligned,
             numCorrectPairs);
    }
  }
}

void CombineContigs::recomputeEdges(int idx,
                                    map<int, map<int, int> >& _graph,
                                    vector<vector<float> >& _edges,
                                    vector<set<int> >& _components,
                                    set<int>& erasedEdges,
                                    vector<int>& edgeRef,
                                    vector<list<TwoValues<int> > >& _root_alignments,
                                    list<TwoValues<float> >& scoredEdges)
{
  map<int, float> newEdgeScores;
  map<int, int>::iterator mapIt;
  list<TwoValues<float> >::iterator scoredEdgeIt;
  float newEdgeScore;
  set<int> seen;
  int oldSize = edgeRef.size();

  for (mapIt = _graph[idx].begin(); mapIt != _graph[idx].end(); mapIt++)
  {
    int node2 = mapIt->first;
    int edgeToNode = mapIt->second;
    seen.insert(node2);
    if (_components[idx].count(node2) > 0 || node2 == idx)
      continue;
    int edgeFromNode = _graph[node2][idx];
    float node2FShift = _edges[edgeToNode][2];
    float node2RShift = _edges[edgeToNode][3];
    TwoValues<float> res = getConsensusOverlap(idx,
                                               node2,
                                               node2FShift,
                                               node2RShift,
                                               _graph,
                                               _edges,
                                               _components,
                                               _root_alignments,
                                               false);
    newEdgeScore = max(res[0], res[1]);
    newEdgeScores[edgeToNode] = newEdgeScore;
    newEdgeScores[edgeFromNode] = newEdgeScore;
  }

  for (scoredEdgeIt = scoredEdges.begin(); scoredEdgeIt != scoredEdges.end(); scoredEdgeIt++)
  {
    int nextEdgeIdx = floatToInt((*scoredEdgeIt)[0]);
    if (nextEdgeIdx == -1)
      continue;
    while (edgeRef[nextEdgeIdx] != -1)
      nextEdgeIdx = edgeRef[nextEdgeIdx];
    if (erasedEdges.count(nextEdgeIdx) > 0)
      continue;
    if (newEdgeScores.count(nextEdgeIdx) > 0)
      (*scoredEdgeIt)[1] = newEdgeScores[nextEdgeIdx];
  }

  scoredEdges.sort(SortEdges());
}

TwoValues<float> CombineContigs::getOverlapScore(int idx1,
                                                 int idx2,
                                                 float fShift,
                                                 float rShift,
                                                 map<int, map<int, int> >& _graph,
                                                 vector<vector<float> >& _edges,
                                                 vector<set<int> >& _components)
{
  float totalScoreF = 0, numScoreF = 0, totalScoreR = 0, numScoreR = 0;
  set<int>::iterator nodeIt, nodeIt2;
  TwoValues<float> res;
  TwoValues<int> mp;
  for (nodeIt = _components[idx1].begin(); nodeIt != _components[idx1].end(); nodeIt++)
  {
    int node1 = *nodeIt;
    float shift1 = 0;
    if (node1 != idx1)
      shift1 = _edges[_graph[idx1][node1]][2];

    for (nodeIt2 = _components[idx2].begin(); nodeIt2
        != _components[idx2].end(); nodeIt2++)
    {
      int node2 = *nodeIt2;
      if (node1 == node2)
      {
        continue;
      }
      float shift2 = 0, shiftr2 = 0;
      if (node2 != idx2)
      {
        shift2 = _edges[_graph[idx2][node2]][2];
        shiftr2 = _edges[_graph[idx2][node2]][5];
      }

      float locFShift = fShift - shift1 + shift2;
      float locRShift = rShift - shift1 + shiftr2;

      bool debug = false;//(idx1 == 35 and idx2 == 560);
      res = getContigOverlapScores(oriented[node1],
                                   oriented[node2],
                                   locFShift,
                                   locRShift,
                                   parentMassTol,
                                   minContigOverlap,
                                   debug);
      mp = getContigOverlapPeaks(oriented[node1],
                                 oriented[node2],
                                 locFShift,
                                 locRShift,
                                 parentMassTol,
                                 minContigOverlap,
                                 debug);
      //if (idx1 == 35 and idx      2 == 560)
      {
        cout << res[0] << " AND " << res[1] << ", " << mp[0] << " AND "
            << mp[1] << "\n";
      }
      /*
       if (idx1 == 51 and idx2 == 409) {
       float protShift1 = getContigShift(node1, prot_match[node1][2] == 1);
       float protShift2 = getContigShift(node2, prot_match[node2][2] == 1);
       float protFFShift = protShift2 - protShift1;
       float protRRShift = (protShift1+contigs[node1].parentMass) - (protShift2+contigs[node2].parentMass);
       cout << "overlapping " << node1 << " and " << node2 << ": " << res[0] << " " << res[1] << ":\n";
       cout << "shifts are " << locFShift << " and " << locRShift << ": one should match " << protFFShift << " or " << protRRShift << "\n\n";
       }
       */
      if (res[0] > -0.01)
      {
        numScoreF += 1.0;
        totalScoreF += res[0];
      }
      if (res[1] > -0.01)
      {
        numScoreR += 1.0;
        totalScoreR += res[1];
      }
    }
  }
  if (numScoreF < 0.01)
  {
    totalScoreF = 1.0;
    numScoreF = -1.0;
  }
  if (numScoreR < 0.01)
  {
    totalScoreR = 1.0;
    numScoreR = -1.0;
  }

  return TwoValues<float> (totalScoreF / numScoreF, totalScoreR / numScoreR);
}

TwoValues<float> CombineContigs::getConsensusShifts(int idx1,
                                                    int idx2,
                                                    float fShift,
                                                    float rShift,
                                                    map<int, map<int, int> >& _graph,
                                                    vector<vector<float> >& _edges,
                                                    vector<set<int> >& _components,
                                                    bool debug)
{

  float totalScoreF, totalScoreR;
  float leftEdgeF = 0, rightEdgeF = 0, rightEdgeR = 0;
  for (set<int>::iterator nodeIt = _components[idx1].begin(); nodeIt
      != _components[idx1].end(); nodeIt++)
  {
    if (*nodeIt == idx1)
      continue;
    leftEdgeF = min(leftEdgeF, _edges[_graph[idx1][*nodeIt]][2]);
  }
  for (set<int>::iterator nodeIt = _components[idx2].begin(); nodeIt
      != _components[idx2].end(); nodeIt++)
  {
    if (*nodeIt == idx2)
      continue;
    rightEdgeF = min(rightEdgeF, _edges[_graph[idx2][*nodeIt]][2]);
    rightEdgeR = min(rightEdgeR, _edges[_graph[idx2][*nodeIt]][5]);
  }
  if (debug)
  {
    cout << "Root shift between " << idx1 << " and " << idx2 << " = " << fShift
        << ", " << rShift << "\n";
    cout << "consensusFShift = " << fShift << " - " << leftEdgeF << " + "
        << rightEdgeF << " = " << fShift - leftEdgeF + rightEdgeF << "\n";
    cout << "consensusRShift = " << rShift << " - " << leftEdgeF << " + "
        << rightEdgeR << " = " << rShift - leftEdgeF + rightEdgeR << "\n";
  }
  float consensusFShift = fShift - leftEdgeF + rightEdgeF;
  float consensusRShift = rShift - leftEdgeF + rightEdgeR;
  return TwoValues<float> (consensusFShift, consensusRShift);
}

TwoValues<float> CombineContigs::getConsensusOverlap(int idx1,
                                                     int idx2,
                                                     float fShift,
                                                     float rShift,
                                                     map<int, map<int, int> >& _graph,
                                                     vector<vector<float> >& _edges,
                                                     vector<set<int> >& _components,
                                                     vector<
                                                         list<TwoValues<int> > >& _root_alignments,
                                                     bool debug)
{

  TwoValues<float> consensusShifts = getConsensusShifts(idx1,
                                                        idx2,
                                                        fShift,
                                                        rShift,
                                                        _graph,
                                                        _edges,
                                                        _components,
                                                        debug);

  Spectrum spec1 = consensus[idx1], spec2 = consensus[idx2];
  TwoValues<float> extraShift1 = extraShifts[idx1], extraShift2 =
      extraShifts[idx2];

  float consensusFShift = consensusShifts[0] - extraShift1[0] + extraShift2[0];
  float consensusRShift = consensusShifts[1] - extraShift1[0] + extraShift2[1];
  /*
   if (isEqual(spec1.parentMass - AAJumps::massHion - spec1.peakList.front()[0], AAJumps::massH2O, 0.01)) {
   for (int i = 0; i < spec1.size(); i++) {
   spec1[i][0] -= AAJumps::massH2O;
   }
   }
   */
  if (debug)
  {
    cout << "\nconsensusFShift = " << consensusFShift << " = "
        << consensusShifts[0] << " - " << extraShift1[0] << " + "
        << extraShift2[0] << "\n";
    cout << "consensusRShift = " << consensusRShift << " = "
        << consensusShifts[1] << " - " << extraShift1[0] << " + "
        << extraShift2[1] << "\n";
    cout << "Consensus " << idx1 << ":\n";
    spec1.output(cout);
    cout << "\nConsensus " << idx2 << ":\n";
    spec2.output(cout);
    cout << "\nConsensus " << idx2 << " (reversed):\n";
    Spectrum tempContig = spec2;
    tempContig.reverse(0.0 - AAJumps::massH2O, 0);
    tempContig.output(cout);
  }
  TwoValues<float> res = getContigOverlapScores(spec1,
                                                spec2,
                                                consensusFShift,
                                                consensusRShift,
                                                parentMassTol,
                                                minContigOverlap,
                                                debug);
  TwoValues<int> mp = getContigOverlapPeaks(spec1,
                                            spec2,
                                            consensusFShift,
                                            consensusRShift,
                                            parentMassTol,
                                            minContigOverlap,
                                            false);

  if (debug)
  {
    cout << "Forward score = " << res[0] << "(" << mp[0]
        << " matching peaks)\n";
    cout << "Reverse score = " << res[1] << "(" << mp[1]
        << " matching peaks)\n";
  }
  //if (res[0] > -0.01) res[0] = 0;
  //if (res[1] > -0.01) res[1] = 0;

  if ((idx1 == 67 && idx2 == 712) || (idx1 == 712 && idx2 == 67))
  {
    cout << "new score between " << idx1 << " and " << idx2 << ": " << res[0]
        * ((float)mp[0]) << ", " << res[1] * ((float)mp[1]) << endl;
  }

  return TwoValues<float> (res[0] * ((float)mp[0]), res[1] * ((float)mp[1]));
}

bool CombineContigs::tryMergeContigs(int idx1,
                                     int idx2,
                                     int edgeIdx,
                                     map<int, map<int, int> >& _graph,
                                     vector<vector<float> >& _edges,
                                     vector<set<int> >& _components,
                                     vector<list<TwoValues<int> > >& _root_alignments,
                                     bool rev,
                                     bool debug)
{
  float conFFShift = (rev) ? _edges[edgeIdx][3] : _edges[edgeIdx][2];

  int most_matched_peaks = 0;
  TwoValues<int> best_node_pair;
  best_node_pair[0] = idx1;
  best_node_pair[1] = idx2;
  set<int>::iterator nodeIt, nodeIt2;
  for (nodeIt = _components[idx1].begin(); nodeIt != _components[idx1].end(); nodeIt++)
  {
    if (debug)
    {
      cout << "idx1 = " << idx1 << ", nodeIt = " << *nodeIt << "\n";
    }
    float shift1 = (*nodeIt == idx1) ? 0 : _edges[_graph[idx1][*nodeIt]][2];
    for (nodeIt2 = _components[idx2].begin(); nodeIt2
        != _components[idx2].end(); nodeIt2++)
    {
      if (*nodeIt2 == *nodeIt)
        continue;
      float shift_use = (rev) ? _edges[_graph[idx2][*nodeIt2]][5]
          : _edges[_graph[idx2][*nodeIt2]][2];
      float shift2 = (*nodeIt2 == idx2) ? 0 : shift_use;
      float locFShift = conFFShift - shift1 + shift2;
      Spectrum spec2 = oriented[*nodeIt2];
      if (rev)
        spec2.reverse(0.0 - AAJumps::massH2O);
      TwoValues<int> mp = getContigOverlapPeaks(oriented[*nodeIt],
                                                spec2,
                                                locFShift,
                                                0,
                                                parentMassTol,
                                                minContigOverlap,
                                                debug);
      if (mp[0] > most_matched_peaks)
      {
        best_node_pair[0] = *nodeIt;
        best_node_pair[1] = *nodeIt2;
        most_matched_peaks = mp[0];
      }
    }
  }
  if (most_matched_peaks == 0)
  {
    cout << "\nCOULD NOT FIND OVERLAPPING CONTIGS BETWEEN COMONENTS " << idx1
        << " AND " << idx2 << "\n";
  }
  if (debug)
  {
    cout << "\nbest node pair: " << best_node_pair[0] << " and "
        << best_node_pair[1] << "(" << most_matched_peaks
        << " matching peaks)\n";
  }
  list<TwoValues<int> > idx2_aligns(_root_alignments[idx1]);
  for (list<TwoValues<int> >::iterator alIt = _root_alignments[idx2].begin(); alIt
      != _root_alignments[idx2].end(); alIt++)
  {
    idx2_aligns.push_back(*alIt);
  }
  idx2_aligns.push_back(best_node_pair);

  if (debug)
  {
    cout << "_root_alignments[idx1].size() = " << _root_alignments[idx1].size()
        << "\n";
    cout << "_components[idx1].size() = " << _components[idx1].size() << "\n";
    cout << "_root_alignments[idx2].size() = " << _root_alignments[idx2].size()
        << "\n";
    cout << "_components[idx2].size() = " << _components[idx2].size() << "\n";
    cout << "Merged " << idx2_aligns.size() << " total alignments for "
        << _components[idx1].size() + _components[idx2].size() << " spectra\n";
  }
  //list<TwoValues<int> > idx2_aligns = root_alignments[idx2];
  //    idx2_aligns.push_back(best_node_pair);
  //    root_alignments[idx1].merge(idx2_aligns);
  //    components[idx1].insert(idx2);

  Results_PA next;
  vector<Results_PA> outAlign;

  SpecSet inspec;
  abinfo_t inabinfo;

  float minLeftEdgeF = 0;
  float maxRightEdge = oriented[idx1].peakList.back()[0];
  map<int, TwoValues<float> > contigIndent;
  masab_parameters m_params;

  m_params.INPUT_SPECS = oriented;

  //debug = idx == 0;

  if (debug)
    cout << "\nInput spectra to masab for " << idx1 << ":\n";
  for (nodeIt = _components[idx1].begin(); nodeIt != _components[idx1].end(); nodeIt++)
  {
    //node_to_index[*nodeIt] = outs.specs.size();
    //index_to_node[outs.specs.size()] = *nodeIt;
    //outs.specs.push_back(oriented[*nodeIt]);
    float FFShift = (*nodeIt == idx1) ? 0 : _edges[_graph[idx1][*nodeIt]][2];
    minLeftEdgeF = min(minLeftEdgeF, FFShift);
    float righEdge = FFShift + oriented[*nodeIt].peakList.back()[0];
    maxRightEdge = max(maxRightEdge, righEdge);
    if (debug)
    {
      cout << "\nContig " << *nodeIt << "(shift=" << FFShift << "):\n";
      oriented[*nodeIt].output(cout);
    }
  }
  for (nodeIt = _components[idx2].begin(); nodeIt != _components[idx2].end(); nodeIt++)
  {
    //node_to_index[*nodeIt] = outs.specs.size();
    //index_to_node[outs.specs.size()] = *nodeIt;
    //outs.specs.push_back(oriented[*nodeIt]);
    if (rev)
    {
      m_params.INPUT_SPECS[*nodeIt].reverse(0.0 - AAJumps::massH2O);
    }
    float shift_use = (rev) ? _edges[_graph[idx2][*nodeIt]][5]
        : _edges[_graph[idx2][*nodeIt]][2];
    float FFShift = (*nodeIt == idx2) ? conFFShift : conFFShift + shift_use;
    minLeftEdgeF = min(minLeftEdgeF, FFShift);
    float righEdge = FFShift + oriented[*nodeIt].peakList.back()[0];
    maxRightEdge = max(maxRightEdge, righEdge);
    if (debug)
    {
      cout << "\nContig " << *nodeIt << "(shift=" << FFShift << "):\n";
      oriented[*nodeIt].output(cout);
    }
  }

  for (set<int>::iterator nodeIt = _components[idx1].begin(); nodeIt
      != _components[idx1].end(); nodeIt++)
  {
    //int locIdx = node_to_index[*nodeIt];
    float FFShift = (*nodeIt == idx1) ? 0 : _edges[_graph[idx1][*nodeIt]][2];
    contigIndent[*nodeIt][0] = FFShift - minLeftEdgeF;
  }

  for (set<int>::iterator nodeIt = _components[idx2].begin(); nodeIt
      != _components[idx2].end(); nodeIt++)
  {
    //int locIdx = node_to_index[*nodeIt];
    float shift_use = (rev) ? _edges[_graph[idx2][*nodeIt]][5]
        : _edges[_graph[idx2][*nodeIt]][2];
    float FFShift = (*nodeIt == idx2) ? conFFShift : conFFShift + shift_use;
    contigIndent[*nodeIt][0] = FFShift - minLeftEdgeF;
  }

  maxRightEdge -= minLeftEdgeF;
  if (debug)
    cout << "Found " << idx2_aligns.size() << " alignments\n";

  for (list<TwoValues<int> >::iterator node_pair_it =

  idx2_aligns.begin(); node_pair_it != idx2_aligns.end(); node_pair_it++)
  {
    int locIdx1 = (*node_pair_it)[0], locIdx2 = (*node_pair_it)[1];
    //int locIdx1 = node_to_index[idx1], locIdx2 = node_to_index[idx2];

    float FFShift;
    if (_components[idx1].count(locIdx1) > 0)
    {
      FFShift = (locIdx1 == idx1) ? 0 : _edges[_graph[idx1][locIdx1]][2];
    }
    else
    {
      float shift_use = (rev) ? _edges[_graph[idx2][locIdx1]][5]
          : _edges[_graph[idx2][locIdx1]][2];
      FFShift = (locIdx1 == idx2) ? conFFShift : conFFShift + shift_use;
    }
    float FFShift2;
    if (_components[idx1].count(locIdx2) > 0)
    {
      FFShift2 = (locIdx2 == idx1) ? 0 : _edges[_graph[idx1][locIdx2]][2];
    }
    else
    {
      float shift_use = (rev) ? _edges[_graph[idx2][locIdx2]][5]
          : _edges[_graph[idx2][locIdx2]][2];
      FFShift2 = (locIdx2 == idx2) ? conFFShift : conFFShift + shift_use;
    }

    float shift = FFShift2 - FFShift;

    next.spec1 = locIdx1;
    next.spec2 = locIdx2;
    next.score1 = 1.0;
    next.score2 = 1.0;
    next.shift1 = shift;
    next.shift2 = oriented[locIdx1].parentMass + 100.0;
    outAlign.push_back(next);
    if (debug)
    {
      cout << "Shift: " << locIdx1 << ", " << locIdx2 << " = " << shift << "\n";
    }
  }

  m_params.INPUT_ALIGNS = outAlign;
  m_params.PENALTY_PTM = -200;
  m_params.PENALTY_SAME_VERTEX = -1000000;
  m_params.GRAPH_TYPE = 2;
  m_params.LABELS = (vector<SpectrumPeakLabels> *)0;
  m_params.MAX_AA_JUMP = 1;
  m_params.MAX_MOD_MASS = 100.0;
  m_params.TOLERANCE_PEAK = merging_params->parent_mass_tol;
  m_params.TOLERANCE_PM = merging_params->parent_mass_tol;
  m_params.EDGE_SCORE_TYPE = 1;
  m_params.MIN_MATCHED_PEAKS = merging_params->min_matched_peaks_contig_align;
  m_params.MIN_EDGES_TO_COMPONENT = 0;
  m_params.PATH_MIN_SPECS = 2;
  m_params.PATH_MIN_PEAKS = merging_params->min_matched_peaks_contig_align;
  m_params.SPEC_TYPE_MSMS = 0;

  m_params.NO_SEQUENCING = false;
  m_params.ADD_ENDPOINTS = false;
  m_params.OUTPUT_COMPLETE_ABRUIJN = false;
  ofstream masab_output;

  if (debug)
    Masab(m_params);
  else
    Masab(m_params, masab_output);

  inspec = m_params.OUTPUT_SPECS;
  inabinfo = m_params.COMPONENT_INFO;

  if (inspec.size() == 0 || inspec[0].size() == 0)
  {
    return false;
    /*cout << "\nMASAB RETURNED 0 COMPONENTS, MERGING FOR COMPONENT " << idx
     << " FAILED!\n";
     inspec[0].output(cout);
     putSpec = oriented[idx];
     return TwoValues<float> (0, 0);*/
  }
  else if (inspec.size() > 1)
  {
    cout << "\nCAUGHT MASAB SPLIT for (idx1=" << idx1 << ",idx2=" << idx2
        << ")!!!\n\n";
    return false;
  }

  return true;
}

TwoValues<float> CombineContigs::mergeContigs(int idx,
                                              map<int, map<int, int> >& _graph,
                                              vector<vector<float> >& _edges,
                                              vector<set<int> >& _components,
                                              vector<list<TwoValues<int> > >& _root_alignments,
                                              Spectrum & putSpec,
                                              pair<pair<vector<int> , vector<
                                                  int> > , vector<pair<vector<
                                                  int> , vector<double> > > >& putAb,
                                              bool debug)
{

  Results_PA next;
  vector<Results_PA> outAlign;

  SpecSet inspec;
  abinfo_t inabinfo;

  float minLeftEdgeF = 0;
  float maxRightEdge = oriented[idx].peakList.back()[0];
  map<int, TwoValues<float> > contigIndent;
  masab_parameters m_params;

  //debug = idx == 0;

  if (debug)
    cout << "\nInput spectra to masab for " << idx << ":\n";
  for (set<int>::iterator nodeIt = _components[idx].begin(); nodeIt
      != _components[idx].end(); nodeIt++)
  {
    //node_to_index[*nodeIt] = outs.specs.size();
    //index_to_node[outs.specs.size()] = *nodeIt;
    //outs.specs.push_back(oriented[*nodeIt]);
    float FFShift = (*nodeIt == idx) ? 0 : _edges[_graph[idx][*nodeIt]][2];
    minLeftEdgeF = min(minLeftEdgeF, FFShift);
    float righEdge = FFShift + oriented[*nodeIt].peakList.back()[0];
    maxRightEdge = max(maxRightEdge, righEdge);
    if (debug)
    {
      cout << "\nContig " << *nodeIt << "(shift=" << FFShift << "):\n";
      oriented[*nodeIt].output(cout);
    }
  }

  for (set<int>::iterator nodeIt = _components[idx].begin(); nodeIt
      != _components[idx].end(); nodeIt++)
  {
    //int locIdx = node_to_index[*nodeIt];
    float FFShift = (*nodeIt == idx) ? 0 : _edges[_graph[idx][*nodeIt]][2];
    contigIndent[*nodeIt][0] = FFShift - minLeftEdgeF;
  }

  maxRightEdge -= minLeftEdgeF;
  if (debug)
    cout << "Found " << _root_alignments[idx].size() << " alignments\n";
  for (list<TwoValues<int> >::iterator node_pair_it =
      _root_alignments[idx].begin(); node_pair_it
      != _root_alignments[idx].end(); node_pair_it++)
  {
    int idx1 = (*node_pair_it)[0], idx2 = (*node_pair_it)[1];
    //int locIdx1 = node_to_index[idx1], locIdx2 = node_to_index[idx2];

    float FFShift = (idx1 == idx) ? 0 : _edges[_graph[idx][idx1]][2];
    float FFShift2 = (idx2 == idx) ? 0 : _edges[_graph[idx][idx2]][2];

    float shift = FFShift2 - FFShift;

    next.spec1 = idx1;
    next.spec2 = idx2;
    next.score1 = 1.0;
    next.score2 = 1.0;
    next.shift1 = shift;
    next.shift2 = oriented[idx1].parentMass + 100.0;
    outAlign.push_back(next);
    if (debug)
    {
      cout << "Shift: " << idx1 << ", " << idx2 << " = " << shift << "\n";
    }
  }

  m_params.INPUT_SPECS = oriented;
  m_params.INPUT_ALIGNS = outAlign;
  m_params.PENALTY_PTM = -200;
  m_params.PENALTY_SAME_VERTEX = -1000000;
  m_params.GRAPH_TYPE = 2;
  m_params.LABELS = (vector<SpectrumPeakLabels> *)0;
  m_params.MAX_AA_JUMP = 1;
  m_params.MAX_MOD_MASS = 100.0;
  m_params.TOLERANCE_PEAK = merging_params->parent_mass_tol;
  m_params.TOLERANCE_PM = merging_params->parent_mass_tol;
  m_params.EDGE_SCORE_TYPE = 1;
  m_params.MIN_MATCHED_PEAKS = merging_params->min_matched_peaks_contig_align;
  m_params.MIN_EDGES_TO_COMPONENT = 0;
  m_params.PATH_MIN_SPECS = 2;
  m_params.PATH_MIN_PEAKS = merging_params->min_matched_peaks_contig_align;
  m_params.SPEC_TYPE_MSMS = 0;

  m_params.NO_SEQUENCING = false;
  m_params.ADD_ENDPOINTS = false;
  m_params.OUTPUT_COMPLETE_ABRUIJN = false;
  ofstream masab_output;

  if (debug)
    Masab(m_params);
  else
    Masab(m_params, masab_output);

  inspec = m_params.OUTPUT_SPECS;
  inabinfo = m_params.COMPONENT_INFO;

  if (inspec.size() == 0 || inspec[0].size() == 0)
  {
    cout << "\nMASAB RETURNED 0 COMPONENTS, MERGING FOR COMPONENT " << idx
        << " FAILED!\n";
    inspec[0].output(cout);
    putSpec = oriented[idx];
    return TwoValues<float> (0, 0);
  }
  else if (inspec.size() > 1)
  {
    cout << "\nMASAB SPLIT COMPONENT " << idx << " INTO " << inspec.size()
        << " CONTIGS!\n";
  }
  putSpec = inspec[0];
  putAb = inabinfo[0];

  putSpec.parentMass = putSpec.peakList.back()[0] + AAJumps::massMH;

  /*
   typedef std::map<
   unsigned,  // contig index
   std::pair<
   std::pair< vector<int>, vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
   vector< std::pair< // ABruijn vertices
   vector<int>, vector<double> > // Spectrum index, peak mass
   >
   >
   > abinfo_t;
   */
  pair<vector<int> , vector<double> > first_mass = inabinfo[0].second.front();
  int locIdxF = first_mass.first.front();
  float gapF = first_mass.second.front() + contigIndent[locIdxF][0];

  //Spectrum revPut = putSpec;
  //revPut.reverse(0.0 - AAJumps::massH2O, 0)
  float gapR = maxRightEdge - putSpec.peakList.back()[0] - gapF;// - ((putSpec.parentMass - AAJumps::massHion - putSpec.peakList.back()[0]) - putSpec.peakList.front()[0]);

  /*
   pair< vector<int>, vector<double> > last_mass = inabinfo[0].second.back();
   int locIdxR = last_mass.first.front();
   float gapR = putSpec.parentMass - putSpec.peakList.back()[0] - AAJumps::massHion;// - last_mass.second.front();
   */
  if (debug)
  {
    cout << "\ngapF = " << gapF << " = " << first_mass.second.front() << " + "
        << contigIndent[locIdxF][0] << "\n";
    cout << "gapR = " << gapR << " = " << maxRightEdge << " - "
        << putSpec.peakList.back()[0] << " - " << gapF << "\n";
    cout << "First mass from " << locIdxF << "\n";
    cout << "minLeftEdgeF " << minLeftEdgeF << "\n";
  }

  return TwoValues<float> (gapF, gapR);
}

void CombineContigs::reverseNode(int idx,
                                 map<int, map<int, int> >& _graph,
                                 vector<vector<float> >& _edges,
                                 vector<set<int> >& _components,
                                 vector<bool>& _reversed)
{
  map<int, int>::iterator mapIt;
  float FFShift, FRShift, RFShift, RRShift;

  consensus[idx].reverse(0.0 - AAJumps::massH2O, 0);
  float tempExtr1 = extraShifts[idx][0], tempExtr2 = extraShifts[idx][1];
  extraShifts[idx][0] = tempExtr2;
  extraShifts[idx][1] = tempExtr1;
  oriented[idx].reverse(0.0 - AAJumps::massH2O, 0);
  _reversed[idx] = !_reversed[idx];

  for (mapIt = _graph[idx].begin(); mapIt != _graph[idx].end(); mapIt++)
  {
    int idx2 = mapIt->first;
    int edge2 = mapIt->second;

    if (idx2 == idx)
      continue;

    FFShift = _edges[edge2][2];
    FRShift = _edges[edge2][3];
    RFShift = _edges[edge2][4];
    RRShift = _edges[edge2][5];
    if (_components[idx].count(idx2) > 0)
    {
      _edges[edge2][2] = RRShift;
      _edges[edge2][3] = RFShift;
      _edges[edge2][4] = FRShift;
      _edges[edge2][5] = FFShift;
      oriented[idx2].reverse(0.0 - AAJumps::massH2O, 0);
      consensus[idx2].reverse(0.0 - AAJumps::massH2O, 0);
      _reversed[idx2] = !_reversed[idx2];
    }
    else
    {
      _edges[edge2][2] = RFShift;
      _edges[edge2][3] = RRShift;
      _edges[edge2][4] = FFShift;
      _edges[edge2][5] = FRShift;
      int edgeFrom = _graph[idx2][idx];
      _edges[edgeFrom][2] = 0.0 - RFShift;
      _edges[edgeFrom][3] = 0.0 - FFShift;
      _edges[edgeFrom][4] = 0.0 - RRShift;
      _edges[edgeFrom][5] = 0.0 - FRShift;
    }
  }
}
/*
 void CombineContigs::mergeContigs(Spectrum& spec1, Spectrum& spec2, float consensusShift, Spectrum& putSpec) {
 Results_PA next;
 vector<Results_PA> outAlign;
 SpecSet outs;
 SpecSet inspec;
 outs.specs.push_back(spec1);
 outs.specs.push_back(spec2);
 next.spec1 = 0;
 next.spec2 = 1;
 next.score1 = 1.0;
 next.score2 = 1.0;
 next.shift1 = consensusShift;
 next.shift2 = spec1.parentMass + 100.0;
 outAlign.push_back(next);
 chdir("/home/aguthals/tem p_masab");
 Save_results_bin("sps_merge_aligns.bin", (unsigned int)outAlign.size(), outAlign.begin());
 outs.SaveSpecSet_pklbin("sps_orient.pklbin");
 system("/home/aguthals/cpplib/masab assembly.params");
 //execl( ".", "/home/aguthals/cpplib/masab assembly.params");
 inspec.LoadSpecSet_pklbin("sps_seqs_merged.pklbin");
 chdir("/home/aguthals/merge_cid_hcd_karl_4mp");
 putSpec = inspec[0];
 }
 */
void CombineContigs::getCondensedContigAlignmentGraph(map<int, map<int, int> >& _graph,
                                                      vector<vector<float> >& _edges,
                                                      list<TwoValues<float> >& scoredEdges,
                                                      float minCombScore,
                                                      int minMatchedPeaks)
{
  list<Results_OCC>::iterator pivot_res;
  _graph.clear();
  scoredEdges.clear();
  int edgeIdxUse = _edges.size();
  _edges.resize((contigs.size() * contigs.size()) - contigs.size());
  edge.resize(6);
  TwoValues<int> mp;
  TwoValues<float> res;

  int total = 0;
  int used = 0;

  float FFShift, FRShift, RFShift, RRShift, shift, mRatio, score, reversedShift;
  for (pivot_res = contig_alignments.begin(); pivot_res
      != contig_alignments.end(); pivot_res++)
  {
    int idx1 = (*pivot_res).spec1, idx2 = (*pivot_res).spec2, mPeaks =
        (*pivot_res).matchedPeaks;
    if (haveRes && (idxNotUsed.count(prot_match[idx1][0]) > 0
        || idxNotUsed.count(prot_match[idx2][0]) > 0))
      continue;
    shift = (*pivot_res).shift1, mRatio = min((*pivot_res).score1,
                                              (*pivot_res).score2);
    bool reverse2 = isEqual((*pivot_res).shift2, 1.0, 0.0001);

    total++;
    if (((float)mPeaks) * mRatio < minCombScore || mPeaks < minMatchedPeaks)
      continue;
    //bool correct = validShift(shift, idx1, idx2, reverse2);
    //if (! correct || prot_match[idx1][0] < 0 || prot_match[idx2][0] < 0) continue;

    score = mRatio;
    reversedShift = contigs[idx1].parentMass - (shift
        + contigs[idx2].parentMass);

    used++;

    FFShift = (reverse2) ? contigs[idx1].parentMass : shift;
    FRShift = (reverse2) ? shift : contigs[idx1].parentMass;
    RFShift = (reverse2) ? reversedShift : contigs[idx1].parentMass;
    RRShift = (reverse2) ? contigs[idx1].parentMass : reversedShift;

    edge[0] = (float)idx1;
    edge[1] = (float)idx2;
    edge[2] = FFShift;
    edge[3] = FRShift;
    edge[4] = RFShift;
    edge[5] = RRShift;
    scoredEdge[0] = (float)edgeIdxUse;
    scoredEdge[1] = (score * ((float)mPeaks));
    if (_graph.count(idx1) == 0)
    {
      edgeMap.clear();
      edgeMap[idx2] = edgeIdxUse;
      _graph[idx1] = edgeMap;
    }
    else
      _graph[idx1][idx2] = edgeIdxUse;

    _edges[edgeIdxUse] = edge;
    scoredEdges.push_back(scoredEdge);

    if ((idx1 == 67 && idx2 == 712) || (idx1 == 712 && idx2 == 67))
    {
      cout << "adding edge " << idx1 << " to " << idx2 << ", edgId x == "
          << edgeIdxUse << ", score = " << scoredEdge[1] << endl;
    }

    edgeIdxUse++;

    FFShift = (reverse2) ? contigs[idx2].parentMass : 0.0 - shift;
    FRShift = (reverse2) ? 0.0 - reversedShift : contigs[idx2].parentMass;
    RFShift = (reverse2) ? 0.0 - shift : contigs[idx2].parentMass;
    RRShift = (reverse2) ? contigs[idx2].parentMass : 0.0 - reversedShift;

    edge[0] = (float)idx2;
    edge[1] = (float)idx1;
    edge[2] = FFShift;
    edge[3] = FRShift;
    edge[4] = RFShift;
    edge[5] = RRShift;
    if (_graph.count(idx2) == 0)
    {
      edgeMap.clear();
      edgeMap[idx1] = edgeIdxUse;
      _graph[idx2] = edgeMap;
    }
    else
      _graph[idx2][idx1] = edgeIdxUse;

    _edges[edgeIdxUse] = edge;
    edgeIdxUse++;

  }
  _edges.resize(edgeIdxUse);
  scoredEdges.sort(SortEdges());
  cout << "\nUsed " << used << " of " << total << " contig-contig shifts\n";
}
/*
 graph: idx1 --> idx2 --> edgeIdx (undirected)
 edges: edgeIdx --> (idx1, idx2, FFShift, FRShift, RFShift, RRShift)
 scoredEdges: list of (edgeIdx, score) sorted in by decreasing score
 */

void CombineContigs::getCondensedAlignmentGraph(vector<set<int> >& _components,
                                                map<int, map<int, int> >& _graph,
                                                vector<vector<float> >& _edges,
                                                vector<bool>& _reversed,
                                                list<TwoValues<float> >& scoredEdges,
                                                vector<list<TwoValues<int> > >& _root_alignments,
                                                short precCharge,
                                                float minRatioConnContig,
                                                float minRatio,
                                                int minNumMatchedPeaks,
                                                int minNumConnectors,
                                                int minEdgesToComponent)
{

  //float minComponentOverlapScore = 0.2;
  //float minContigOverlapScore = 0.3;
  //int minComponentSize = 3;

  int edgeIdxUse = _edges.size();
  _edges.resize((contigs.size() * contigs.size()) - contigs.size());
  edge.resize(6);

  vector<float> edgeHold(7);
  edgeNum.resize(_edges.size());
  scoredEdges.clear();

  map<unsigned int, map<unsigned int, Results_OCC> > connContRes;
  map<unsigned int, map<unsigned int, list<unsigned int> > > graphCont;
  map<unsigned int, list<unsigned int> >::iterator graphContIt;
  list<unsigned int>::iterator nodeIt;
  TwoValues<float> res;
  TwoValues<int> mp;

  getAlignmentGraph(&graphCont,
                    &connContRes,
                    precCharge,
                    minNumMatchedPeaks,
                    minRatio,
                    minRatioConnContig);

  vector<int> componentRef(contigs.size());
  //int numberComponents = 0;
  for (int i = 0; i < _components.size(); i++)
  {
    componentRef[i] = -1;
  }
  for (int i = 0; i < _components.size(); i++)
  {
    //if (_components[i].size() >= minComponentSize) {
    //numberComponents ++;
    for (set<int>::iterator nodeIt = _components[i].begin(); nodeIt
        != _components[i].end(); nodeIt++)
    {
      componentRef[*nodeIt] = i;
    }
    //}
  }
  map<int, vector<float> >* tempEdges = new map<int, vector<float> > ;
  map<int, set<int> >* edgeMPs = new map<int, set<int> > ;
  int tempEdgeIdx = 0;
  map<int, map<int, list<int> > >* mpConnect = new map<int,
      map<int, list<int> > > ;
  map<int, map<int, list<int> > >::iterator mpConnectIt;
  map<int, list<int> > mpConnectDest;
  map<int, list<int> >::iterator mpConnectDestIt;
  list<int> mpConnectEdges;
  list<int>::iterator mpConnectEdgesIt;

  map<unsigned int, set<unsigned int> > used;
  map<int, float> shiftScore;
  map<int, vector<float> > shiftCount;
  map<int, set<int> > shiftMPs;
  vector<float> pmCount;
  map<int, float>::iterator shiftScoreIt;
  //int minSize = 10;

  //FILE* output = fopen("debug_connect.txt", "w");
  int numConn = 0;
  for (int contig = 0; contig < contigs.size(); contig++)
  {
    if (graphCont.count(contig) == 0)
      continue;
    if (haveRes && (idxNotUsed.count(prot_match[contig][0]) > 0))
      continue;
    for (graphContIt = graphCont[contig].begin(); graphContIt
        != graphCont[contig].end(); graphContIt++)
    {
      int idx1 = contig, idx2 = graphContIt->first;
      int comp1 = componentRef[idx1], comp2 = componentRef[idx2];
      if (haveRes && (idxNotUsed.count(prot_match[idx2][0]) > 0))
        continue;
      if (used.count(idx1) > 0 && used[idx1].count(idx2) > 0)
        continue;
      if (used.count(idx2) > 0 && used[idx2].count(idx1) > 0)
        continue;
      if (used.count(idx1) == 0)
      {
        set<unsigned int> myS;
        myS.insert(idx2);
        used[idx1] = myS;
      }
      else
        used[idx1].insert(idx2);
      bool debug = false;//(idx1 == 46 && idx2 == 517) || (idx1 == 517 && idx2 == 46);//(comp1 == 454 && comp2 == 344) || (comp1 == 344 && comp2 == 454)
      if (debug)
        cout << "found indicies " << idx1 << "(" << comp1 << ") and " << idx2
            << "(" << comp2 << ")\n";
      if (comp1 == comp2 || comp1 < 0 || comp2 < 0)
        continue;
      if (debug)
        cout << "proceeding\n";
      //if (_components[comp1].size() < minComponentSize || _components[comp2].size() < minComponentSize) continue;
      //if ((graphContIt->second).size() < minNumConnectors) continue;
      //if (contigs[idx1].size() < minSize || contigs[idx2].size() < minSize) continue;
      //if (_graph.count(idx1) > 0 && _graph[idx1].count(idx2) > 0) continue;
      shiftScore.clear();
      shiftCount.clear();
      shiftMPs.clear();
      map<int, list<float> > shiftScore2;
      map<int, list<float> > shiftScore1;

      float bestShift;
      float score;
      int numShifts;
      set<int> connIdxs;
      for (nodeIt = (graphContIt->second).begin(); nodeIt
          != (graphContIt->second).end(); nodeIt++)
      {
        Results_OCC res1 = connContRes[*nodeIt][idx1];
        Results_OCC res2 = connContRes[*nodeIt][idx2];
        int connIdx = *nodeIt;
        float minScore = min(res1.score1, res1.score2);
        minScore = min(minScore, res2.score1);
        minScore = min(minScore, res2.score2);
        float fShift = res2.shift1 - res1.shift1;
        bestShift = fShift;
        score = minScore;
        int intShift = floatToInt(fShift / InputParams::Resolution);
        for (int s = intShift - intParentMassTol; s <= intShift
            + intParentMassTol; s++)
        {
          if (shiftScore.count(s) == 0)
          {
            shiftScore[s] = minScore - (0.0001 * (float)abs(s - intShift));
            pmCount.clear();
            pmCount.push_back(connectors[connIdx].parentMass);
            shiftCount[s] = pmCount;
            set<int> mpIdxs;
            mpIdxs.insert(connIdx);
            shiftMPs[s] = mpIdxs;
            list<float> Score2;
            Score2.push_back(min(res1.score2, res2.score2));
            shiftScore2[s] = Score2;
            list<float> Score1;
            Score1.push_back(min(res1.score1, res2.score1));
            shiftScore1[s] = Score1;
          }
          else
          {
            shiftScore[s] += minScore - (0.0001 * (float)abs(s - intShift));
            shiftCount[s].push_back(connectors[connIdx].parentMass);
            shiftMPs[s].insert(connIdx);
            shiftScore2[s].push_back(min(res1.score2, res2.score2));
            shiftScore1[s].push_back(min(res1.score1, res2.score1));
          }
        }
      }
      score = 0;
      list<float> bestScore2s;
      list<float> bestScore1s;
      for (shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt++)
      {
        float shiftScore = shiftScoreIt->second;
        if (shiftScore > score)
        {
          score = shiftScore;
          bestShift = ((float)shiftScoreIt->first) * InputParams::Resolution;
          pmCount = shiftCount[shiftScoreIt->first];
          connIdxs = shiftMPs[shiftScoreIt->first];
          bestScore2s = shiftScore2[shiftScoreIt->first];
          bestScore1s = shiftScore1[shiftScoreIt->first];
        }
      }
      sort(pmCount.begin(), pmCount.end());
      bool first = true;
      float lastPM;
      numShifts = 0;
      for (vector<float>::iterator Lit = pmCount.begin(); Lit != pmCount.end(); Lit++)
      {
        float pm = *Lit;
        if (first || lastPM + minParentMassSep <= pm)
        {
          first = false;
          lastPM = pm;
          numShifts++;
        }
      }
      //if (numShifts < minNumConnectors) continue;
      //if (idx1 == 53 && idx2 == 1) cout << idx1 << " " << idx2 << " " << bestShift << "\n";
      first = true;
      float FFShift, FRShift, RFShift, RRShift, dist;
      for (nodeIt = (graphContIt->second).begin(); nodeIt
          != (graphContIt->second).end(); nodeIt++)
      {
        Results_OCC res1 = connContRes[*nodeIt][idx1];
        Results_OCC res2 = connContRes[*nodeIt][idx2];
        float fShift = res2.shift1 - res1.shift1;
        if (first || abs(fShift - bestShift) < dist)
        {
          FFShift = fShift;
          FRShift = res2.shift2 - res1.shift1;
          RFShift = res2.shift1 - res1.shift2;
          RRShift = res2.shift2 - res1.shift2;
          dist = abs(fShift - bestShift);
          first = false;
        }
      }
      float contig2PM = contigs[idx2].parentMass;

      float TempFFShift = FFShift, TempFRShift = FRShift,
          TempRFShift = RFShift, TempRRShift = RRShift;
      if (_reversed[idx1] && !_reversed[idx2])
      {
        FFShift = TempRFShift;
        FRShift = TempRRShift;
        RFShift = TempFFShift;
        RRShift = TempFRShift;
      }
      else if (!_reversed[idx1] && _reversed[idx2])
      {
        FFShift = TempFRShift;
        FRShift = TempFFShift;
        RFShift = TempRRShift;
        RRShift = TempRFShift;
      }
      else if (_reversed[idx1] && _reversed[idx2])
      {
        FFShift = TempRRShift;
        FRShift = TempRFShift;
        RFShift = TempFRShift;
        RRShift = TempFFShift;
      }
      /*
       TwoValues<float> ratios = getContigOverlapScores(oriented[idx1], oriented[idx2], FFShift, FRShift, parentMassTol, minContigOverlap, false);
       TwoValues<int> mps = getContigOverlapPeaks(oriented[idx1], oriented[idx2], FFShift, FRShift, parentMassTol, minContigOverlap, false);
       int fr = (ratios[0] > ratios[1]) ? 0 : 1;
       if (mps[fr] >= consecMP && ratios[fr] > -0.01 && ratios[fr] < minContigOverlapScore) continue;
       */

      if (debug)
      {
        bool reverse1 = prot_match[idx1][2] == 1;
        bool reverse2 = prot_match[idx2][2] == 1;
        float protShift1 = getContigShift(idx1, reverse1);
        float protShift2 = getContigShift(idx2, reverse2);
        float protFFShift = protShift2 - protShift1;
        float protRRShift = (protShift1 + contigs[idx1].parentMass)
            - (protShift2 + contigs[idx2].parentMass);
        //protFFShift = (reverse2 == _reversed[idx2]) ? protFFShift : protRRShift;
        cout << "\ncorrect shift between " << idx1 << " and " << idx2 << ": "
            << protFFShift << ", " << protRRShift << "\n";
        cout << "all shifts observed: " << TempFFShift << ", " << TempFRShift
            << ", " << TempRFShift << ", " << TempRRShift << "\n";
        for (list<float>::iterator sc2it = bestScore1s.begin(); sc2it
            != bestScore1s.end(); sc2it++)
        {
          cout << *sc2it << ", ";
        }
        cout << "\n";
        for (list<float>::iterator sc2it = bestScore2s.begin(); sc2it
            != bestScore2s.end(); sc2it++)
        {
          cout << *sc2it << ", ";
        }
        cout << "\n";
        for (set<int>::iterator idxit = connIdxs.begin(); idxit
            != connIdxs.end(); idxit++)
        {
          cout << *idxit << ", ";
        }
        cout << "\n";
        cout << "adding edge between " << comp1 << "(" << prot_match[comp1][2]
            << _reversed[comp1] << ") and " << comp2 << "("
            << prot_match[comp2][2] << _reversed[comp2] << ") - (" << idx1
            << "(" << prot_match[idx1][2] << _reversed[idx1] << "), " << idx2
            << "(" << prot_match[idx2][2] << _reversed[idx2] << ")): "
            << FFShift << ", " << FRShift << "\n";
        cout << FFShift << " + " << _edges[_graph[comp1][idx1]][2] << " - "
            << _edges[_graph[comp2][idx2]][2] << " = " << FFShift
            + _edges[_graph[comp1][idx1]][2] - _edges[_graph[comp2][idx2]][2]
            << "\n";
        cout << FRShift << " + " << _edges[_graph[comp1][idx1]][2] << " - "
            << _edges[_graph[comp2][idx2]][5] << " = " << FRShift
            + _edges[_graph[comp1][idx1]][2] - _edges[_graph[comp2][idx2]][5]
            << "\n";
      }
      float locFFShift = (idx1 == comp1) ? 0 : _edges[_graph[comp1][idx1]][2];
      float locRRShift = (idx1 == comp1) ? 0 : _edges[_graph[comp1][idx1]][5];
      float locFFShift2 = (idx2 == comp2) ? 0 : _edges[_graph[comp2][idx2]][2];
      float locRRShift2 = (idx2 == comp2) ? 0 : _edges[_graph[comp2][idx2]][5];

      FFShift += locFFShift - locFFShift2;
      FRShift += locFFShift - locRRShift2;
      RFShift += locRRShift - locFFShift2;
      RRShift += locRRShift - locRRShift2;
      //TwoValues<float> overlapScores = getOverlapScore(comp1, comp2, FFShift, FRShift, consecMP, _graph, _edges, _components);
      //float overlapScore = max(overlapScores[0], overlapScores[1]);
      //if (overlapScore < minComponentOverlapScore) continue;

      edgeHold[0] = (float)comp1;
      edgeHold[1] = (float)comp2;
      edgeHold[2] = FFShift;
      edgeHold[3] = FRShift;
      edgeHold[4] = RFShift;
      edgeHold[5] = RRShift;
      edgeHold[6] = score;

      (*edgeMPs)[tempEdgeIdx] = connIdxs;
      (*tempEdges)[tempEdgeIdx] = edgeHold;

      mpConnectEdges.clear();
      mpConnectEdges.push_back(tempEdgeIdx);
      mpConnectDest.clear();
      mpConnectDest[comp2] = mpConnectEdges;

      if ((*mpConnect).count(comp1) == 0)
      {
        (*mpConnect)[comp1] = mpConnectDest;
      }
      else if ((*mpConnect)[comp1].count(comp2) == 0)
      {
        (*mpConnect)[comp1][comp2] = mpConnectEdges;
      }
      else
      {
        (*mpConnect)[comp1][comp2].push_back(tempEdgeIdx);
      }

      tempEdgeIdx++;

      edgeHold[0] = (float)comp2;
      edgeHold[1] = (float)comp1;
      edgeHold[2] = 0.0 - FFShift;
      edgeHold[3] = 0.0 - RFShift;
      edgeHold[4] = 0.0 - FRShift;
      edgeHold[5] = 0.0 - RRShift;
      edgeHold[6] = score;
      (*edgeMPs)[tempEdgeIdx] = connIdxs;
      (*tempEdges)[tempEdgeIdx] = edgeHold;

      mpConnectEdges.clear();
      mpConnectEdges.push_back(tempEdgeIdx);
      mpConnectDest.clear();
      mpConnectDest[comp1] = mpConnectEdges;

      if ((*mpConnect).count(comp2) == 0)
      {
        (*mpConnect)[comp2] = mpConnectDest;
      }
      else if ((*mpConnect)[comp2].count(comp1) == 0)
      {
        (*mpConnect)[comp2][comp1] = mpConnectEdges;
      }
      else
      {
        (*mpConnect)[comp2][comp1].push_back(tempEdgeIdx);
      }

      tempEdgeIdx++;
    }
  }
  used.clear();
  //shiftScore
  map<int, int> shiftCountMP;
  map<int, set<int> > shiftCountMPIdx;
  map<int, vector<float> > shiftEdge;

  map<int, float> shiftScoreR;
  map<int, int> shiftCountMPR;
  map<int, set<int> > shiftCountMPIdxR;
  map<int, vector<float> > shiftEdgeR;

  for (mpConnectIt = (*mpConnect).begin(); mpConnectIt != (*mpConnect).end(); mpConnectIt++)
  {
    int comp1 = mpConnectIt->first;
    for (mpConnectDestIt = mpConnectIt->second.begin(); mpConnectDestIt
        != mpConnectIt->second.end(); mpConnectDestIt++)
    {
      int comp2 = mpConnectDestIt->first;
      if (used.count(comp1) > 0 && used[comp1].count(comp2) > 0)
        continue;
      if (used.count(comp2) > 0 && used[comp2].count(comp1) > 0)
        continue;
      if (used.count(comp1) == 0)
      {
        set<unsigned int> myS;
        myS.insert(comp2);
        used[comp1] = myS;
      }
      else
        used[comp1].insert(comp2);

      shiftScore.clear();
      shiftCountMP.clear();
      shiftCountMPIdx.clear();
      shiftEdge.clear();
      shiftScoreR.clear();
      shiftCountMPR.clear();
      shiftCountMPIdxR.clear();
      shiftEdgeR.clear();

      bool debug = false;//(comp1 == 454 && comp2 == 344) || (comp1 == 344 && comp2 == 454);//(comp1 == 454 && comp2 == 344) || (comp1 == 344 && comp2 == 454);//(comp1 == 16 && comp2 == 26) || (comp1 == 26 && comp2 == 16);

      if (debug)
      {
        cout << "Counting edges between " << comp1 << " and " << comp2
            << "\nFFShifts: ";
      }
      for (mpConnectEdgesIt = mpConnectDestIt->second.begin(); mpConnectEdgesIt
          != mpConnectDestIt->second.end(); mpConnectEdgesIt++)
      {
        int edgeIdx = *mpConnectEdgesIt;
        edgeHold = (*tempEdges)[edgeIdx];
        if (debug)
        {
          cout << edgeHold[2] << ", ";
        }
        int intShift = floatToInt(edgeHold[2] / InputParams::Resolution);
        for (int s = intShift - intParentMassTol; s <= intShift
            + intParentMassTol; s++)
        {
          if (shiftScore.count(s) == 0)
          {
            shiftScore[s] = edgeHold[6] - (0.0001 * (float)abs(s - intShift));
            shiftCountMP[s] = 1;
            shiftEdge[s] = edgeHold;
            shiftCountMPIdx[s] = (*edgeMPs)[edgeIdx];
          }
          else
          {
            shiftScore[s] += edgeHold[6] - (0.0001 * (float)abs(s - intShift));
            shiftCountMP[s] += 1;
            shiftCountMPIdx[s].insert((*edgeMPs)[edgeIdx].begin(),
                                      (*edgeMPs)[edgeIdx].end());
          }
        }
      }
      if (debug)
      {
        cout << "\nFRShifts: ";
      }
      for (mpConnectEdgesIt = mpConnectDestIt->second.begin(); mpConnectEdgesIt
          != mpConnectDestIt->second.end(); mpConnectEdgesIt++)
      {
        int edgeIdx = *mpConnectEdgesIt;
        edgeHold = (*tempEdges)[edgeIdx];
        if (debug)
        {
          cout << edgeHold[3] << ", ";
        }
        int intShift = floatToInt(edgeHold[3] / InputParams::Resolution);
        for (int s = intShift - intParentMassTol; s <= intShift
            + intParentMassTol; s++)
        {
          if (shiftScoreR.count(s) == 0)
          {
            shiftScoreR[s] = edgeHold[6] - (0.0001 * (float)abs(s - intShift));
            shiftCountMPR[s] = 1;
            shiftEdgeR[s] = edgeHold;
            shiftCountMPIdxR[s] = (*edgeMPs)[edgeIdx];
          }
          else
          {
            shiftScoreR[s] += edgeHold[6] - (0.0001 * (float)abs(s - intShift));
            shiftCountMPR[s] += 1;
            shiftCountMPIdxR[s].insert((*edgeMPs)[edgeIdx].begin(),
                                       (*edgeMPs)[edgeIdx].end());
          }
        }
      }

      float bestScore = 0;
      int bestCount = 0;
      int bestMPCount = 0;
      vector<float> bestEdge;
      for (shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt++)
      {
        float shiftScore = shiftScoreIt->second;
        int shiftCount = shiftCountMP[shiftScoreIt->first];
        int shiftMPCount = shiftCountMPIdx[shiftScoreIt->first].size();
        if (shiftScore > bestScore && shiftMPCount >= bestMPCount)
        {
          bestScore = shiftScore;
          bestCount = shiftCount;
          bestMPCount = shiftMPCount;
          bestEdge = shiftEdge[shiftScoreIt->first];
        }
      }

      for (shiftScoreIt = shiftScoreR.begin(); shiftScoreIt
          != shiftScoreR.end(); shiftScoreIt++)
      {
        float shiftScore = shiftScoreIt->second;
        int shiftCount = shiftCountMPR[shiftScoreIt->first];
        int shiftMPCount = shiftCountMPIdxR[shiftScoreIt->first].size();
        if (shiftScore > bestScore && shiftMPCount >= bestMPCount)
        {
          bestScore = shiftScore;
          bestCount = shiftCount;
          bestMPCount = shiftMPCount;
          bestEdge = shiftEdgeR[shiftScoreIt->first];
        }
      }

      if (bestMPCount < minNumConnectors)
      {
        continue;
      }
      if (bestCount < minEdgesToComponent)
      {
        continue;
      }

      numConn++;
      edge[0] = bestEdge[0];
      edge[1] = bestEdge[1];
      edge[2] = bestEdge[2];
      edge[3] = bestEdge[3];
      edge[4] = bestEdge[4];
      edge[5] = bestEdge[5];
      if (debug)
      {
        cout << "\nbest edge: " << edge[2] << ", " << edge[3] << ", "
            << edge[4] << ", " << edge[5] << "\n";
        cout << "num connectors: " << bestMPCount << ", edges to component: "
            << bestCount << "\n";
      }
      //getConsensusOverlap(idx1, idx2, fShift, rShift, consecMP, _graph, _edges, _components, _root_alignments, debug);
      //TwoValues<float> overlapScores = getOverlapScore(comp1, comp2, edge[2], edge[3], consecMP, _graph, _edges, _components);
      TwoValues<float> overlapScores = getConsensusOverlap(comp1,
                                                           comp2,
                                                           edge[2],
                                                           edge[3],
                                                           _graph,
                                                           _edges,
                                                           _components,
                                                           _root_alignments);
      float overlapScore = max(overlapScores[0], overlapScores[1]);

      scoredEdge[0] = (float)edgeIdxUse;
      scoredEdge[1] = overlapScore;
      if (_graph.count(comp1) == 0)
      {
        edgeMap.clear();
        edgeMap[comp2] = edgeIdxUse;
        _graph[comp1] = edgeMap;
      }
      else
        _graph[comp1][comp2] = edgeIdxUse;

      _edges[edgeIdxUse] = edge;

      edgeNum[edgeIdxUse] = 2 * bestMPCount;
      scoredEdges.push_back(scoredEdge);

      edgeIdxUse++;

      edge[0] = bestEdge[1];
      edge[1] = bestEdge[0];
      edge[2] = 0.0 - bestEdge[2];
      edge[3] = 0.0 - bestEdge[4];
      edge[4] = 0.0 - bestEdge[3];
      edge[5] = 0.0 - bestEdge[5];
      if (_graph.count(comp1) == 0)
      {
        edgeMap.clear();
        edgeMap[comp1] = edgeIdxUse;
        _graph[comp2] = edgeMap;
      }
      else
        _graph[comp2][comp1] = edgeIdxUse;

      _edges[edgeIdxUse] = edge;
      edgeNum[edgeIdxUse] = 2 * bestMPCount;
      edgeIdxUse++;
    }
  }
  //fclose(output);
  _edges.resize(edgeIdxUse);
  edgeNum.resize(edgeIdxUse);
  scoredEdges.sort(SortEdges());
  delete mpConnect;
  delete edgeMPs;
  delete tempEdges;
  cout << "Added " << numConn
      << " contig-contig edges from connector-contig shifts.\n";
}

bool CombineContigs::validShift(float FFShift,
                                int idx1,
                                int idx2,
                                bool isReversed)
{
  bool reverse1 = prot_match[idx1][2] == 1;
  bool reverse2 = prot_match[idx2][2] == 1;
  float protShift1 = getContigShift(idx1, reverse1);
  float protShift2 = getContigShift(idx2, reverse2);
  float protFFShift = protShift2 - protShift1;
  float protRRShift = (protShift1 + contigs[idx1].parentMass) - (protShift2
      + contigs[idx2].parentMass);
  protFFShift = (reverse2 == isReversed) ? protFFShift : protRRShift;
  bool correct = isEqual(FFShift, protFFShift, AAJumps::massHion
      + merging_params->parent_mass_tol) && prot_match[idx1][0]
      == prot_match[idx2][0];
  return correct;
}

bool CombineContigs::validShiftMod(float FFShift,
                                   int idx1,
                                   int idx2,
                                   bool isReversed)
{
  bool reverse1 = prot_match[idx1][2] == 1;
  bool reverse2 = prot_match[idx2][2] == 1;
  float protShift1 = getContigShift(idx1, reverse1);
  float protShift2 = getContigShift(idx2, reverse2);
  float protFFShift = protShift2 - protShift1;
  float protRRShift = (protShift1 + contigs[idx1].parentMass) - (protShift2
      + contigs[idx2].parentMass);
  protFFShift = (reverse2 == isReversed) ? protFFShift : protRRShift;
  bool correct = isEqual(FFShift, protFFShift, 58.0) && prot_match[idx1][0]
      == prot_match[idx2][0];
  return correct;
}

float CombineContigs::getContigShift(int index, bool reverse)
{
  if (!haveRes)
    return 0;
  unsigned int pepIndex = prot_match[index][0];
  char* peptide = fasta.getSequence(pepIndex);
  Spectrum masses = fasta.getMassesSpec(pepIndex);
  float pmTol = parentMassTol;
  Spectrum overlap = overlaps[index];
  Spectrum contig_spec = contigs[index];
  if (reverse)
    contig_spec.reverse(0.0 - AAJumps::massH2O, 0);

  int end = (int)(overlap[overlap.size() - 1][1] + 0.01);
  list<int>::iterator start_other_pivot;
  map<int, float> shiftScore;
  map<int, float>::iterator shiftScoreIt;
  float best_shift1;
  float shift;
  float intensity;
  char sequence[strlen(peptide)];
  strcpy(sequence, peptide);
  int intShift;
  int intPmTol = intParentMassTol;

  for (int i = 0; i < overlap.size(); i++)
  {
    int cIdx = floatToInt(overlap[i][0]);
    int pIdx = floatToInt(overlap[i][1]);
    shift = masses[pIdx][0] - contig_spec[cIdx][0];
    //if (index == 98) cout << cIdx << " : " << contig_spec[cIdx][0] << ", " << pIdx << " : " << masses[pIdx][0] << " = " << shift << "\n";
    intShift = floatToInt(shift / InputParams::Resolution);
    intensity = contig_spec[cIdx][1];

    for (int idx = intShift - intPmTol; idx <= intShift + intPmTol; idx++)
    {
      if (shiftScore.count(idx) == 0)
        shiftScore[idx] = intensity - (0.0001 * (float)abs(intShift - idx));
      else
        shiftScore[idx] += intensity - (0.0001 * (float)abs(intShift - idx));
    }
  }
  float maxScore = 0.0;
  for (shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt++)
  {
    if (shiftScoreIt->second > maxScore)
    {
      maxScore = shiftScoreIt->second;
      best_shift1 = ((float)shiftScoreIt->first) * InputParams::Resolution;
    }
  }
  return best_shift1;
}

void CombineContigs::getContigDistances(map<int, map<int, float> >& contigContigShifts)
{
  if (!haveRes)
    return;
  contigContigShifts.clear();
  map<int, float> contigShifts;
  map<int, float> computedShifts;
  for (int i = 0; i < contigs.size(); i++)
  {
    int pepIndex = prot_match[i][0];

    float shift1;
    if (computedShifts.count(i) == 0)
    {
      shift1 = getContigShift(i, prot_match[i][2] == 1);
      computedShifts[i] = shift1;
    }
    else
      shift1 = computedShifts[i];

    contigShifts.clear();
    contigContigShifts[i] = contigShifts;

    for (int j = 0; j < contigs.size(); j++)
    {
      if (i == j || pepIndex != prot_match[j][0])
        continue;
      float shift2;
      if (computedShifts.count(j) == 0)
      {
        shift2 = getContigShift(j, prot_match[j][2] == 1);
        computedShifts[j] = shift2;
      }
      else
        shift2 = computedShifts[j];

      contigContigShifts[i][j] = shift2 - shift1;
    }
  }
}
