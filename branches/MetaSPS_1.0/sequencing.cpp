#include "sequencing.h"

// Items that need fixing to avoid having to split aligns into ASP/PA
// - MSGraph g;   g.build(alignsASP);   g.add(cAligns[cIdx]); - should be able to initialize using alignsPA
// - VertexSet::addEndpointEdges - only input is alignsASP?
// - SetMerger::splitSet(cIdx,usedSpecs) - ASP/PA cAligns/cAligns_idx are class members

void SplitASPPA_v2(vector<Results_PA> &pairs, float pmTol, vector<Results_ASP> &pairsASP, vector<Results_PA> &pairsPA) {
	pairsASP.resize(pairs.size());   pairsPA.resize(pairs.size());
	unsigned int idxPair, idxASP=0, idxPA=0;
	for(idxPair=0; idxPair<pairs.size(); idxPair++)
		if(abs(pairs[idxPair].shift1)<=pmTol or abs(pairs[idxPair].shift2)<=pmTol)
			pairsASP[idxASP++]=pairs[idxPair]; else pairsPA[idxPA++]=pairs[idxPair];
	pairsASP.resize(idxASP);   pairsPA.resize(idxPA);
}


/*
struct masab_parameters {
	SpecSet INPUT_SPECS;
	vector<Results_PA> INPUT_ALIGNS;
	ostream STDOUT;
	
	float PENALTY_PTM;
	float PENALTY_SAME_VERTEX;
	int GRAPH_TYPE;

	vector<SpectrumPeakLabels>* LABELS;
	int MAX_AA_JUMP;
	float TOLERANCE_PEAK;
	float TOLERANCE_PM;
	short EDGE_SCORE_TYPE;
	unsigned int MIN_MATCHED_PEAKS;
	unsigned int MIN_EDGES_TO_COMPONENT;
	unsigned int PATH_MIN_SPECS;
	short PATH_MIN_PEAKS;
	int SPEC_TYPE_MSMS;
	float ION_OFFSET;

	bool NO_SEQUENCING;
	bool ADD_ENDPOINTS;

	SetMerger COMPONENTS;
	vector<vector<float> > COMPONENT_STATS;
	vector<list<int> > COMPONENT_SPECTRA;
	abinfo_t COMPONENT_INFO;
	abinfo_t COMPLETE_ABRUIJN;
	Clusters PATH_SPECTRA_AS_CLUSTER;
	vector<vector<short> > ABCOUNTS;
	SpecSet OUTPUT_SPECS;
};
*/

void Masab(masab_parameters& m_parameters, ostream& STDOUT) {
	
	vector<SpectrumPeakLabels>* labelsP = m_parameters.LABELS;
	float penalty_ptm = m_parameters.PENALTY_PTM;
	float penalty_sameVert = m_parameters.PENALTY_SAME_VERTEX;
	int graphType = m_parameters.GRAPH_TYPE;
	
	int maxAAjump = m_parameters.MAX_AA_JUMP;
	float maxModMass = m_parameters.MAX_MOD_MASS;
	float peakTol = m_parameters.TOLERANCE_PEAK;
	float pmTol = m_parameters.TOLERANCE_PM;
    short edgeScoreType = m_parameters.EDGE_SCORE_TYPE;
    unsigned int minMatchedPeaks = m_parameters.MIN_MATCHED_PEAKS;
    unsigned int minEdgesToComponent = m_parameters.MIN_EDGES_TO_COMPONENT;
    unsigned int pathMinSpecs = m_parameters.PATH_MIN_SPECS;
    short pathMinPeaks = m_parameters.PATH_MIN_PEAKS;
    int specType = m_parameters.SPEC_TYPE_MSMS;
	float ionOffset = specType?AAJumps::massHion:0;

	bool noSequencing = m_parameters.NO_SEQUENCING;
	bool addEndpoints = m_parameters.ADD_ENDPOINTS;
	bool wholeABFN = m_parameters.OUTPUT_COMPLETE_ABRUIJN;

	SpecSet specSet = m_parameters.INPUT_SPECS;

	vector<Results_PA> aligns = m_parameters.INPUT_ALIGNS;
	vector<Results_PA> tmpAlignsPA;
	vector<Results_ASP> tmpAlignsASP;
	
	//STDOUT << penalty_ptm << " " << penalty_sameVert << " " << graphType << " " << maxAAjump << " " << peakTol << " " << pmTol << " " << edgeScoreType << " " << minMatchedPeaks << " " << minEdgesToComponent << " " << pathMinSpecs << " " << pathMinPeaks << " " << specType << " " << ionOffset << " " << noSequencing << " " << addEndpoints << " " << wholeABFN << "\n";
	
	if (addEndpoints) {
	    // Make sure every spectrum has zero/parentMass-19 nodes
	    if(labelsP) for(unsigned int i=0; i<specSet.size(); i++) specSet[i].addZPMpeaks(peakTol,true,&(*labelsP)[i]);
	    else for(unsigned int i=0; i<specSet.size(); i++) specSet[i].addZPMpeaks(peakTol,ionOffset,true);
	  }

	  // Separate aligns into components
	//  vector<vector<Results_PA> > cAligns;    vector<vector<int> > cAligns_idx;  // To allow going back from components to original order
	  SetMerger components(specSet.size());
	  components.createSets(specSet.size(),2,tmpAlignsASP,aligns);
	  components.splitAligns(tmpAlignsASP,aligns);
	//  components.splitAligns(aligns,cAligns,cAligns_idx);
	  STDOUT<<"Got "<<components.size()<<" component(s).\n";

	  vector<MSGraph> spectrumGraphs(specSet.size());   // Used to build an MSGraph for each spectrum
	  AAJumps jumps2(2);

	  // Split aligns per component and process each component
	  char sBuf[1024]; AAJumps jumps(1);
	  unsigned int numElemsInSets=components.numElemsInSets();  // Maximum possible number of components (all singletons)
	  vector<vector<float> > cStats(numElemsInSets);  for(unsigned int cIdx=0;cIdx<cStats.size();cIdx++) { cStats[cIdx].resize(9); for(unsigned int i=0;i<9;i++) cStats[cIdx][i]=0; }  // Absolute value 9 depends on MSGraph::info_heaviestPath
	  Clusters           pathSpectra;   pathSpectra.resize(numElemsInSets);  // Keep the de-novo reconstructed heaviestPath sequences as spectra in a Cluster variable
	  vector<list<int> > cSpectra(numElemsInSets);  // Lists of spectrum indices per component
	  vector<bool>       specFlipped(specSet.size()); for(unsigned int i=0;i<specFlipped.size();i++) specFlipped[i]=false;

	  // Keeps track of which spectrum peaks were matched (dim.3) in each ABruijn vertex (dim.2)
	  //   from the de novo sequence (heaviest path) in each component (dim.1)
	  vector<vector<list<TwoValues<int> > > > abVertices(numElemsInSets);
	  for(unsigned int i=0;i<abVertices.size();i++) abVertices[i].resize(0);

	  // Keeps track of which spectrum peaks were matched (dim.3) for all
	  //   ABruijn vertices (dim.2) in each component (dim.1)
	  vector<vector<list<TwoValues<int> > > > abVerticesAll(numElemsInSets);
	  for(unsigned int i=0;i<abVertices.size();i++) abVerticesAll[i].resize(0);

	    vector<vector<short> > abCounts(numElemsInSets);   // Records info on number of vertices and edges per ABruijn graph
	    for(unsigned int i=0; i<abCounts.size(); i++) { abCounts[i].resize(2); abCounts[i][0]=0; abCounts[i][1]=0; }

	//  for(unsigned int cIdx=20; cIdx<21; cIdx++) {
	  unsigned int prevNumSpecs;
	    for(unsigned int cIdx=0; cIdx<components.size(); cIdx++) {
	    STDOUT << "Processing component "<<cIdx<<endl;

	    prevNumSpecs = 0;
	    while(prevNumSpecs!=components.sets[cIdx].size()) {  // Iterates ABruijn/sequencing/split until all
	                                                       //   remaining spectra match the best ABruijn path
	      prevNumSpecs = components.sets[cIdx].size();

	      STDOUT << "  - Spectrum indices: "; for(list<int>::iterator iter=components.sets[cIdx].begin();iter!=components.sets[cIdx].end();iter++) STDOUT<<*iter<<" "; STDOUT<<endl; STDOUT.flush();
	      STDOUT << "  - Component defined by "<< components.cAlignsPA[cIdx].size() << " pairs...\n"; STDOUT.flush();
	      vector<vector<TwoValues<int> > > matches;

	      //
	      // Choose consensus orientations for pairwise alignments - a spectrum can only be
	      //   used as-is or reversed, not both. Also determines set of matched peaks for
	      //   the consensus orientations.
	      //
	      vector<float> modPos;
	  //    SplitPairs3(specSet, cAlignsASP[cIdx], components.cAlignsPA[cIdx], peakTol, maxAAjump, penalty_sameVert, penalty_ptm, matches, matchesPA, specFlipped, modPos, false, labelsP);
	      SplitPairs(specSet, components.cAlignsPA[cIdx], peakTol, pmTol, maxAAjump, maxModMass, penalty_sameVert, penalty_ptm, matches, specFlipped, modPos, minMatchedPeaks, minEdgesToComponent, false, NULL, labelsP);
	      SplitASPPA_v2(components.cAlignsPA[cIdx],pmTol,tmpAlignsASP,tmpAlignsPA);

	/*      vector<bool> present(specSet.size()); for(unsigned int i=0; i<present.size();i++) present[i]=false;
	      for(unsigned int i=0;i<components.cAlignsPA[cIdx].size();i++)
	        if(matches[i].size()>0) {
	          if(!present[components.cAlignsPA[cIdx][i].spec1]) { cSpectra[cIdx].push_back(components.cAlignsPA[cIdx][i].spec1); present[components.cAlignsPA[cIdx][i].spec1]=true; }
	          if(!present[components.cAlignsPA[cIdx][i].spec2]) { cSpectra[cIdx].push_back(components.cAlignsPA[cIdx][i].spec2); present[components.cAlignsPA[cIdx][i].spec2]=true; }
	        }
	*/
	      //
	      // Maximize endpoint scores for the matched spectra
	      //
	      list<int>::iterator sicIter;
	      if(addEndpoints)
	        for(sicIter=cSpectra[cIdx].begin();sicIter!=cSpectra[cIdx].end();sicIter++)
	          specSet[*sicIter].maximizeZPMpeaks(peakTol,true);

	  #ifdef DBG_MASAB
	      // Add labels to the split pairs graph for graphviz output
	      MSGraph g;   g.build(tmpAlignsASP);   g.add(tmpAlignsPA);
	      g.vLabels.resize(specSet.size());
	      if(specSet.size()==specSet.size())
	        for(unsigned int i=0; i<specSet.size(); i++) {
	          if(specFlipped[i]) sprintf(sBuf,"v%d_R",i); else sprintf(sBuf,"v%d",i);
	          g.vLabels[i] = string((const char *)sBuf);
	        }
	      else for(unsigned int i=0; i<specSet.size(); i++) {
	          sprintf(sBuf,"v%d",i);    g.vLabels[2*i] = string((const char *)sBuf);
	          sprintf(sBuf,"v%d_R",i);  g.vLabels[2*i+1] = string((const char *)sBuf);
	         }
	      sprintf(sBuf,"split_pairs_graph_%d.txt",cIdx);  g.output_graphviz(sBuf);
	  #endif

	      //
	      // Build spectrum graphs for the spectra in this component
	      //
	      if(graphType>0)
	        for(unsigned int i=0; i<components.cAlignsPA[cIdx].size(); i++) {
	          if(spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].numVerts()==0) {
	            if(graphType==1) spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].ConnectConsecutive(specSet[components.cAlignsPA[cIdx][i].spec1]);
	            else spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].ConnectJumps(specSet[components.cAlignsPA[cIdx][i].spec1],jumps2,peakTol);
	//sprintf(sBuf,"spectrum_graph_%d.txt",components.cAlignsPA[cIdx][i].spec1);
	//spectrumGraphs[components.cAlignsPA[cIdx][i].spec1].output_graphviz(sBuf);
	          }
	          if(spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].numVerts()==0) {
	            if(graphType==1) spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].ConnectConsecutive(specSet[components.cAlignsPA[cIdx][i].spec2]);
	            else spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].ConnectJumps(specSet[components.cAlignsPA[cIdx][i].spec2],jumps2,peakTol);
	//sprintf(sBuf,"spectrum_graph_%d.txt",components.cAlignsPA[cIdx][i].spec2);
	//spectrumGraphs[components.cAlignsPA[cIdx][i].spec2].output_graphviz(sBuf);
	          }
	        }

	      //
	      // Build A-Bruijn graph
	      //
	      VertexSet vSet(specSet,2048);
	      vector<bool> usedSpectra(specSet.size());  for(unsigned int i=0; i<specSet.size(); i++) usedSpectra[i]=false;
	      for(unsigned int i=0; i<matches.size(); i++)
	        if(matches[i].size()>0) {
	          if(graphType==0)  // adds edges for all consecutive matched spectrum peaks
	            vSet.addGlues(components.cAlignsPA[cIdx][i].spec1,components.cAlignsPA[cIdx][i].spec2,matches[i]);
	          else  // add edges later
	            vSet.addGlues(components.cAlignsPA[cIdx][i].spec1,components.cAlignsPA[cIdx][i].spec2,matches[i],&spectrumGraphs);
	          usedSpectra[components.cAlignsPA[cIdx][i].spec1] = true;
	          usedSpectra[components.cAlignsPA[cIdx][i].spec2] = true;
	        }
	      list<int> includedSpecs;     includedSpecs.clear();
	      for(unsigned int i=0; i<specSet.size(); i++)
	        if(usedSpectra[i]) includedSpecs.push_back(i);

	      int specCount=0; for(unsigned int i=0;i<usedSpectra.size();i++) if(usedSpectra[i]) specCount++;
	      if(specCount==0) { STDOUT<<"  - component is empty, ABruijn graph not built\n"; STDOUT.flush(); continue; }
	      STDOUT<<"  - ABruijn graph built on "<<specCount<<" spectra.\n";

	      if(specCount>10000) {
	        STDOUT << "Spectrum count per component exceeds maximum of 7500 spectra! Skipping component...\n";
	        continue;
	      }

	      //
	      // Find/count/list/split composite vertices
	      //
	      int compositeVertexCount=0;  list<int> compositeSet;
	      for(unsigned int i=0; i<vSet.vertices.size(); i++)
	        if(vSet.vertices[i].size()>0 and vSet.vertices[i].compositeVertex)
	          { compositeVertexCount++; compositeSet.push_front(i);
	STDOUT<<">>>>> composite vertex "<<i<<": ";
	for(list<TwoValues<int> >::iterator iter=vSet.vertices[i].specPeaks.begin(); iter!=vSet.vertices[i].specPeaks.end(); iter++)
	  STDOUT<<"("<<(*iter)[0]<<","<<(*iter)[1]<<")";
	STDOUT<<endl;
	          }
	      STDOUT<<"  - Abruijn graph contains "<<compositeVertexCount<<" composite vertices: ";
	      list<int>::iterator iter=compositeSet.begin();
	      for(; iter!=compositeSet.end(); iter++) STDOUT<<(*iter)<<" ";
	      if(compositeVertexCount>0) {
	        STDOUT<<"-> splitting...";STDOUT.flush();
	        vSet.splitComposite(spectrumGraphs,peakTol,&usedSpectra);
	        STDOUT<<"done.\n";
	      } else STDOUT<<endl;
	      STDOUT.flush();

	      //
	      // Add spectrum graph edges to ABruijn graph (may connect non-aligned peaks, e.g. if
	      //   missing in all but one of the aligned spectra)
	      //
	      if(graphType==1) { vSet.addEdges(spectrumGraphs,&usedSpectra);  vSet.consolidatePaths(); }
	      if(graphType==2) { vSet.addEdges(spectrumGraphs,&usedSpectra); }
	// Changed 2010/09/14     if(not addEndpoints) vSet.addEndpointEdges(tmpAlignsASP,matches,modPos,jumps2,peakTol);
	      if(addEndpoints) vSet.addEndpointEdges(tmpAlignsASP,matches,modPos,jumps2,peakTol);

	      //
	      //  ABruijn path-finishing procedures should go here.
	      //

	  #ifdef DBG_MASAB
	      vSet.outputGraph(labelsP);
	  #endif

	      //
	      // Create ABruijn MSGraph from VertexSet for b-ion and y-ion endpoints; choose option with highest-scoring path
	      //
	      VertexSet copy(vSet), *vSetP;
	      MSGraph *abg, abgB, abgY, path;
	      vector<int> vSet_indexB, vSet_indexY, *vSet_indexP; // Correspondences between vertex indices in the ABruijn graph and simplified graph
	      vector<int> pathVertsIdxB, pathVertsIdxY, *pathVertsIdx;  // Indices of the vertices in the heaviest path
	      Spectrum specWithBep, specWithYep;
	      float scoreWithBep, scoreWithYep;

	      if(addEndpoints) vSet.removeEndpoints(false,peakTol);  // Remove y-mass endpoints added above
	      //
	      //  ABruijn endpoint-finishing procedures should go here.
	      //
	      vSet.buildGraph(abgB,jumps,peakTol,vSet_indexB,edgeScoreType,labelsP);
	//char filename[2048];   sprintf(filename,"component_%d.txt",cIdx+1);
	//abgB.output_graphviz(filename);
	      if(noSequencing) {
	        pathVertsIdxB.resize(0);
	        if(wholeABFN) vSet.getMatchedPeaks(pathVertsIdxB,abVerticesAll[cIdx]);
	        continue;
	      }
	      scoreWithBep = abgB.heaviestPath(path,false,&specWithBep,&pathVertsIdxB);

	      if(addEndpoints) copy.removeEndpoints(true,peakTol);  // Remove b-mass endpoints added above
	      //
	      //  ABruijn endpoint-finishing procedures should go here.
	      //
	      copy.buildGraph(abgY,jumps,peakTol,vSet_indexY,edgeScoreType,labelsP);
	//sprintf(filename,"component_%d_y.txt",cIdx+1);
	//abgY.output_graphviz(filename);
	      scoreWithYep = abgY.heaviestPath(path,false,&specWithYep,&pathVertsIdxY);

	//      scoreWithBep = max(scoreWithBep,scoreWithYep)+1;  // Force B endpoints

	      if(scoreWithBep>=scoreWithYep) { abg=&abgB; pathSpectra.consensus[cIdx]=specWithBep; vSet_indexP = &vSet_indexB; vSetP = &vSet; pathVertsIdx=&pathVertsIdxB; STDOUT<<"  - Selected B endpoints.\n";}
	        else { abg=&abgY; pathSpectra.consensus[cIdx]=specWithYep; vSet_indexP = &vSet_indexY; vSetP = &copy; pathVertsIdx=&pathVertsIdxY; STDOUT<<"  - Selected Y endpoints.\n";}
	      for(unsigned int vIdx=0; vIdx<pathVertsIdx->size(); vIdx++) (*pathVertsIdx)[vIdx] = (*vSet_indexP)[(*pathVertsIdx)[vIdx]];  // Convert simplified graph vertex indices to ABruijn vertex indices.

	  #ifdef DBG_MASAB
	  //    sprintf(sBuf,"graph_ma_%d.txt",cIdx);
	  //    vSetP->output_graphviz_ma(sBuf, *pathVertsIdx);
	  #endif
	      vSetP->getMatchedPeaks(*pathVertsIdx,abVertices[cIdx]);
	      pathSpectra.set_endpoints(cIdx,specSet,abVertices[cIdx],peakTol,pmTol);

	      if(wholeABFN) {
	        pathVertsIdx->resize(0);
	        vSetP->getMatchedPeaks(*pathVertsIdx,abVerticesAll[cIdx]);
	      }

	      abg->info_heaviestPath(cStats[cIdx]);   abCounts[cIdx][0]=abg->numVertices();   abCounts[cIdx][1]=abg->numEdges();
	  //    STDOUT<<"  - Heaviest path stats: ["; for(unsigned int i=0;i<cStats[cIdx].size();i++) {STDOUT<<cStats[cIdx][i]; if(i<cStats[cIdx].size()-1) STDOUT<<", "; } STDOUT<<"]\n";

	  #ifdef DBG_MASAB
	      sprintf(sBuf,"graph_%d.txt",cIdx); abg->output_graphviz(sBuf);
	  #endif

	      //
	      // Find spectra with not enough ABruijn vertices on the heaviest path (if any)
	      //   and remove them from the current component. Unused spectra may define new (leftover)
	      //   components if connected by at least one edge (ie, at least two spectra in a
	      //   connected component).
	      //
	      list<int> usedSpecs;     usedSpecs.clear();
	      list<int> unusedSpecs;   unusedSpecs.clear();
	      vector<short> numPathPeaks(specSet.size());   for(unsigned int i=0; i<specSet.size(); i++) numPathPeaks[i]=0;
	      for(unsigned int i=0; i<abVertices[cIdx].size(); i++) {
	//STDOUT<<"AB vertex "<<i<<": ";
	        list<TwoValues<int> >::iterator vNext; // Used to remove remaining composite vertices
	        for(list<TwoValues<int> >::iterator vIter=abVertices[cIdx][i].begin(); vIter!=abVertices[cIdx][i].end(); ) {
	//STDOUT<<"("<<(*vIter)[0]<<","<<(*vIter)[1]<<")";
	          vNext=vIter;  if(vNext!=abVertices[cIdx][i].end()) vNext++;
	          if(vNext!=abVertices[cIdx][i].end() and (*vIter)[0]==(*vNext)[0]) {
	            STDOUT<<"  - ERROR: inconsistent contig vertex containing ("<<(*vIter)[0]<<","<<(*vIter)[1]<<") and ("<<(*vNext)[0]<<","<<(*vNext)[1]<<")"; STDOUT.flush();
	            vIter=abVertices[cIdx][i].erase(vIter);  // Remove inconsistent spectrum/peaks
	            vIter=abVertices[cIdx][i].erase(vIter);
	          } else {
	            if(++numPathPeaks[(*vIter)[0]]>=pathMinPeaks) usedSpecs.push_back((*vIter)[0]);
	            vIter++;
	          }
	        }
	//STDOUT<<endl;
	      }
	      usedSpecs.sort();     usedSpecs.unique();

	      // Prune ABruijn vertices: remove peaks from spectra without enough matches to the consensus path
	      unsigned int maxPeaksPerVertex=0;
	      for(unsigned int i=0; i<abVertices[cIdx].size(); i++) {
	        for(list<TwoValues<int> >::iterator vIter=abVertices[cIdx][i].begin(); vIter!=abVertices[cIdx][i].end(); )
	          if(numPathPeaks[(*vIter)[0]]>=pathMinPeaks) vIter++;
	          else { unusedSpecs.push_back((*vIter)[0]); vIter=abVertices[cIdx][i].erase(vIter); }
	        if(abVertices[cIdx][i].size()>maxPeaksPerVertex) maxPeaksPerVertex=abVertices[cIdx][i].size();
	        if(abVertices[cIdx][i].size()>usedSpecs.size()) {
	          STDOUT<<"  - ERROR: inconsistent contig vertex: "; STDOUT.flush();
	          for(list<TwoValues<int> >::iterator vIter=abVertices[cIdx][i].begin(); vIter!=abVertices[cIdx][i].end(); vIter++)
	            { STDOUT<<"("<<(*vIter)[0]<<","<<(*vIter)[1]<<")"; STDOUT.flush(); }
	          STDOUT<<"\n";
	        }
	      }
	      unusedSpecs.sort();   unusedSpecs.unique();

	      // Output complete set of vertices
	      if(wholeABFN) {
	        for(unsigned int i=0; i<abVerticesAll[cIdx].size(); i++)
	          for(list<TwoValues<int> >::iterator vIter=abVerticesAll[cIdx][i].begin(); vIter!=abVerticesAll[cIdx][i].end(); )
	            if(numPathPeaks[(*vIter)[0]]>=pathMinPeaks) vIter++; else vIter=abVerticesAll[cIdx][i].erase(vIter);
	      }

	      if(usedSpecs.size()>=pathMinSpecs)  {
	        if (usedSpecs.size()<components.sets[cIdx].size()-1) {  // Whenever there are at least 2 unused spectra
	          STDOUT<<"  - Keeping "<<usedSpecs.size()<<" spectra; number of components: "<<components.size()<<" -> ";
	          components.splitSet(cIdx,usedSpecs);
	          STDOUT<<components.size()<<"\n"; STDOUT.flush();
	        } else if (usedSpecs.size()==components.sets[cIdx].size()-1) components.removeElements(cIdx,unusedSpecs);
	        if(maxPeaksPerVertex>usedSpecs.size()) {
	          STDOUT<<"  - ERROR: inconsistent contig containing a vertex with "<<maxPeaksPerVertex<<" spectrum peaks but only "<<usedSpecs.size()<<" used spectra (contig deleted)\n";
	          pathSpectra.consensus[cIdx].resize(0);
	          abVertices[cIdx].resize(0);
	        }
	      } else {
	        if (includedSpecs.size()>0) {  // Whenever there is at least 1 included spectrum
	          STDOUT<<"  - Keeping "<<components.sets[cIdx].size()-includedSpecs.size()<<" spectra; number of components: "<<components.size()<<" -> ";
	          components.splitSet(cIdx,includedSpecs);
	          STDOUT<<components.size()<<"\n"; STDOUT.flush();
	        }
	        STDOUT<<"  - Only "<<usedSpecs.size()<<" spectra with at least "<<pathMinPeaks<<" peaks in the consensus path - component deleted.\n";
	        pathSpectra.consensus[cIdx].resize(0);  // Remove de novo sequences for poor contigs (not enough spectra with enough matched peaks)
	      }

	      STDOUT.flush();   STDOUT.flush();

	      if(components.cAlignsPA[cIdx].size()==0 and components.cAlignsASP[cIdx].size()==0) { // No pairs left
	        prevNumSpecs = components.sets[cIdx].size();
	        STDOUT<<"  --> No pairs left for second iteration - keeping ABruijn path\n";
	      }
	    } // while (prevNumSpecs!=components.sets[cIdx].size())
	  }
	  // Resize down to the final number of resulting connected components
	  abCounts.resize(components.size());
	  cStats.resize(components.size());
	  pathSpectra.resize(components.size());
	  cSpectra.resize(components.size());
	  abVertices.resize(components.size());
	  abVerticesAll.resize(components.size());
	
	/*
	SetMerger COMPONENTS
	vector<vector<float> > COMPONENT_STATS;
	vector<list<int> > COMPONENT_SPECTRA;
	abinfo_t COMPONENT_INFO;
	abinfo_t COMPLETE_ABRUIJN;
	Clusters PATH_SPECTRA_AS_CLUSTER;
	vector<vector<short> > ABCOUNTS;
	*/
	m_parameters.COMPONENTS = components;
	m_parameters.COMPONENT_STATS = cStats;
	m_parameters.COMPONENT_SPECTRA = cSpectra;
	m_parameters.PATH_SPECTRA_AS_CLUSTER = pathSpectra;
	m_parameters.ABCOUNTS = abCounts;
	m_parameters.OUTPUT_SPECS = pathSpectra.consensus;
	
	Merge_abinfo_v1_0(specSet, components.sets, specFlipped, abVertices, m_parameters.COMPONENT_INFO);
	if (wholeABFN) Merge_abinfo_v1_0(specSet, components.sets, specFlipped, abVerticesAll, m_parameters.COMPLETE_ABRUIJN);
}
