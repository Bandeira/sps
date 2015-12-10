#include "SpsMergeHandler.h"

bool SpsMergeHandler::AlignContigs(const char* contigs_file,
                                   const char* results_file,
                                   ContigContigAlignParams& align_params,
                                   unsigned int num_threads)
{

  SpecSet* contigs = new SpecSet;
  ContigContigAlignParams params;

  printf("Loading contigs from %s ...\n", contigs_file);
  fflush( stdout);
  if (!contigs->LoadSpecSet_pklbin(contigs_file))
  {
    cerr << "ERROR: Failed to load " << contigs_file << "\n";
    delete contigs;
    return false;
  }

  printf("\nAligning contigs...\n");
  fflush(stdout);

  copy(align_params, params);
  params.contigs = contigs;

  getContigContigAlignment(&params, num_threads);

  printf("Saving contig alignments to %s ...\n", results_file);
  fflush(stdout);
  if (!Save_resultsOCC(results_file, params.results))
  {
    cerr << "ERROR: Failed to save to " << results_file << "\n";
    return false;
  }

  return true;
}

bool SpsMergeHandler::MergeContigsContigs(const char* in_contigs_file,
                                          const char* in_stars_file,
                                          const char* in_abinfo_file,
                                          const char* in_contig_aligns_file,
                                          const char* out_contigs_file,
                                          const char* out_abinfo_file,
                                          const char* out_components_file,
                                          const char* out_reversed_file,
                                          CombineContigsParams& merging_params)
{

  SpecSet* contigs = new SpecSet;
  SpecSet* stars = new SpecSet;
  abinfo_t* in_abinfo = new abinfo_t;
  list<Results_OCC>* contig_alignments = new list<Results_OCC> ;
  CombineContigsParams params;
  CombineContigs comp;

  printf("Loading contigs from %s ...\n", in_contigs_file);
  fflush( stdout);
  if (!contigs->LoadSpecSet_pklbin(in_contigs_file))
  {
    cerr << "ERROR: Failed to load " << in_contigs_file << "\n";
    goto exit;
  }

  printf("Loading star spectra from %s ...\n", in_stars_file);
  fflush(stdout);
  if (!stars->LoadSpecSet_pklbin(in_stars_file))
  {
    cerr << "ERROR: Failed to load " << in_stars_file << "\n";
    goto exit;
  }

  printf("Loading abinfo from %s ...\n", in_abinfo_file);
  fflush(stdout);
  if (!Load_abinfo(in_abinfo_file, *in_abinfo))
  {
    cerr << "ERROR: Failed to load " << in_abinfo_file << "\n";
    goto exit;
  }

  printf("Loading alignments from %s ...\n", in_contig_aligns_file);
  fflush(stdout);
  if (!Load_resultsOCC(in_contig_aligns_file, *contig_alignments))
  {
    cerr << "ERROR: Failed to load " << in_contig_aligns_file << "\n";
    goto exit;
  }

  copy(merging_params, params);
  params.contigs = contigs;
  params.star_spectra = stars;
  params.in_abinfo = in_abinfo;
  params.contig_alignments = contig_alignments;
  params.contig_prot_idx = NULL;
  params.matepairs = NULL;

  printf("\nMerging contigs with contig-contig shifts ...\n");
  fflush(stdout);

  comp.construct(&params);
  comp.combineEdges(0);

  printf("Saving consensus spectra to %s ...\n", out_contigs_file);
  fflush(stdout);
  if (!(params.consensus).SaveSpecSet_pklbin(out_contigs_file))
  {
    cerr << "ERROR: Failed to save to " << out_contigs_file << "\n";
    return false;
  }

  printf("Saving consensus abinfo to %s ...\n", out_abinfo_file);
  fflush(stdout);
  if (!Save_abinfo_v1_0(out_abinfo_file, params.out_abinfo))
  {
    cerr << "ERROR: Failed to save to " << out_abinfo_file << "\n";
    return false;
  }

  printf("Saving contig components to %s ...\n", out_components_file);
  fflush(stdout);
  if (!comp.saveComponents(out_components_file))
  {
    cerr << "ERROR: Failed to save to " << out_components_file << "\n";
    return false;
  }

  printf("Saving reversed contig indicies to %s ...\n", out_reversed_file);
  fflush(stdout);
  if (!Save_binArray<int> (out_reversed_file, params.out_reversed))
  {
    cerr << "ERROR: Failed to save to " << out_reversed_file << "\n";
    return false;
  }

  return true;

  exit:

  delete contigs;
  delete stars;
  delete in_abinfo;
  delete contig_alignments;
  return false;
}

bool SpsMergeHandler::MergeContigsContigsMatchma(const char* in_contigs_file,
                                                 const char* in_stars_file,
                                                 const char* in_abinfo_file,
                                                 const char* in_contig_aligns_file,
                                                 const char* in_matched_prots_file,
                                                 const char* in_matched_peaks_file,
                                                 const char* in_fasta_file,
                                                 const char* out_contigs_file,
                                                 const char* out_abinfo_file,
                                                 const char* out_components_file,
                                                 const char* out_reversed_file,
                                                 CombineContigsParams& merging_params)
{

  SpecSet* contigs = new SpecSet;
  SpecSet* stars = new SpecSet;
  abinfo_t* in_abinfo = new abinfo_t;
  list<Results_OCC>* contig_alignments = new list<Results_OCC> ;
  SpecSet* overlaps = new SpecSet;
  vector < vector<int> > *prot_match = new vector<vector<int> > ;
  DB_fasta* fasta = new DB_fasta;
  CombineContigsParams params;
  CombineContigs comp;

  printf("Loading contigs from %s ...\n", in_contigs_file);
  fflush( stdout);
  if (!contigs->LoadSpecSet_pklbin(in_contigs_file))
  {
    cerr << "ERROR: Failed to load " << in_contigs_file << "\n";
    goto exit;
  }

  printf("Loading star spectra from %s ...\n", in_stars_file);
  fflush(stdout);
  if (!stars->LoadSpecSet_pklbin(in_stars_file))
  {
    cerr << "ERROR: Failed to load " << in_stars_file << "\n";
    goto exit;
  }

  printf("Loading abinfo from %s ...\n", in_abinfo_file);
  fflush(stdout);
  if (!Load_abinfo(in_abinfo_file, *in_abinfo))
  {
    cerr << "ERROR: Failed to load " << in_abinfo_file << "\n";
    goto exit;
  }

  printf("Loading alignments from %s ...\n", in_contig_aligns_file);
  fflush(stdout);
  if (!Load_resultsOCC(in_contig_aligns_file, *contig_alignments))
  {
    cerr << "ERROR: Failed to load " << in_contig_aligns_file << "\n";
    goto exit;
  }

  printf("Loading overlaps from %s ...\n", in_matched_peaks_file);
  fflush(stdout);
  if (!overlaps->LoadSpecSet_pklbin(in_matched_peaks_file))
  {
    cerr << "ERROR: Failed to load " << in_matched_peaks_file << "\n";
    goto exit;
  }

  printf("Loading matched protein idxs from %s ...\n", in_matched_prots_file);
  fflush(stdout);
  if (!Load_binArray(in_matched_prots_file, *prot_match))
  {
    cerr << "ERROR: Failed to load " << in_matched_prots_file << "\n";
    goto exit;
  }

  printf("Loading fasta proteins from %s ...\n", in_fasta_file);
  fflush(stdout);
  if (!fasta->Load(in_fasta_file))
  {
    cerr << "ERROR: Failed to load " << in_fasta_file << "\n";
    goto exit;
  }

  copy(merging_params, params);
  params.contigs = contigs;
  params.star_spectra = stars;
  params.in_abinfo = in_abinfo;
  params.contig_alignments = contig_alignments;
  params.contig_peak_idx = overlaps;
  params.contig_prot_idx = prot_match;
  params.proteins = fasta;
  params.matepairs = NULL;

  printf("\nMerging contigs with contig-contig shifts and matchma results ...\n");
  fflush(stdout);

  comp.construct(&params);
  comp.combineEdges(0);

  printf("Saving consensus spectra to %s ...\n", out_contigs_file);
  fflush(stdout);
  if (!(params.consensus).SaveSpecSet_pklbin(out_contigs_file))
  {
    cerr << "ERROR: Failed to save to " << out_contigs_file << "\n";
    return false;
  }

  printf("Saving consensus abinfo to %s ...\n", out_abinfo_file);
  fflush(stdout);
  if (!Save_abinfo_v1_0(out_abinfo_file, params.out_abinfo))
  {
    cerr << "ERROR: Failed to save to " << out_abinfo_file << "\n";
    return false;
  }

  printf("Saving contig components to %s ...\n", out_components_file);
  fflush(stdout);
  if (!comp.saveComponents(out_components_file))
  {
    cerr << "ERROR: Failed to save to " << out_components_file << "\n";
    return false;
  }

  printf("Saving reversed contig indicies to %s ...\n", out_reversed_file);
  fflush(stdout);
  if (!Save_binArray(out_reversed_file, params.out_reversed))
  {
    cerr << "ERROR: Failed to save to " << out_reversed_file << "\n";
    return false;
  }

  return true;

  exit:

  delete contigs;
  delete stars;
  delete in_abinfo;
  delete contig_alignments;
  delete overlaps;
  delete prot_match;
  delete fasta;
  return false;
}

bool SpsMergeHandler::GetContigShiftStats(const char* in_contigs_file,
                                          const char* in_abinfo_file,
                                          const char* in_contig_aligns_file,
                                          const char* in_matched_prots_file,
                                          const char* in_matched_peaks_file,
                                          const char* in_fasta_file,
                                          const char* out_stats_file,
                                          CombineContigsParams& merging_params)
{

  SpecSet* contigs = new SpecSet;
  abinfo_t* in_abinfo = new abinfo_t;
  list<Results_OCC>* contig_alignments = new list<Results_OCC> ;
  SpecSet* overlaps = new SpecSet;
  vector < vector<int> > *prot_match = new vector<vector<int> > ;
  DB_fasta* fasta = new DB_fasta;
  CombineContigsParams params;
  CombineContigs comp;
  set<int> protein_idx_ignore;

  printf("Loading contigs from %s ...\n", in_contigs_file);
  fflush( stdout);
  if (!contigs->LoadSpecSet_pklbin(in_contigs_file))
  {
    cerr << "ERROR: Failed to load " << in_contigs_file << "\n";
    goto exit;
  }

  printf("Loading abinfo from %s ...\n", in_abinfo_file);
  fflush(stdout);
  if (!Load_abinfo(in_abinfo_file, *in_abinfo))
  {
    cerr << "ERROR: Failed to load " << in_abinfo_file << "\n";
    goto exit;
  }

  printf("Loading alignments from %s ...\n", in_contig_aligns_file);
  fflush(stdout);
  if (!Load_resultsOCC(in_contig_aligns_file, *contig_alignments))
  {
    cerr << "ERROR: Failed to load " << in_contig_aligns_file << "\n";
    goto exit;
  }

  printf("Loading overlaps from %s ...\n", in_matched_peaks_file);
  fflush(stdout);
  if (!overlaps->LoadSpecSet_pklbin(in_matched_peaks_file))
  {
    cerr << "ERROR: Failed to load " << in_matched_peaks_file << "\n";
    goto exit;
  }

  printf("Loading matched protein idxs from %s ...\n", in_matched_prots_file);
  fflush(stdout);
  if (!Load_binArray(in_matched_prots_file, *prot_match))
  {
    cerr << "ERROR: Failed to load " << in_matched_prots_file << "\n";
    goto exit;
  }

  printf("Loading fasta proteins from %s ...\n", in_fasta_file);
  fflush(stdout);
  if (!fasta->Load(in_fasta_file))
  {
    cerr << "ERROR: Failed to load " << in_fasta_file << "\n";
    goto exit;
  }

  copy(merging_params, params);
  params.contigs = contigs;
  params.in_abinfo = in_abinfo;
  params.contig_alignments = contig_alignments;
  params.contig_peak_idx = overlaps;
  params.contig_prot_idx = prot_match;
  params.proteins = fasta;
  params.matepairs = NULL;

  params.protein_idx_ignore = &protein_idx_ignore;

  printf("\nGenerating contig-contig shift statistics ...\n");
  fflush(stdout);

  comp.construct(&params);
  if (comp.getContigShiftStats(out_stats_file,
                               merging_params.min_matched_peaks_contig_align,
                               merging_params.min_matched_peaks_contig_align
                                   + 3))
  {
    printf("\nWrote statistics to file %s ...\n", out_stats_file);
    fflush(stdout);
    return true;
  }
  return false;

  exit:

  delete contigs;
  delete in_abinfo;
  delete contig_alignments;
  delete overlaps;
  delete prot_match;
  delete fasta;
  return false;
}

bool SpsMergeHandler::AlignMatepairs(const char* contigs_file,
                                     const char* connectors_file,
                                     const char* results_file,
                                     MatepairContigAlignParams& align_params,
                                     unsigned int num_threads)
{

  SpecSet* contigs = new SpecSet;
  SpecSet* connectors = new SpecSet;
  MatepairContigAlignParams params;
  int totalpeaks;
  SpecSet* tempSpecs;
  SpecSet* tempPtr;

  printf("Loading contigs from %s ...\n", contigs_file);
  fflush( stdout);
  if (!contigs->LoadSpecSet_pklbin(contigs_file))
  {
    cerr << "ERROR: Failed to load " << contigs_file << "\n";
    goto exit;
  }

  printf("Loading matepairs from %s ...\n", connectors_file);
  fflush(stdout);
  if (!connectors->LoadSpecSet_pklbin(connectors_file))
  {
    cerr << "ERROR: Failed to load " << connectors_file << "\n";
    goto exit;
  }
  /*
   if (startIdx != 0 && endIdx != 0) {
   if (startIdx > endIdx) {
   cerr << "ERROR: Second slice index must be greater or equal to first in CONNECTOR_IDXS\n";
   delete contigs;
   delete connectors;
   return false;
   }

   if (startIdx < 0) {
   cerr << "ERROR: First slice index " << startIdx << " is < 0\n";
   delete contigs;
   delete connectors;
   return false;
   }
   if (endIdx > connectors->size()) {
   cerr << "ERROR: Second slice index " << endIdx << " is greater than size of connectors (" << connectors->size() << ")\n";
   delete contigs;
   delete connectors;
   return false;
   }

   tempSpecs = new SpecSet(endIdx - startIdx + 1);
   for (int i = startIdx; i < endIdx + 1; i++) {
   (*tempSpecs)[i - startIdx] = (*connectors)[i];
   }
   tempPtr = connectors;
   connectors = tempSpecs;
   delete tempPtr;
   }
   */
  totalpeaks = 0;
  for (int i = 0; i < connectors->size(); i++)
  {
    if ((*connectors)[i].parentCharge >= align_params.precCharge
        && (*connectors)[i].size() > 0)
    {
      totalpeaks += (*connectors)[i].size();
    }
  }

  //printf("\nAligning connectors %d - %d to all contigs (%d connectors, %d peaks)\n", startIdx, endIdx,  connectors->size(), totalpeaks ); fflush(stdout);

  printf("\nAligning matepairs to contigs ...\n");
  fflush(stdout);

  copy(align_params, params);
  params.connectors = connectors;
  params.contigs = contigs;
  //params.offset = startIdx;

  getMatepairContigAlignment(&params, num_threads);

  printf("\nSaving mate-pair alignments to %s ...\n", results_file);
  fflush(stdout);

  if (!Save_resultsOCC(results_file, params.results))
  {
    cerr << "ERROR: Failed to save to " << results_file << "\n";
    return false;
  }

  return true;

  exit:

  delete contigs;
  delete connectors;
  return false;
}
