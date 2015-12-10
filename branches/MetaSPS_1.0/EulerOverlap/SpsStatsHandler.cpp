/*
 * SpsStatsHandler.cpp
 *
 *  Created on: Oct 1, 2010
 *      Author: aguthals
 */

#include "SpsStatsHandler.h"

bool SpsStatsHandler::OutputCummulativeStatistics(const char* in_contigs_file,
                                                  const char* in_abinfo_file,
                                                  const char* in_stars_file,
                                                  const char* in_matched_peaks_file,
                                                  const char* in_matched_prots_file,
                                                  const char* in_fasta_file,
                                                  const char* in_inspect_file,
                                                  const char* in_iontypes_file,
                                                  const char* out_stats_file,
                                                  const char* out_contig_stats_root,
                                                  set<int>& target_proteins,
                                                  float peak_tol,
                                                  float pm_tol)
{

  SpecSet* contigs = new SpecSet;
  SpecSet* stars = new SpecSet;
  abinfo_t* in_abinfo = new abinfo_t;
  SpecSet* overlaps = new SpecSet;
  vector < vector<int> > *prot_match = new vector<vector<int> > ;
  DB_fasta* fasta = new DB_fasta;
  MS2ScoringModel* model = new MS2ScoringModel;
  InspectResultsSet* in_results = new InspectResultsSet;
  unsigned int textLength = strlen(in_iontypes_file);
  char* textCpy = (char*)malloc(textLength + 1);
  vector < string > target_prots;
  vector < string > peptides;
  SpecnetsStatistics stats;
  string file_str = out_stats_file;
  string root_str = out_contig_stats_root;

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
  peptides.resize(stars->size(), "");

  printf("Loading abinfo from %s ...\n", in_abinfo_file);
  fflush(stdout);
  if (!Load_abinfo(in_abinfo_file, *in_abinfo))
  {
    cerr << "ERROR: Failed to load " << in_abinfo_file << "\n";
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

  printf("Loading ion types model from %s ...\n", in_iontypes_file);
  fflush(stdout);

  if (strncpy(textCpy, in_iontypes_file, textLength) == NULL)
  {
    cerr << "ERROR COPYING TEXT\n";
    free(textCpy);
    goto exit;
  }

  textCpy[textLength] = '\0';

  if (!model->LoadModel(textCpy))
  {
    cerr << "ERROR: Failed to load " << textCpy << "\n";
    free(textCpy);
    goto exit;
  }
  free(textCpy);

  printf("Loading Inspect results from %s ...\n", in_inspect_file);
  fflush(stdout);
  if (!in_results->loadInspectResultsFile(in_inspect_file))
  {
    cerr << "ERROR: Failed to load " << in_inspect_file << "\n";
    goto exit;
  }

  FilterFastaProteins(*fasta, target_proteins, target_prots);

  FilterInspectIds(*in_results, peptides);

  delete fasta;
  delete in_results;

  stats.reset(*contigs,
              *in_abinfo,
              *stars,
              *overlaps,
              *prot_match,
              target_prots,
              peptides,
              *model,
              peak_tol,
              pm_tol);

  cout << endl;
  if (!stats.outputCummulativeCSVReport(file_str))
    return false;

  if (!stats.outputAllContigsCSVReport(root_str))
    return false;

  return true;

  exit:

  delete contigs;
  delete stars;
  delete in_abinfo;
  delete overlaps;
  delete prot_match;
  delete fasta;
  delete model;
  delete in_results;
  return false;
}

void SpsStatsHandler::FilterFastaProteins(DB_fasta& fasta,
                                          set<int>& target_proteins,
                                          vector<string>& put_proteins)
{
  put_proteins.resize(fasta.size());
  for (int i = 0; i < put_proteins.size(); i++)
  {
    if (target_proteins.count(i) > 0 || target_proteins.size() == 0)
      put_proteins[i] = fasta.getSequence(i);
    else
      put_proteins[i] = "";
  }
}

void SpsStatsHandler::FilterInspectIds(InspectResultsSet& results, vector<
    string>& put_peptides)
{

  for (int i = 0; i < put_peptides.size(); i++)
  {
    if (results.results.count(i) == 0)
      put_peptides[i] = "";
    else
      results.inspectToSpecNets(results.results[i].Annotation, put_peptides[i]);
  }
}
