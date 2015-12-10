/*
 * SpsStatsHandler.h
 *
 *  Created on: Oct 1, 2010
 *      Author: aguthals
 */

#ifndef SPSSTATSHANDLER_H_
#define SPSSTATSHANDLER_H_

#include <cstring>
#include <vector>
#include <set>
#include <string>

#include "../spectrum.h"
#include "../inspect_parse.h"
#include "../specnets_statistics.h"
#include "../db_fasta.h"

/*

 #include "../inspect_parse.h"
 #include "../spectrum.cpp"
 #include "../utils.h"
 */
using namespace std;

class SpsStatsHandler
{
public:

  static bool OutputCummulativeStatistics(const char* in_contigs_file,
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
                                          float pm_tol);

  static void FilterFastaProteins(DB_fasta& fasta,
                                  set<int>& target_proteins,
                                  vector<string>& put_proteins);

  static void FilterInspectIds(InspectResultsSet& results,
                               vector<string>& put_peptides);
};

#endif /* SPSSTATSHANDLER_H_ */
