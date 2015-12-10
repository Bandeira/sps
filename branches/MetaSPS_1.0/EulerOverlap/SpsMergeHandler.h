#ifndef SPSMERGEHANDLER_H
/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 09/02/10
 */
#define SPSMERGEHANDLER_H

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>
#include <list>
#include <vector>
#include <iomanip>

#include "CombineContigs.h"
#include "../spectrum.h"
#include "../aminoacid.h"
#include "../twovalues.h"
#include "../db_fasta.h"
#include "../utils.h"
#include "../batch.h"

class SpsMergeHandler
{
public:

  static bool AlignContigs(const char* contigs_file,
                           const char* results_file,
                           ContigContigAlignParams& align_params,
                           unsigned int num_threads);

  static bool MergeContigsContigs(const char* in_contigs_file,
                                  const char* in_stars_file,
                                  const char* in_abinfo_file,
                                  const char* in_contig_aligns_file,
                                  const char* out_contigs_file,
                                  const char* out_abinfo_file,
                                  const char* out_components_file,
                                  const char* out_reversed_file,
                                  CombineContigsParams& merging_params);

  static bool MergeContigsContigsMatchma(const char* in_contigs_file,
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
                                         CombineContigsParams& merging_params);

  static bool GetContigShiftStats(const char* in_contigs_file,
                                  const char* in_abinfo_file,
                                  const char* in_contig_aligns_file,
                                  const char* in_matched_prots_file,
                                  const char* in_matched_peaks_file,
                                  const char* in_fasta_file,
                                  const char* out_stats_file,
                                  CombineContigsParams& merging_params);

  static bool AlignMatepairs(const char* contigs_file,
                             const char* connectors_file,
                             const char* results_file,
                             MatepairContigAlignParams& align_params,
                             unsigned int num_threads);

};

#endif
