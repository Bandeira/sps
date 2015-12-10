/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 10/27/09
 */

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
#include <sys/stat.h>

#include "SpsMergeHandler.h"
#include "parameter_names.h"
#include "../utils.h"
#include "../inputParams.h"

using namespace std;

int main(int argc, char *argv[], char **envp)
{

  set < string > flags;
  for (unsigned int i = 0; i < NUM_FLAGS; i++)
  {
    string flag(FLAGS[i]);
    flags.insert(flag);
  }
  set < string > stages;
  for (unsigned int i = 0; i < NUM_STAGES; i++)
  {
    string stage(STAGES[i]);
    stages.insert(stage);
  }
  map < string, string > parsedArgs;
  parseArguments(argv, argc, flags, parsedArgs);

  string break_flag(BREAK_FLAG);
  string param_flag(PARAM_FILE_FLAG);
  string stage_flag(STAGE_FLAG);
  string stats_flag(STATS_FLAG);

  const char* paramFile = (parsedArgs.count(param_flag) == 0)
      ? "spsmerge.params" : parsedArgs[param_flag].c_str();

  if (parsedArgs.count(stage_flag) == 0)
  {
    cerr << "ERROR: Missing " << STAGE_FLAG << " option\n";
    cout << endl << SPSMERGE_USAGE;
    return -1;
  }
  const char* stage = parsedArgs[stage_flag].c_str();

  AAJumps jumps(0);

  bool go_to_next_stage = false;
  struct stat st;

  //Contig-contig alignment stage
  if (strcmp(stage, START) == 0 || strcmp(stage, CONTIG_ALIGN) == 0)
  {

    InputParams params;

    if (!params.readParams(paramFile))
    {
      cerr << "ERROR: Cannot open parameters file " << paramFile << endl;
      cout << endl << SPSMERGE_USAGE;
      return -1;
    }

    ContigContigAlignParams contig_align_params;

    if (parsedArgs.count(stats_flag) > 0)
    {

      CombineContigsParams contig_stats_params;

      printf("\nPROCESSING %s FLAG FOR STAGE %s ...\n\n",
             STATS_FLAG,
             CONTIG_ALIGN);
      fflush( stdout);
      vector<const char *> paramStrs(9);

      paramStrs[0] = INPUT_CONTIG_SPECTRA;
      paramStrs[1] = MIN_MATCHED_PEAKS_CONTIG_ALIGN;
      paramStrs[2] = TOLERANCE_PEAK;
      paramStrs[3] = TOLERANCE_PM;
      paramStrs[4] = RESOLUTION;
      paramStrs[5] = MIN_CONTIG_OVERLAP_SCORE;
      paramStrs[6] = INPUT_FASTA;
      paramStrs[7] = INPUT_MATCHED_PEAKS_IDX;
      paramStrs[8] = INPUT_MATCHED_PROTS;

      if (!params.confirmParams(paramStrs))
      {
        cerr << "ERROR: Parameters file " << paramFile
            << " is incomplete for stage " << stage << " with " << stats_flag
            << " option. One of the following is missing: " << paramStrs[0];
        for (int i = 1; i < paramStrs.size(); i++)
          cerr << ", " << paramStrs[i];
        cerr << endl;
        return -1;
      }

      contig_stats_params.min_matched_peaks_contig_align
          = params.getValueInt(MIN_MATCHED_PEAKS_CONTIG_ALIGN);

      vector<bool> paramCnfrm;
      vector<TwoValues<string> > additonalParams;
      TwoValues<string> nextP(INPUT_CONTIG_ALIGNMENT, DEFAULT_CONTIG_ALIGNMENT);
      additonalParams.push_back(nextP);
      nextP[0] = OUTPUT_CONTIG_ALIGN_STATS;
      nextP[1] = DEFAULT_CONTIG_ALIGN_STATS;
      additonalParams.push_back(nextP);

      bool one_added = params.addParams(additonalParams, paramCnfrm);
      if (one_added)
      {
        bool added = mkdir_if_not_exist(CONTIG_ASSEMBLY_DIR);

        if (!added && stat(CONTIG_ASSEMBLY_DIR, &st) != 0)
        {
          cerr << "ERROR: Failed to create directory " << CONTIG_ASSEMBLY_DIR
              << endl;
          return -1;
        }

        FILE* params_out = fopen(DEFAULT_CONTIG_ALIGN_STATS_PARAMS, "w");
        if (params_out != NULL)
        {
          printf("Saving parameters file for contig-contig alignment to %s. To re-run this stage do:\n\nspsmerge --p %s --s %s %s\n\n",
                 DEFAULT_CONTIG_ALIGN_STATS_PARAMS,
                 DEFAULT_CONTIG_ALIGN_STATS_PARAMS,
                 CONTIG_ALIGN,
                 STATS_FLAG);
          params.dumpParams(params_out);
          fclose(params_out);
        }
      }

      if (!SpsMergeHandler::GetContigShiftStats(params.getValue(INPUT_CONTIG_SPECTRA),
                                                params.getValue(INPUT_CONTIG_ABINFO),
                                                params.getValue(INPUT_CONTIG_ALIGNMENT),
                                                params.getValue(INPUT_MATCHED_PROTS),
                                                params.getValue(INPUT_MATCHED_PEAKS_IDX),
                                                params.getValue(INPUT_FASTA),
                                                params.getValue(OUTPUT_CONTIG_ALIGN_STATS),
                                                contig_stats_params))
        return -1;

      return 0;
    }

    printf("\nBEGINNING STAGE %s ...\n\n", CONTIG_ALIGN);
    fflush( stdout);

    vector<const char *> paramStrs(6);
    paramStrs[0] = INPUT_CONTIG_SPECTRA;
    paramStrs[1] = MIN_MATCHED_PEAKS_CONTIG_ALIGN;
    paramStrs[2] = TOLERANCE_PEAK;
    paramStrs[3] = TOLERANCE_PM;
    paramStrs[4] = RESOLUTION;
    paramStrs[5] = MIN_CONTIG_OVERLAP_SCORE;

    if (!params.confirmParams(paramStrs))
    {
      cerr << "ERROR: Parameters file " << paramFile
          << " is incomplete for stage " << stage
          << ". One of the following is missing: " << paramStrs[0];
      for (int i = 1; i < paramStrs.size(); i++)
        cerr << ", " << paramStrs[i];
      cerr << endl;
      return -1;
    }

    contig_align_params.minScore
        = (float)params.getValueDouble(MIN_CONTIG_OVERLAP_SCORE);
    contig_align_params.minNumMatchedPeaks
        = params.getValueInt(MIN_MATCHED_PEAKS_CONTIG_ALIGN);
    contig_align_params.allowedJumps = &jumps;

    unsigned int num_threads = (params.paramPresent(NUM_THREADS))
        ? (unsigned int)params.getValueInt(NUM_THREADS) : 1;

    bool added = params.addParam(OUTPUT_CONTIG_ALIGNMENT,
                                 DEFAULT_CONTIG_ALIGNMENT);
    if (added)
    {
      added = mkdir_if_not_exist(CONTIG_ASSEMBLY_DIR);

      if (!added && stat(CONTIG_ASSEMBLY_DIR, &st) != 0)
      {
        cerr << "ERROR: Failed to create directory " << CONTIG_ASSEMBLY_DIR
            << endl;
        return -1;
      }

      FILE* params_out = fopen(DEFAULT_CONTIG_ALIGN_PARAMS, "w");
      if (params_out != NULL)
      {
        printf("Saving parameters file for contig-contig alignment to %s. To re-run this stage do:\n\nspsmerge --p %s --s %s %s\n\n",
               DEFAULT_CONTIG_ALIGN_PARAMS,
               DEFAULT_CONTIG_ALIGN_PARAMS,
               CONTIG_ALIGN,
               BREAK_FLAG);
        params.dumpParams(params_out);
        fclose(params_out);
      }
    }

    if (!SpsMergeHandler::AlignContigs(params.getValue(INPUT_CONTIG_SPECTRA),
                                       params.getValue(OUTPUT_CONTIG_ALIGNMENT),
                                       contig_align_params,
                                       num_threads))
      return -1;

    printf("\nCOMPLETED STAGE %s ...\n\n", CONTIG_ALIGN);
    fflush(stdout);

    if (parsedArgs.count(break_flag) > 0)
      go_to_next_stage = false;
    else
      go_to_next_stage = true;
  }

  //Contig-contig merging stage
  if (go_to_next_stage || strcmp(stage, CONTIG_MERGE) == 0)
  {

    printf("\nBEGINNING STAGE %s ...\n\n", CONTIG_MERGE);
    fflush( stdout);

    InputParams params;

    if (!params.readParams(paramFile))
    {
      cerr << "ERROR: Cannot open parameters file " << paramFile << endl;
      cout << endl << SPSMERGE_USAGE;
      return -1;
    }

    CombineContigsParams contig_merge_params;

    vector<const char *> paramStrs(9);
    paramStrs[0] = INPUT_CONTIG_SPECTRA;
    paramStrs[1] = INPUT_STAR_SPECTRA;
    paramStrs[2] = MIN_MATCHED_PEAKS_CONTIG_ALIGN;
    paramStrs[3] = TOLERANCE_PEAK;
    paramStrs[4] = TOLERANCE_PM;
    paramStrs[5] = RESOLUTION;
    paramStrs[6] = MIN_CONTIG_OVERLAP_SCORE;
    paramStrs[7] = INPUT_CONTIG_ABINFO;
    paramStrs[8] = MIN_COMPONENT_SIZE;

    if (!params.confirmParams(paramStrs))
    {
      cerr << "ERROR: Parameters file " << paramFile
          << " is incomplete for stage " << stage
          << ". One of the following is missing: " << paramStrs[0];
      for (int i = 1; i < paramStrs.size(); i++)
        cerr << ", " << paramStrs[i];
      cerr << endl;
      return -1;
    }

    contig_merge_params.contig_overlap_score
        = (float)params.getValueDouble(MIN_CONTIG_OVERLAP_SCORE);
    contig_merge_params.min_matched_peaks_contig_align
        = params.getValueInt(MIN_MATCHED_PEAKS_CONTIG_ALIGN);
    contig_merge_params.peak_tol = (float)params.getValueDouble(TOLERANCE_PEAK);
    contig_merge_params.parent_mass_tol
        = (float)params.getValueDouble(TOLERANCE_PM);
    contig_merge_params.resolution = (float)params.getValueDouble(RESOLUTION);
    contig_merge_params.min_component_size = params.getValueInt(MIN_COMPONENT_SIZE);

    vector<bool> paramCnfrm(5);
    vector<TwoValues<string> > additonalParams;
    TwoValues<string> nextP(INPUT_CONTIG_ALIGNMENT, DEFAULT_CONTIG_ALIGNMENT);
    additonalParams.push_back(nextP);
    nextP[0] = OUTPUT_CONTIG_ASSEMBLY_SPECTRA;
    nextP[1] = DEFAULT_CONTIG_ASSEMBLY_SPECTRA;
    additonalParams.push_back(nextP);
    nextP[0] = OUTPUT_CONTIG_ASSEMBLY_ABINFO;
    nextP[1] = DEFAULT_CONTIG_ASSEMBLY_ABINFO;
    additonalParams.push_back(nextP);
    nextP[0] = OUTPUT_CONTIG_ASSEMBLY_COMPONENTS;
    nextP[1] = DEFAULT_CONTIG_ASSEMBLY_COMPONENTS;
    additonalParams.push_back(nextP);
    nextP[0] = OUTPUT_CONTIG_ASSEMBLY_REVERSED;
    nextP[1] = DEFAULT_CONTIG_ASSEMBLY_REVERSED;
    additonalParams.push_back(nextP);

    bool one_added = params.addParams(additonalParams, paramCnfrm);
    if (one_added)
    {
      bool added = mkdir_if_not_exist(CONTIG_ASSEMBLY_DIR);

      if (!added && stat(CONTIG_ASSEMBLY_DIR, &st) != 0)
      {
        cerr << "ERROR: Failed to create directory " << CONTIG_ASSEMBLY_DIR
            << endl;
        return -1;
      }

      FILE* params_out = fopen(DEFAULT_CONTIG_MERGE_PARAMS, "w");
      if (params_out != NULL)
      {
        printf("Saving parameters file for contig-contig merging to %s. To re-run this stage do:\n\nspsmerge --p %s --s %s %s\n\n",
               DEFAULT_CONTIG_MERGE_PARAMS,
               DEFAULT_CONTIG_MERGE_PARAMS,
               CONTIG_MERGE,
               BREAK_FLAG);
        params.dumpParams(params_out);
        fclose(params_out);
      }
    }

    vector<const char*> matchmaParams(4);
    matchmaParams[0] = INPUT_FASTA;
    matchmaParams[1] = INPUT_MATCHED_PEAKS_IDX;
    matchmaParams[2] = INPUT_MATCHED_PROTS;
    matchmaParams[3] = CONTAMINANT_PROTEINS;

    if (params.confirmParams(matchmaParams))
    {
      set<int> protein_idx_ignore;
      list < string > strIdxs;

      if (!splitText(params.getValue(CONTAMINANT_PROTEINS), strIdxs, ":"))
        return -1;

      for (list<string>::iterator strIt = strIdxs.begin(); strIt
          != strIdxs.end(); strIt++)
        protein_idx_ignore.insert(atoi((*strIt).c_str()));

      contig_merge_params.protein_idx_ignore = &protein_idx_ignore;

      if (!SpsMergeHandler::MergeContigsContigsMatchma(params.getValue(INPUT_CONTIG_SPECTRA),
                                                       params.getValue(INPUT_STAR_SPECTRA),
                                                       params.getValue(INPUT_CONTIG_ABINFO),
                                                       params.getValue(INPUT_CONTIG_ALIGNMENT),
                                                       params.getValue(INPUT_MATCHED_PROTS),
                                                       params.getValue(INPUT_MATCHED_PEAKS_IDX),
                                                       params.getValue(INPUT_FASTA),
                                                       params.getValue(OUTPUT_CONTIG_ASSEMBLY_SPECTRA),
                                                       params.getValue(OUTPUT_CONTIG_ASSEMBLY_ABINFO),
                                                       params.getValue(OUTPUT_CONTIG_ASSEMBLY_COMPONENTS),
                                                       params.getValue(OUTPUT_CONTIG_ASSEMBLY_REVERSED),
                                                       contig_merge_params))
        return -1;
    }
    else
    {
      contig_merge_params.protein_idx_ignore = NULL;

      if (!SpsMergeHandler::MergeContigsContigs(params.getValue(INPUT_CONTIG_SPECTRA),
                                                params.getValue(INPUT_STAR_SPECTRA),
                                                params.getValue(INPUT_CONTIG_ABINFO),
                                                params.getValue(INPUT_CONTIG_ALIGNMENT),
                                                params.getValue(OUTPUT_CONTIG_ASSEMBLY_SPECTRA),
                                                params.getValue(OUTPUT_CONTIG_ASSEMBLY_ABINFO),
                                                params.getValue(OUTPUT_CONTIG_ASSEMBLY_COMPONENTS),
                                                params.getValue(OUTPUT_CONTIG_ASSEMBLY_REVERSED),
                                                contig_merge_params))
        return -1;
    }

    printf("\nCOMPLETED STAGE %s ...\n\n", CONTIG_MERGE);
    fflush(stdout);

    if (parsedArgs.count(break_flag) > 0)
      go_to_next_stage = false;
    else
      go_to_next_stage = true;
  }

  //Matepair-contig alignment stage
  if (strcmp(stage, MP_ALIGN) == 0)
  {

    printf("\nBEGINNING STAGE %s ...\n\n", MP_ALIGN);
    fflush( stdout);

    InputParams params;

    if (!params.readParams(paramFile))
    {
      cerr << "ERROR: Cannot open parameters file " << paramFile << endl;
      cout << endl << SPSMERGE_USAGE;
      return -1;
    }

    MatepairContigAlignParams mp_align_params;

    vector<const char *> paramStrs(7);
    paramStrs[0] = INPUT_MATEPAIR_SPECTRA;
    paramStrs[1] = TOLERANCE_PEAK;
    paramStrs[2] = TOLERANCE_PM;
    paramStrs[3] = RESOLUTION;
    paramStrs[4] = MINIMUM_MATEPAIR_CHARGE;
    paramStrs[5] = MIN_RATIO_MATEPAIR_ALIGN;
    paramStrs[6] = MIN_MATCHED_PEAKS_MATEPAIR_ALIGN;

    if (!params.confirmParams(paramStrs))
    {
      cerr << "ERROR: Parameters file " << paramFile
          << " is incomplete for stage " << stage
          << ". One of the following is missing: " << paramStrs[0];
      for (int i = 1; i < paramStrs.size(); i++)
        cerr << ", " << paramStrs[i];
      cerr << endl;
      return -1;
    }

    mp_align_params.minRatio
        = (float)params.getValueDouble(MIN_RATIO_MATEPAIR_ALIGN);
    mp_align_params.minNumMatchedPeaks
        = params.getValueInt(MIN_MATCHED_PEAKS_MATEPAIR_ALIGN);
    mp_align_params.precCharge
        = (short)params.getValueInt(MINIMUM_MATEPAIR_CHARGE);
    mp_align_params.checkConsPeaks = false;
    mp_align_params.minAbsShift = 0.0;
    mp_align_params.allowedJumps = &jumps;

    unsigned int num_threads = (params.paramPresent(NUM_THREADS))
        ? (unsigned int)params.getValueInt(NUM_THREADS) : 1;

    vector<bool> paramCnfrm(5);
    vector<TwoValues<string> > additonalParams;
    TwoValues<string> nextP(OUTPUT_MP_ALIGNMENT, DEFAULT_MP_ALIGNMENT);
    additonalParams.push_back(nextP);
    nextP[0] = INPUT_MP_ASSEMBLY_SPECTRA;
    nextP[1] = DEFAULT_CONTIG_ASSEMBLY_SPECTRA;
    additonalParams.push_back(nextP);

    bool one_added = params.addParams(additonalParams, paramCnfrm);
    if (one_added)
    {
      bool added = mkdir_if_not_exist(MP_ASSEMBLY_DIR);

      if (!added && stat(MP_ASSEMBLY_DIR, &st) != 0)
      {
        cerr << "ERROR: Failed to create directory " << CONTIG_ASSEMBLY_DIR
            << endl;
        return -1;
      }

      FILE* params_out = fopen(DEFAULT_MP_ALIGN_PARAMS, "w");
      if (params_out != NULL)
      {
        printf("Saving parameters file for contig-matepair alignment to %s. To re-run this stage do:\n\nspsmerge --p %s --s %s %s\n\n",
               DEFAULT_MP_ALIGN_PARAMS,
               DEFAULT_MP_ALIGN_PARAMS,
               MP_ALIGN,
               BREAK_FLAG);
        params.dumpParams(params_out);
        fclose(params_out);
      }
    }

    if (!SpsMergeHandler::AlignMatepairs(params.getValue(INPUT_MP_ASSEMBLY_SPECTRA),
                                         params.getValue(INPUT_MATEPAIR_SPECTRA),
                                         params.getValue(OUTPUT_MP_ALIGNMENT),
                                         mp_align_params,
                                         num_threads))
      return -1;

    printf("\nCOMPLETED STAGE %s ...\n\n", MP_ALIGN);
    fflush(stdout);

    if (parsedArgs.count(break_flag) > 0)
      go_to_next_stage = false;
    else
      go_to_next_stage = true;
  }

  return 0;
}
