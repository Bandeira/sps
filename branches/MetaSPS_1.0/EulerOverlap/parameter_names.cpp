/*
 * parameter_names.cpp
 *
 *  Created on: Oct 4, 2010
 *      Author: aguthals
 */

#include "parameter_names.h"

//SPSMERGE PARAMETERS

//Command-line flags and arguments
const char* PARAM_FILE_FLAG = "--p";
const char* STAGE_FLAG = "--s";
const char* BREAK_FLAG = "--b";
const char* STATS_FLAG = "--stats";

const char* FLAGS[] = { PARAM_FILE_FLAG, STAGE_FLAG, BREAK_FLAG, STATS_FLAG };

const char* START = "start";
const char* RESUME = "resume";
const char* CONTIG_ALIGN = "contig-align";
const char* CONTIG_MERGE = "contig-merge";
const char* MP_ALIGN = "mp-align";
const char* MP_MERGE = "mp-merge";

const char* STAGES[] = { START, RESUME, CONTIG_ALIGN, CONTIG_MERGE, MP_ALIGN,
                         MP_MERGE };

//Usage
const char
    * SPSMERGE_USAGE =
        "spsmerge --p [param-file] --s [stage]\n\n    --p : input parameters \
        file (required) (see trunk/EulerOverlap/spsmerge.params for template)\
        , default value is \"spsmerge.params\" if not passed\n    --s : stage\
         to run (required), options are \"contig-align\" (alias=\"start\"), \"\
        contig-merge\", \"mp-align\", or \"mp-merge\" (alias=\"resume\")\n    \
        --b : break execution after finishing specified stage\n    --stats : \
        just output statistics on stage's results\n\n";

//Global input parameters
const char* INPUT_CONTIG_SPECTRA = "INPUT_CONTIG_SPECTRA (*.pklbin)";
const char* INPUT_MATEPAIR_SPECTRA = "INPUT_MATEPAIR_SPECTRA (*.pklbin)";
const char* INPUT_CONTIG_ABINFO = "INPUT_CONTIG_ABINFO";

const char* INPUT_FASTA = "INPUT_FASTA";
const char* INPUT_MATCHED_PEAKS_IDX = "INPUT_MATCHED_PEAKS (*midx.pklbin)";
const char* INPUT_MATCHED_PROTS = "INPUT_MATCHED_PROTS (*mp.bin)";
const char* CONTAMINANT_PROTEINS = "CONTAMINANT_PROTEINS (*:*)";

const char* TOLERANCE_PEAK = "TOLERANCE_PEAK";
const char* TOLERANCE_PM = "TOLERANCE_PM";
const char* RESOLUTION = "RESOLUTION";
const char* MIN_CONTIG_OVERLAP_SCORE = "MIN_CONTIG_OVERLAP_SCORE";
const char* MIN_MATCHED_PEAKS_CONTIG_ALIGN = "MIN_MATCHED_PEAKS_CONTIG_ALIGN";
const char* MINIMUM_MATEPAIR_CHARGE = "MINIMUM_MATEPAIR_CHARGE";
const char* MIN_RATIO_MATEPAIR_ALIGN = "MIN_RATIO_MATEPAIR_ALIGN";
const char* MIN_RATIO_MATEPAIR_CONTIG_ALIGN = "MIN_RATIO_MATEPAIR_CONTIG_ALIGN";
const char* MIN_MATCHED_PEAKS_MATEPAIR_ALIGN =
    "MIN_MATCHED_PEAKS_MATEPAIR_ALIGN";

const char* NUM_THREADS = "NUM_THREADS";

const char* GRID_NUMNODES = "GRID_NUMNODES";
const char* GRID_EXE_DIR = "GRID_EXE_DIR";

const char* MIN_NUM_MATEPAIRS_MERGE = "MIN_NUM_MATEPAIRS_MERGE";
const char* MIN_EDGES_TO_COMPONENT = "MIN_EDGES_TO_COMPONENT";
const char* MIN_COMPONENT_SIZE = "MIN_COMPONENT_SIZE";

//Local Parameters
const char* CONTIG_ASSEMBLY_DIR = "contig-assembly";
const char* MP_ASSEMBLY_DIR = "mp-assembly";
const char* MP_GRID_DIR = "mp-assembly/grid";

const char* OUTPUT_CONTIG_ALIGNMENT = "OUTPUT_CONTIG_ALIGNMENT (results_OCC)";
const char* INPUT_CONTIG_ALIGNMENT = "INPUT_CONTIG_ALIGNMENT (results_OCC)";
const char* DEFAULT_CONTIG_ALIGNMENT = "contig-assembly/aligns_contigs.txt";

const char* OUTPUT_MP_ALIGNMENT = "OUTPUT_MP_ALIGNMENT (results_OCC)";
const char* INPUT_MP_ALIGNMENT = "INPUT_MP_ALIGNMENT (results_OCC)";
const char* DEFAULT_MP_ALIGNMENT = "mp-assembly/aligns_mate-pairs.txt";

const char* DEFAULT_CONTIG_ALIGN_PARAMS = "contig-assembly/spsalign.params";
const char* DEFAULT_CONTIG_ALIGN_STATS_PARAMS =
    "contig-assembly/spsalign_stats.params";
const char* DEFAULT_CONTIG_MERGE_PARAMS = "contig-assembly/spsmerge.params";

const char* DEFAULT_MP_ALIGN_PARAMS = "mp-assembly/spsalign.params";
const char* DEFAULT_MP_MERGE_PARAMS = "mp-assembly/spsmerge.params";

const char* MATEPAIR_IDXS = "MATEPAIR_IDXS";

const char* OUTPUT_CONTIG_ALIGN_STATS = "OUTPUT_CONTIG_ALIGN_STATS (*.csv)";
const char* DEFAULT_CONTIG_ALIGN_STATS =
    "contig-assembly/contig_align_stats.csv";

const char* OUTPUT_MP_ALIGN_STATS = "OUTPUT_MP_ALIGN_STATS (*.csv)";
const char* DEFAULT_MP_ALIGN_STATS = "mp-assembly/mp_align_stats.csv";

const char* OUTPUT_CONTIG_ASSEMBLY_COMPONENTS =
    "OUTPUT_CONTIG_ASSEMBLY_COMPONENTS (*.csv)";
const char* DEFAULT_CONTIG_ASSEMBLY_COMPONENTS =
    "contig-assembly/final_components.csv";

const char* OUTPUT_CONTIG_ASSEMBLY_SPECTRA =
    "OUTPUT_CONTIG_ASSEMBLY_SPECTRA (*.pklbin)";
const char* INPUT_MP_ASSEMBLY_SPECTRA = "INPUT_MP_ASSEMBLY_SPECTRA (*.pklbin)";
const char* DEFAULT_CONTIG_ASSEMBLY_SPECTRA =
    "contig-assembly/meta_sps_seqs.pklbin";

const char* OUTPUT_CONTIG_ASSEMBLY_REVERSED =
    "OUTPUT_CONTIG_ASSEMBLY_REVERSED (*.bin)";
const char* DEFAULT_CONTIG_ASSEMBLY_REVERSED =
    "contig-assembly/sps_seqs_reversed.bin";

const char* OUTPUT_CONTIG_ASSEMBLY_ABINFO = "OUTPUT_CONTIG_ASSEMBLY_ABINFO";
const char* INPUT_MP_ASSEMBLY_ABINFO = "INPUT_MP_ASSEMBLY_ABINFO (*.pklbin)";
const char* DEFAULT_CONTIG_ASSEMBLY_ABINFO =
    "contig-assembly/component_info.bin";

const char* OUTPUT_MP_ASSEMBLY_COMPONENTS =
    "OUTPUT_MP_ASSEMBLY_COMPONENTS (*.csv)";
const char* DEFAULT_MP_ASSEMBLY_COMPONENTS = "mp-assembly/final_components.csv";

const char* OUTPUT_MP_ASSEMBLY_SPECTRA =
    "OUTPUT_MP_ASSEMBLY_SPECTRA (*.pklbin)";
const char* DEFAULT_MP_ASSEMBLY_SPECTRA = "mp-assembly/meta_sps_seqs.pklbin";

const char* OUTPUT_MP_ASSEMBLY_REVERSED = "OUTPUT_MP_ASSEMBLY_REVERSED (*.bin)";
const char* DEFAULT_MP_ASSEMBLY_REVERSED = "mp-assembly/sps_seqs_reversed.bin";

const char* OUTPUT_MP_ASSEMBLY_ABINFO = "OUTPUT_MP_ASSEMBLY_ABINFO";
const char* DEFAULT_MP_ASSEMBLY_ABINFO = "mp-assembly/component_info.bin";

//SPSSTATS PARAMETERS

const char* SPSSTATS_USAGE = "spsstats [parameters file]\n\n";

const char* INPUT_STAR_SPECTRA = "INPUT_STAR_SPECTRA (*.pklbin)";
const char* INPUT_ION_TYPES = "INPUT_ION_TYPES";
const char* INPUT_INSPECT_RESULTS = "INPUT_INSPECT_RESULTS";

const char* TARGET_PROTEINS = "TARGET_PROTEINS (*:*)";

const char* OUTPUT_STATS_FILE = "OUTPUT_STATS_FILE";
const char* OUTPUT_CONTIG_STATS_FILEROOT = "OUTPUT_CONTIG_STATS_FILEROOT";
