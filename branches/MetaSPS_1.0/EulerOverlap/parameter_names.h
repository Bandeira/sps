/*
 * parameter_names.h
 *
 *  Created on: Oct 4, 2010
 *      Author: aguthals
 */

#ifndef PARAMETER_NAMES_H_
#define PARAMETER_NAMES_H_

//SPSMERGE PARAMETERS

//Command-line flags and arguments
extern const char* PARAM_FILE_FLAG;
extern const char* STAGE_FLAG;
extern const char* BREAK_FLAG;
extern const char* STATS_FLAG;

const unsigned int NUM_FLAGS = 4;
extern const char* FLAGS[NUM_FLAGS];

extern const char* START;
extern const char* RESUME;
extern const char* CONTIG_ALIGN;
extern const char* CONTIG_MERGE;
extern const char* MP_ALIGN;
extern const char* MP_MERGE;

const unsigned int NUM_STAGES = 6;
extern const char* STAGES[NUM_STAGES];

//Usage
extern const char* SPSMERGE_USAGE;

//Global input parameters
extern const char* INPUT_CONTIG_SPECTRA;
extern const char* INPUT_MATEPAIR_SPECTRA;
extern const char* INPUT_CONTIG_ABINFO;

extern const char* INPUT_FASTA;
extern const char* INPUT_MATCHED_PEAKS_IDX;
extern const char* INPUT_MATCHED_PROTS;
extern const char* CONTAMINANT_PROTEINS;

extern const char* TOLERANCE_PEAK;
extern const char* TOLERANCE_PM;
extern const char* RESOLUTION;
extern const char* MIN_CONTIG_OVERLAP_SCORE;
extern const char* MIN_MATCHED_PEAKS_CONTIG_ALIGN;
extern const char* MINIMUM_MATEPAIR_CHARGE;
extern const char* MIN_RATIO_MATEPAIR_ALIGN;
extern const char* MIN_RATIO_MATEPAIR_CONTIG_ALIGN;
extern const char* MIN_MATCHED_PEAKS_MATEPAIR_ALIGN;

extern const char* NUM_THREADS;

extern const char* GRID_NUMNODES;
extern const char* GRID_EXE_DIR;

extern const char* MIN_NUM_MATEPAIRS_MERGE;
extern const char* MIN_EDGES_TO_COMPONENT;
extern const char* MIN_COMPONENT_SIZE;

//Local Parameters
extern const char* CONTIG_ASSEMBLY_DIR;
extern const char* MP_ASSEMBLY_DIR;
extern const char* MP_GRID_DIR;

extern const char* OUTPUT_CONTIG_ALIGNMENT;
extern const char* INPUT_CONTIG_ALIGNMENT;
extern const char* DEFAULT_CONTIG_ALIGNMENT;

extern const char* OUTPUT_MP_ALIGNMENT;
extern const char* INPUT_MP_ALIGNMENT;
extern const char* DEFAULT_MP_ALIGNMENT;

extern const char* DEFAULT_CONTIG_ALIGN_PARAMS;
extern const char* DEFAULT_CONTIG_ALIGN_STATS_PARAMS;
extern const char* DEFAULT_CONTIG_MERGE_PARAMS;

extern const char* DEFAULT_MP_ALIGN_PARAMS;
extern const char* DEFAULT_MP_MERGE_PARAMS;

extern const char* MATEPAIR_IDXS;

extern const char* OUTPUT_CONTIG_ALIGN_STATS;
extern const char* DEFAULT_CONTIG_ALIGN_STATS;

extern const char* OUTPUT_MP_ALIGN_STATS;
extern const char* DEFAULT_MP_ALIGN_STATS;

extern const char* OUTPUT_CONTIG_ASSEMBLY_COMPONENTS;
extern const char* DEFAULT_CONTIG_ASSEMBLY_COMPONENTS;

extern const char* OUTPUT_CONTIG_ASSEMBLY_SPECTRA;
extern const char* INPUT_MP_ASSEMBLY_SPECTRA;
extern const char* DEFAULT_CONTIG_ASSEMBLY_SPECTRA;

extern const char* OUTPUT_CONTIG_ASSEMBLY_REVERSED;
extern const char* DEFAULT_CONTIG_ASSEMBLY_REVERSED;

extern const char* OUTPUT_CONTIG_ASSEMBLY_ABINFO;
extern const char* INPUT_MP_ASSEMBLY_ABINFO;
extern const char* DEFAULT_CONTIG_ASSEMBLY_ABINFO;

extern const char* OUTPUT_MP_ASSEMBLY_COMPONENTS;
extern const char* DEFAULT_MP_ASSEMBLY_COMPONENTS;

extern const char* OUTPUT_MP_ASSEMBLY_SPECTRA;
extern const char* DEFAULT_MP_ASSEMBLY_SPECTRA;

extern const char* OUTPUT_MP_ASSEMBLY_REVERSED;
extern const char* DEFAULT_MP_ASSEMBLY_REVERSED;

extern const char* OUTPUT_MP_ASSEMBLY_ABINFO;
extern const char* DEFAULT_MP_ASSEMBLY_ABINFO;

//SPSSTATS PARAMETERS

extern const char* SPSSTATS_USAGE;

extern const char* INPUT_STAR_SPECTRA;
extern const char* INPUT_ION_TYPES;
extern const char* INPUT_INSPECT_RESULTS;

extern const char* TARGET_PROTEINS;

extern const char* OUTPUT_STATS_FILE;
extern const char* OUTPUT_CONTIG_STATS_FILEROOT;

#endif /* PARAMETER_NAMES_H_ */
