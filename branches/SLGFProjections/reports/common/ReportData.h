////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_DATA_H__
#define __REPORT_DATA_H__
////////////////////////////////////////////////////////////////////////////////
#include "aminoacid.h"
#include "Defines.h"

////////////////////////////////////////////////////////////////////////////////
#define DEFAULT_ROOT_DIRECTORY      "."
#define DEFAULT_TABLES_DIR          "./report"
#define DEFAULT_OUTPUT_DIR          "./report"
#define FILE_CLUSTERMS              "spectra/pklbin_files.txt"
#define FILE_SCAN_FILES             "spectra/bin_files.txt"
#define FILE_CONSENSUS_SPECTRA      "spectra/specs_ms.pklbin"
#define FILE_PROTEINS               "protid.fasta"
#define FILE_CONTIG_INDICES         "spectra/contigs_indices.bin"
#define FILE_HOMGLUE_REF_MIDX       "homology/homglue_ref_midx.pklbin"
#define FILE_HOMGLUE_REF_MP         "homology/homglue_ref_mp.bin"
#define FILE_CONTIGS_MIDX           "homology/contigs_midx.pklbin"
#define FILE_CONTIGS_MP             "homology/contigs_mp.bin"
#define FILE_SPS_SEQS               "assembly/sps_seqs.pklbin"
#define FILE_ABRUIJN                "assembly/component_info.bin"
#define FILE_STARS_SPECTRA          "spectra/stars.pklbin"
#define FILE_CONTIG_SPECTRA         "spectra/contigs.pklbin"
#define FILE_HOMGLUE_MATCHES_MP     "homology/homglue_matches_mp.bin"
#define FILE_HOMGLUE_MATCHES_MIDX   "homology/homglue_matches_midx.pklbin"
#define FILE_HOMGLUE_MATCHES        "homology/homglue_matches.pklbin"
#define FILE_INPUT_SPECTRA_LIST     "spectra/input_index.txt"
#define FILE_CONTIG_NAMES           "homology/ref_sps_names.txt"
#define FILE_INPUT_MAPPING          "spectra/input_mapping.bin"

#define DEFAULT_ANNOTATION_FILE_LOCATION "."

#define DEFAULT_MASS_SHIFT            specnets::AAJumps::massHion
#define DEFAULT_PEAK_MASS_TOLERANCE   0.45
#define DEFAULT_PARENT_MASS_TOLERANCE 1.0
#define DEFAULT_RESOLUTION            0.01
#define DEFAULT_CELLS_PER_LINE        20

#define DEFAULT_TOOL                  1

#define DEFAULT_JOB_NAME          "MyJob"
#define DEFAULT_USER_NAME         "Undefined"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
class ReportData {

 public:

  string  exeDir;

  // Directories and data file names
  string  projectDir;
  string  tablesDir;
  string  outDir;
  string  server;

  // directories for report. Used for relocation.
  string  targetProjectDir;

  int     cellsPerLine;
  int     displayLevel;

  int     cpu;

  int     tool;

  bool    noClusters;

  bool    allowRealign;

  string  annotationModelDir;
  string  annotationModel;
  float   massShift;
  string  annotationModelPrm;
  float   massShiftPrm;
  string  annotationModelDirPrm;

  float   peakMassTol;
  float   parentMassTol;
  float   resolution;

  int     tableID;
  int     filterField;
  string  filterData;
  int     updateField;
  string  updateData;

  int     sortColumnIdx;
  int     sortDirection;

  int     startRow;
  int     rows;


  ReportData()
  {
    projectDir          = DEFAULT_ROOT_DIRECTORY;
    tablesDir           = DEFAULT_TABLES_DIR;
    outDir              = DEFAULT_OUTPUT_DIR;
    cellsPerLine        = DEFAULT_CELLS_PER_LINE;

    annotationModelDir  = DEFAULT_ANNOTATION_FILE_LOCATION;
    annotationModel     = DEFAULT_ANNOTATION_MODEL;

    massShift           = DEFAULT_MASS_SHIFT;
    massShiftPrm        = DEFAULT_MASS_SHIFT;

    peakMassTol         = DEFAULT_PEAK_MASS_TOLERANCE;
    parentMassTol       = DEFAULT_PARENT_MASS_TOLERANCE;
    resolution          = DEFAULT_RESOLUTION;

    tableID             = -1;
    filterField         = -1;
    updateField         = -1;
    displayLevel        = 100;

    sortColumnIdx       = -1;
    sortDirection       = 0;

    startRow            = 0;
    rows                = -1;

    cpu                 = 1;

    tool                = DEFAULT_TOOL;

    noClusters          = false;

    allowRealign        = true;

  }

};

////////////////////////////////////////////////////////////////////////////////
class ReportGeneratorData : public ReportData {

 public:


  string  job;
  string  user;

  string  filenameClusterMS;
  string  filenameScanFiles;
  string  filenameConsensusSpectra;
  string  filenameProteins;
  string  filenameContigIndices;
  string  filenameHomglueRefMidx;
  string  filenameHomglueRefMp;
  string  filenameContigsMidx;
  string  filenameContigsMp;
  string  filenameSpsSeqs;
  string  filenameAbruijn;
  string  filenameStarSpectra;
  string  filenameContigSpectra;
  string  filenameHomglueMatchesMp;
  string  filenameHomglueMatchesMidx;
  string  filenameHomglueMatches;
  string  filenameInputSpectraList;
  string  filenameContigNames;
  string  filenameInputMapping;


  // constructor
  ReportGeneratorData() :
      filenameClusterMS(FILE_CLUSTERMS),
      filenameScanFiles(FILE_SCAN_FILES),
      filenameConsensusSpectra(FILE_CONSENSUS_SPECTRA),
      filenameProteins(FILE_PROTEINS),
      filenameContigIndices(FILE_CONTIG_INDICES),
      filenameHomglueRefMidx(FILE_HOMGLUE_REF_MIDX),
      filenameHomglueRefMp(FILE_HOMGLUE_REF_MP),
      filenameContigsMidx(FILE_CONTIGS_MIDX),
      filenameContigsMp(FILE_CONTIGS_MP),
      filenameSpsSeqs(FILE_SPS_SEQS),
      filenameAbruijn(FILE_ABRUIJN),
      filenameStarSpectra(FILE_STARS_SPECTRA),
      filenameContigSpectra(FILE_CONTIG_SPECTRA),
      filenameHomglueMatchesMp(FILE_HOMGLUE_MATCHES_MP),
      filenameHomglueMatchesMidx(FILE_HOMGLUE_MATCHES_MIDX),
      filenameHomglueMatches(FILE_HOMGLUE_MATCHES),
      filenameInputSpectraList(FILE_INPUT_SPECTRA_LIST),
      filenameContigNames(FILE_CONTIG_NAMES),
      filenameInputMapping(FILE_INPUT_MAPPING),

      job(DEFAULT_JOB_NAME),
      user(DEFAULT_USER_NAME)

  {};

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
