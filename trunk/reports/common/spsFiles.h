////////////////////////////////////////////////////////////////////////////////
#ifndef __SPS_FILES_H__
#define __SPS_FILES_H__
////////////////////////////////////////////////////////////////////////////////
// includes section
#include <string>
#include <vector>

#include "spectrum.h"
#include "ClusterData.h"
#include "db_fasta.h"
#include "abruijn.h"
#include "aminoacid.h"

#include "Defines.h"
#include "ReportDefines.h"

#include "ReportData.h"
////////////////////////////////////////////////////////////////////////////////
// SPS Files
//
// This class contains all the files needed by sps
//
////////////////////////////////////////////////////////////////////////////////
// namespaces section
using namespace std;

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// Definitions section
typedef enum {
  SPS_FILE_FASTA,
  SPS_FILE_SPECSET,
  SPS_FILE_BIN_ARRAY,
  SPS_FILE_MSCLUSTER,
  SPS_FILE_STARS,
  SPS_FILE_SEQS,
  SPS_FILE_ABRUIJN,
  SPS_FILE_CONSENSUS_SPECTRA
} SpsFileID;
////////////////////////////////////////////////////////////////////////////////

class SpsFiles {

 public:

  SpsFiles();
  ~SpsFiles();


  void *getData(SpsFileID what, int index);

  //////////////////////////////////////////////////////////////////////////////
  // Spectrum data needed for several table fields.

  // input spectra filenames (pklbin)
  vector<string> m_inputSpectraPklbin;

  // Where to load input spectra data into (pklbin)     (1)
  // from <sps>/spectra/specs_ms_*.pklbin
  vector<specnets::SpecSet> m_specSet;


  // scan number files (text list)
  vector<string> m_scanNumberFiles;

  // Where to load input spectra scan numbers
  // from <sps>/spectra/specs_ms_*.pklbin
  vector<vector<vector<int> > > m_specScan;



  // Where to load consensus spectra data into (pklbin) (2)
  // from <sps>/spectra/specs_ms.pklbin
  specnets::SpecSet *m_consensusSpecSet;

  // where to store clusters information                (3)
  // from <sps>/clusterData.bin
  ClusterData *m_clusterData;

  // where to store the protein sequences               (16)
  specnets::DB_fasta  *m_fasta;

  // where to store the abruijn graph                   (8)
  specnets::abinfo_t  *m_abruijn;

  // where to store the star spectra                    (4)
  specnets::SpecSet *m_starSpectra;

  // where to load contig indices to                    (6)
  // from <sps>/spectra/contig_indices.bin
  vector<vector<int> >  *m_contigIndices;

  // where to load contig names into                    (17)
  vector<string>        *m_contigNames;

  // where to load input spectra file names into        (18)
  vector<string>        *m_input_index;

  // where to load input spectra mapping information to
  // from <sps>/spectra/input_mapping.bin
  vector<vector<int> >  *m_inputMapping;


  //////////////////////////////////////////////////////////////////////////////
  // Data needed to generate Reference sequences

  // where to store the spectra for contigs             (5)
  specnets::SpecSet *m_contigsSpectra;

  // where to load contig indices to                    (12)
  // from <sps>/homology/homglue_ref_midx.pklbin
  specnets::SpecSet *m_homglue_ref_midx;

  // where to load contig indices to                    (11)
  // from <sps>/homology/homglue_ref_mp.bin
  vector<vector<int> > *m_homglue_ref_mp;

  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate Homolog sequences

  // where to load contig indices to                    (7)
  // from <sps>/assembly/sps_seqs.pklbin
  specnets::SpecSet *m_sps_seqs;

  // where to load contig indices to                    (10)
  // from <sps>/homology/contigs_midx_all.pklbin
  specnets::SpecSet *m_contigs_midx;

  // where to load contig indices to                    (9)
  // from <sps>/homology/contigs_mp_all.bin
  vector<vector<int> > *m_contigs_mp;

  //////////////////////////////////////////////////////////////////////////////
  // data related CSPS contigs

  // where to store the CSPS spectra                    (15)
  specnets::SpecSet *m_homglueMatches;

  // where to load CSPS contig indices to               (14)
  // from <sps>/homology/homglue_matches_midx.pklbin
  specnets::SpecSet *m_homglue_matches_midx;

  // where to load CSPS contig indices to               (13)
  // from <sps>/homology/homglue_matches_mp.bin
  vector<vector<int> > *m_homglue_matches_mp;


  //////////////////////////////////////////////////////////////////////////////
  // Methods to load data files


  // General load data method
  virtual int loadData(const ReportGeneratorData &reportGeneratorData);


  // pklbin file load
  virtual int readPklbin(specnets::SpecSet *&pklbin, const string &fileName);
  // binArray file load
  virtual int readBinArray(vector<vector<int> >  *&binArray, const string &fileName);

  // filename composition helper method
  virtual string composeFileName(const string &projectDir, const string &fileName);

  // load the input spectra files list file
  virtual int loadInputSpectraList(const string &projectDir, const string &fileName);
  // (1) Load input spectra
  virtual int loadInputSpectraFiles(const string &projectDir, const string &fileName);
  // load scan numbers files
  virtual int loadScanFiles(const string &projectDir, const string &fileName);
  // (3) load MsCluster data
  virtual int loadMsClusterData(const string &projectDir, const string &fileName);
  // (2) Load cluster consensus spectra file
  virtual int loadConsensusSpectraFile(const string &projectDir, const string &fileName);

  // (6) Load contig indices
  virtual int loadContigIndices(const string &projectDir, const string &fileName);

  // (4) Load Star Spectra
  virtual int loadStarSpectra(const string &projectDir, const string &fileName);
  // (8) Load Abruijn graph
  virtual int loadAbruijn(const string &projectDir, const string &fileName);

  // (5) Load Spectra for contigs that matched a protein
  virtual int loadContigSpectra(const string &projectDir, const string &fileName);
  // (11) load <sps>/homology/homglue_ref_mp.bin
  virtual int loadHomglueRefMp(const string &projectDir, const string &fileName);
  // (12) load <sps>/homology/homglue_ref_midx.pklbin
  virtual int loadHomglueRefMidx(const string &projectDir, const string &fileName);

  // (7) load <sps>/assembly/sps_seqs.pklbin
  virtual int loadSpsSeqs(const string &projectDir, const string &fileName);
  // (9) load <sps>/homology/contigs_mp_all.bin
  virtual int loadContigsMpAll(const string &projectDir, const string &fileName);
  // (10) load <sps>/homology/contigs_midx_all.pklbin
  virtual int loadContigsMidxAll(const string &projectDir, const string &fileName);

  // (13) load <sps>/homology/homglue_matches_mp.pklbin
  virtual int loadCspsMatchesMp(const string &projectDir, const string &fileName);
  // (14) load <sps>/homology/homglue_matches_midx.bin
  virtual int loadCspsMatchesMidx(const string &projectDir, const string &fileName);
  // (15) load <sps>/homology/homglue_matches.pklbin
  virtual int loadCspsSpectra(const string &projectDir, const string &fileName);

  // (16) Load proteins file (FASTA format)
  virtual int loadProteinsFile(const string &projectDir, const string &fileName);

  // (17) Load contig names
  virtual int loadContigNames(const string &projectDir, const string &fileName);

  virtual int loadInputMapping(const string &projectDir, const string &fileName);




  // debug - dump abruijn graph
  void dump_abruijn(ostream &sout, bool web = false) {dump_abruijn(sout, m_abruijn, web);};
  void dump_abruijn(ostream &sout, specnets::abinfo_t *abruijn, bool web);
  void dump_clusterData(ostream &sout, ClusterData *cd, char *title, bool web);
  void dump_specset(ostream &sout, specnets::SpecSet *set, char *title, bool web);
  void dump_binArray(ostream &sout, vector<vector<int> > *a, char *title, bool web);

  void dump(ostream &sout, bool web = false) {
    dump_clusterData(sout, m_clusterData, "clusterData.bin", web);
    dump_binArray(sout, m_inputMapping, "input_mapping", web);
    dump_binArray(sout, m_contigIndices, "contigIndices", web);

    dump_binArray(sout, m_contigs_mp, "contigs_mp", web);
    dump_specset( sout, m_contigs_midx, "contigs_midx_all", web);

    dump_binArray(sout, m_homglue_ref_mp, "homglue_ref_mp", web);
    dump_specset( sout, m_homglue_ref_midx, "homglue_ref_midx", web);

    dump_specset(sout, m_homglueMatches, "homglueMatches", web);
    dump_binArray(sout, m_homglue_matches_mp, "homglue_matches_mp", web);
    dump_specset( sout, m_homglue_matches_midx, "homglue_matches_midx", web);


    dump_specset( sout, m_sps_seqs, "sps_seqs", web);
    dump_specset( sout, m_contigsSpectra, "contig_spectra", web);
    //dump_specset( sout, m_starSpectra, "star spectra", web);
  };

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
