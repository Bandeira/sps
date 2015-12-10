////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_GENERATOR_H__
#define __REPORT_TABLE_GENERATOR_H__
////////////////////////////////////////////////////////////////////////////////
#include <iterator>
#include <vector>

#include "spectrum.h"
#include "ClusterData.h"
#include "db_fasta.h"
#include "abruijn.h"
#include "aminoacid.h"

#include "ReportTableInputSpectra.h"
#include "ReportTableContig.h"
#include "ReportTableClusterConsensus.h"
#include "ReportTableProteinCoverage.h"
#include "ReportTableProtein.h"
#include "ReportTableHeader.h"

#include "Defines.h"

#include "ReportData.h"

////////////////////////////////////////////////////////////////////////////////
#define OK      1
#define ERROR  -1

#define FILE_SEPARATOR  ";"
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {

using namespace specnets;
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Helper methods for sequence generation
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Holds a contig/protein mass index pair unit
struct PairContigProteinMassIdx {
  int proteinMassIdx;
  int contigMassIdx;
};


// holds data for CSPS contigs
struct ContigMatchData {
  std::vector<float>                    contigMass;
  std::vector<PairContigProteinMassIdx> pair;
  int                                   contigIndex;
};

////////////////////////////////////////////////////////////////////////////////
// Ordering vector element. Used to estabelish contig order when rendering protein coverage report
struct ContigOrdering {
  int contigIndex;
  int startIndex;
  int endIndex;

  // Used by sort method. First order by beggining index, then by ending index.
  bool operator<(const ContigOrdering &o) const
  {return (startIndex == o.startIndex ? endIndex < o.endIndex : startIndex < o.startIndex);};

};
////////////////////////////////////////////////////////////////////////////////
struct aaCell {
  std::string   aa;             // aa sequence for this interval
  double        delta;          // mass difference from the expected
  int           colspan;        // for how many mass intervals it spans
  int           startPosition;  // relative index position
  bool          covered;        // used when calculating protein sequence coverage
  double        massPoint;      // star mass value at this point

  aaCell() : delta(0.0), covered(false), colspan(0) {};
};
////////////////////////////////////////////////////////////////////////////////
// Structure to hold a contig/protein info in protein details context
struct ProteinDetaisContigInfo {
	// Cell start position
  int                       startPosition;
  // Cell end position
  int                       endPosition;
  // Processed AAs for page generation
  std::vector<aaCell>       processedAA;
	// Contig name
  std::string								name;
  //base 1 contig index
  int                       base1Idx;
};


// contig information indexed by contig index
typedef std::map<int, ProteinDetaisContigInfo> PdProteinDetail;



// Holds information for an entire protein and it's contigs
struct PdProteinInfo {

	// protein sequence
	ProteinDetaisContigInfo proteinDetail;

	// Holds information for sps contigs.
	PdProteinDetail spsDetails;

	// Holds information for csps contigs.
	PdProteinDetail cspsDetails;

	// actual protein index
	int ProteinIndex;
	// protein length
	int proteinLength;
  // protein name
  string proteinName;
  // Structure to hold protein table section
  std::string               proteinSequenceEntryData;
  // Structure to hold csps generated table section
  std::string 							cspsEntryData;
  // Structure to hold sps generated table section
  std::string 							spsEntryData;
};
////////////////////////////////////////////////////////////////////////////////
struct SequenceMapping {
	// Cell start position
  int                       startPosition;
  // Cell end position
  int                       endPosition;
  // Processed AAs for page generation
  std::vector<aaCell>       processedAA;
	// Contig name
  std::string								name;
  //base 1 contig index
  int                       base1Idx;
  // sequence 1st mass point
  double                    startMass;

  void clear() {base1Idx=startPosition=endPosition=-1;name.clear();processedAA.clear();};
};
////////////////////////////////////////////////////////////////////////////////
struct DiffData {
  vector<pair<int, double> >  abruijnData;
  vector<double>              contigData;

  vector<double>              abDiff;
  vector<double>              contigDiff;
  vector<double>              diffDiff;
};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class ReportTableGenerator {

  // Structure to hold contig info -- intermidiate step
  //ContigMatchData  m_contigInfo;
  // Sequence mapping between protein and contig
  //SequenceMapping m_sequenceMapping;
  // protein index
  int m_proteinIdx;
  // flipped state
  bool m_flipped;

  // contig IDs that match a protein for a protein
  vector<int> m_contigsMatchProtein;
  // contig IDs for a protein
  vector<int> m_contigs;
  // consensus spectra for a contig
  vector<int> m_consensus;

  // consensus counter
  int m_consensusAcc;

  // sequence structure to hold reference and homolog post-processed sequences
  SequenceMapping sequenceMappingReference, sequenceMappingHomolog, sequenceMappingDeNovo;
  // sequence structure to hold reference and homolog post-processed sequences at cluster level
  SequenceMapping sequenceMappingReferenceCluster, sequenceMappingHomologCluster, sequenceMappingDeNovoCluster;

  // reference sequence
  string m_sequenceReference;
  // homolog sequence
  string m_sequenceHomolog;
  // deNovo sequence
  string m_sequenceDeNovo;

  // referentes masses
  vector<float> m_referenceMasses;
  // referentes masses
  vector<float> m_homologMasses;
  // referentes masses
  vector<float> m_deNovoMasses;
  // referentes masses
  vector<float> m_userMasses;


  // mass prefix and suffix
  double m_prefixMass;
  double m_suffixMass;


  vector<string> m_inputSpectraFiles;


  // initial data
  string m_fn_abruijn;
  string m_fn_star;
  string m_fn_seqs;


  // protein name
  string m_proteinName;
  // protein description
  string m_proteinDesc;

  // annotation model directory
  string m_annotationModelDirectory;
  // annotation model filename
  string m_annotationModel;

  // mass shift value
  float m_massShift;
  // peak mass tolerance
  float m_peakMassTol;
  // parent mass tolerance
  float m_parentMassTol;

  // job information
  string m_jobName;
  string m_userName;


  // statitics data
  float m_B;
  float m_Y;
  float m_BYint;


 protected:


  //////////////////////////////////////////////////////////////////////////////
  // Table pointers.

  // Header table
  ReportTableHeader *tableHeader;
  // Proteins table
  ReportTableProtein *tableProteins;
  // Proteins table
  ReportTableProteinCoverage *tableProteinCoverage;
  //Contigs table
  ReportTableContig *tableContigs;
  // CLuster consensus table
  ReportTableClusterConsensus *tableClusterConsensus;
  // Input spectra table
  ReportTableInputSpectra *tableInputSpectra;


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

  // filename composition helper method
  virtual string composeFileName(const string &projectDir, const string &fileName);

  // load the input spectra files list file
  virtual int loadInputSpectraList(const ReportGeneratorData &reportGeneratorData);
  // (1) Load input spectra
  virtual int loadInputSpectraFiles(const ReportGeneratorData &reportGeneratorData);
  // load scan numbers files
  virtual int loadScanFiles(const ReportGeneratorData &reportGeneratorData);
  // (3) load MsCluster data
  virtual int loadMsClusterData(const ReportGeneratorData &reportGeneratorData);
  // (2) Load cluster consensus spectra file
  virtual int loadConsensusSpectraFile(const ReportGeneratorData &reportGeneratorData);

  // (6) Load contig indices
  virtual int loadContigIndices(const ReportGeneratorData &reportGeneratorData);

  // (4) Load Star Spectra
  virtual int loadStarSpectra(const ReportGeneratorData &reportGeneratorData);
  // (8) Load Abruijn graph
  virtual int loadAbruijn(const ReportGeneratorData &reportGeneratorData);

  // (5) Load Spectra for contigs that matched a protein
  virtual int loadContigSpectra(const ReportGeneratorData &reportGeneratorData);
  // (11) load <sps>/homology/homglue_ref_mp.bin
  virtual int loadHomglueRefMp(const ReportGeneratorData &reportGeneratorData);
  // (12) load <sps>/homology/homglue_ref_midx.pklbin
  virtual int loadHomglueRefMidx(const ReportGeneratorData &reportGeneratorData);

  // (7) load <sps>/assembly/sps_seqs.pklbin
  virtual int loadSpsSeqs(const ReportGeneratorData &reportGeneratorData);
  // (9) load <sps>/homology/contigs_mp_all.bin
  virtual int loadContigsMpAll(const ReportGeneratorData &reportGeneratorData);
  // (10) load <sps>/homology/contigs_midx_all.pklbin
  virtual int loadContigsMidxAll(const ReportGeneratorData &reportGeneratorData);

  // (13) load <sps>/homology/homglue_matches_mp.pklbin
  virtual int loadCspsMatchesMp(const ReportGeneratorData &reportGeneratorData);
  // (14) load <sps>/homology/homglue_matches_midx.bin
  virtual int loadCspsMatchesMidx(const ReportGeneratorData &reportGeneratorData);
  // (15) load <sps>/homology/homglue_matches.pklbin
  virtual int loadCspsSpectra(const ReportGeneratorData &reportGeneratorData);

  // (16) Load proteins file (FASTA format)
  virtual int loadProteinsFile(const ReportGeneratorData &reportGeneratorData);

  // (17) Load contig names
  virtual int loadContigNames(const ReportGeneratorData &reportGeneratorData);



  //////////////////////////////////////////////////////////////////////////////
  // Methods to build and edit the table

  //  Build the initial page table
  virtual int buildTableHeader(void);
  // Build proteins table
  virtual int buildTableProteins(void);
  // Build protein detais (civerage) table
  virtual int buildTableProteinCoverage(int protein);
  // build contigs table - add contigs that map to a protein
  virtual int buildTableContigs(int protein);
  // build contigs table - add contigs that do not map to proteins
  virtual int buildTableContigsOrphan(void);
  // add a single contig to the table
  virtual int buildTableContig(int protein, int contig, int matchedContig, int allContigsContig);
  // build clustar consensus table
  virtual int buildTableCluster(int contig, int protein, int homolog);
  // build input spectra table
  virtual int buildTableInputSpectra(int contig, int consensus, int homolog);


  //////////////////////////////////////////////////////////////////////////////
  // Methods for data access

  // Get the cluster consensus index, given the file index and input spectra index
  int  getConsensusFromInputSpectra(int fileIndex, int inputSpectra);
  // get the file index and input spectra index pair list given the cluster consensus index
  list<pair<unsigned,unsigned> > *getInputSpectraFromConsensus(int consensus);
  // get the contig index given the cluster consensus index
  int  getContigFromConsensus(int consensus);
  // Get the cluster consensus index vector given the contig index
  void getConsensusFromContig(int contig, vector<int>&);
  // Get contig index given file index and input spectra index pair
  int  getContigFromInputSpectra(int fileIndex, int inputSpectra);
  // Get the 'contig that maps to a protein' index given the 'all contigs' index
  int  getContigFromAllContig(int contig);
  // Get the 'all contigs' index list given a 'contig that maps to a protein' index list
  void getAllContigFromContig(vector<int> &contigs, vector<int> &allContigs);
  // Get the 'all contigs' index given the 'contig that maps to a protein' index
  int  getAllContigFromContig(int contig);

  // Get a protein index given a 'contig that maps to a protein' index
  int  getProteinFromHomolog(int homolog);
  // Get a 'contig that maps to a protein' index given a protein index
  void getHomologFromProtein(int protein, vector<int> &ret);
  // get a protein given a reference index
  int getProteinFromReference(int reference);
  // Get a protein name given it's index
  string getProteinName(int protein);
  // replaces | char by space
  void cleanProteinName(string &proteinName);

  //////////////////////////////////////////////////////////////////////////////
  // Protein Coverage related methods

  // put protein data in internal structure for processing
  bool processProteinsFile(int protein, PdProteinInfo &proteinData);
  // put contig data in internal structure for processing
  bool processContigFiles(SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, SpecSet &contigSpectra,  std::map<int, std::vector<ContigMatchData> >  &contigs);
  // process contig data
  void processContigs(int ProteinIndex, PdProteinInfo &proteinData, PdProteinDetail &target, std::map<int, std::vector<ContigMatchData> > &source);
  // replace contig IDs by contig names, from contig names file
  void populateContigNames(PdProteinDetail &contig);
  // get a specific contig's name
  string getContigName(int i);
  // get a specific's contig data
  ContigMatchData *getContigData( std::map<int, std::vector<ContigMatchData> > &contig, int contigIndex);
  // generate table field for all contigs
  int generateCoverageOutput(PdProteinInfo &proteinData);
  // generate table section for a contig
  int generateOutputContig(PdProteinInfo &proteinData, vector<int> &vectorID, std::string &page, PdProteinDetail &contig, bool link);
  // contig sorting
  void getOrder(PdProteinDetail &contig, vector<int> &order);
  // get ID from name
  string getIntFromSeqName(string seq);


  //////////////////////////////////////////////////////////////////////////////
  // Sequence and statistics related methods

  // abruijn contig reverse
  void getReversedAbruijn(int contig, abContigData_t &nodes);
  // abruijn star state (reverse)
  bool getAbruijnStarState(int contig, int star);
  bool getContigState(unsigned contig);
  bool getStarIndex(vector<pair<int, double> > &data, int star, double &value);
  int  getSequenceIndex(SequenceMapping &data, int index);


  // generate deNovo sequence, given a spectrum
  string getSequenceDenovo(specnets::Spectrum &s);
  // Generate homolog sequence, given the contig index
  int getSequenceHomolog(int index, double &leadMass);
  // Generate Reference sequence, given the contig index
  int getSequenceReference(int index, double &leadMass);
  // translate a sequence in SequenceMapping structure into a string sequence
  string translateSequence(SequenceMapping &sequenceMapping);
  // translate a sequence form a string to the SequenceMapping structure
  void translateSequenceReverse(SequenceMapping &oSequenceMapping, string &iSequence);

  // get a sequence for clusters consensus spectrum given a contig sequence
  void propagateSequence(SequenceMapping &oSequenceMapping, int contig, int consensus, int homolog, SequenceMapping &iSequenceMapping);
  // sequence propagation helper method
  void getSequenceBetween(aaCell &ret, SequenceMapping &data, int start, int end);

  // protein coverage calculation
  int getProteinCoverage(int protein, specnets::SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, specnets::SpecSet &contigSpectra, string &sequence, int &proteinSize, int &covered);

  // get mass prefix and suffix for sequences
  void getMassPrefix(int contig, int consensus, double &prefix, double &suffix);

  // process mass intervals as a single string
  void processMassIntervals(vector<float> &mass, string &intervals);


  // Base sequence generation method
  int getSequence(int contigIdx, specnets::SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, specnets::SpecSet &contigSpectra, SequenceMapping &sequenceMapping, double &leadMass);
  // First step in common sequence generation
  int generateSequenceStep0( SequenceMapping &proteinData, int proteinIndex);
  // Second step in common sequence generation
  int generateSequenceStep1(specnets::SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, specnets::SpecSet &contigSpectra, ContigMatchData  &contigInfo, int contigIdx);
  // Third step in common sequence generation
  int generateSequenceStep2(int ProteinIndex, SequenceMapping &proteinData, ContigMatchData &source, SequenceMapping &m_sequenceMapping);
  // Third step in protein coverage generation
  int generateSequenceStep2B(int contigIdx, SequenceMapping &proteinData, ContigMatchData &source);

  // Generate statistics (B%, Y% and BY intensity %)
  void generateStatistics(specnets::Spectrum &spectrum, string &peptide);


 public:


  // Constructors and destructor
  ReportTableGenerator();
  // Destructor deletes colTypes vector
  ~ReportTableGenerator();

  // Load data to build tables. Specific per table; should be overloaded
  virtual int loadData(const ReportGeneratorData &reportGeneratorData);

  // Builds table from input data
  virtual int buildTables(ReportTableHeader *,ReportTableProtein *, ReportTableProteinCoverage *, ReportTableContig *, ReportTableClusterConsensus *, ReportTableInputSpectra *);


  // debug - dump abruijn graph
  void dump_abruijn(ostream &sout, bool web = false) {dump_abruijn(sout, m_abruijn, web);};
  void dump_abruijn(ostream &sout, specnets::abinfo_t *abruijn, bool web);
  void dump_clusterData(ostream &sout, ClusterData *cd, char *title, bool web);
  void dump_specset(ostream &sout, specnets::SpecSet *set, char *title, bool web);
  void dump_binArray(ostream &sout, vector<vector<int> > *a, char *title, bool web);
  void dump_contigData(ostream &sout, int contig, int star, int homolog, DiffData &diffData, SequenceMapping &iSequenceMapping);

  void dump(ostream &sout, bool web = false) {dump_clusterData(sout, m_clusterData, "clusterData.bin", web);
                                              dump_binArray(sout, m_contigs_mp, "contigs_mp", web);
                                              dump_binArray(sout, m_homglue_ref_mp, "homglue_ref_mp", web);
                                              dump_binArray(sout, m_contigIndices, "contigIndices", web);
                                              dump_specset( sout, m_sps_seqs, "sps_seqs", web);
                                              dump_specset( sout, m_contigsSpectra, "contig_spectra", web);
                                              dump_specset( sout, m_contigs_midx, "contigs_midx_all", web);
                                             };

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
