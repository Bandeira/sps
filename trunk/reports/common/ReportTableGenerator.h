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

#include "spsFiles.h"

////////////////////////////////////////////////////////////////////////////////
#define OK      1
#define ERROR  -1

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
// Table data
//
struct ReportInternalData {

  int protein;
  int contig;
  int cluster;
  int inputSpectra;

  string proteinText;
  string contigText;
  string clusterText;
  string inputSpectraText;

  int consensus;

  // spectra (consensus or input spectra, depending on the presence of cluster layer) for a contig
  vector<int> spectraOfContig;

  // number of contigs that match a protein
  int matchedContigs;

  // consensus counter (per protein)
  int consensusAcc;

  int matchedContig;
  int allContigsContig;

  int homolog;

  // protein data

  string proteinName;
  string proteinDesc;

  // contig data

   // sequence structure to hold reference and homolog post-processed sequences
  SequenceMapping sequenceMappingReference, sequenceMappingHomolog, sequenceMappingDeNovo;

  // cluster and spectra sequences (final string)
  string clusterSequenceReference, clusterSequenceHomolog, clusterSequenceDeNovo;

  // Reference sequence to be used in tables
  string clusterSequenceReferenceEfective;

  // offset masses
  double offsetHomolog, offsetReference;

  // referentes masses
  vector<float> m_referenceMasses;
  // referentes masses
  vector<float> m_homologMasses;
  // referentes masses
  vector<float> m_deNovoMasses;
  // referentes masses
  vector<float> m_userMasses;

  // true if reference and homolog are equal at cluster level.
  bool auxRefCmp;

  // statitics data
  float cluster_B;
  float cluster_Y;
  float cluster_BYint;

  float spectra_B;
  float spectra_Y;
  float spectra_BYint;


  string tool;

};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class ReportTableGenerator {

  // input spectra files
  vector<string> m_inputSpectraFiles;

  // initial data
  string m_fn_abruijn;
  string m_fn_star;
  string m_fn_seqs;

  //project directory
  string m_projectDir;
  // annotation model directory
  string m_annotationModelDirectory;
  // annotation model filename
  string m_annotationModel;
  // annotation model directory (Prm)
  string m_annotationModelDirectoryPrm;
  // annotation model filename (Prm)
  string m_annotationModelPrm;

  // mass shift value
  float m_massShift;
  // mass shift value (Prm)
  float m_massShiftPrm;
  // peak mass tolerance
  float m_peakMassTol;
  // parent mass tolerance
  float m_parentMassTol;
  // resolution
  float m_resolution;

  // job information
  string m_jobName;
  string m_userName;

  // tool used
  int m_tool;

  // if cluster layer is used
  bool m_noClusters;

  // total number of clusters in clusterData.bin
  int m_clusterCount;

  // total number of spectra mapped in clusterData.bin
  int m_spectraCount;

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
  // Methods to load data files

  // filename subtraction helper method - converts absolute path to relative
  virtual string pathAbsoluteToRelative(const string &projectDir, const string &fileName);


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build and edit the table

  // Add the table headings
  virtual void buildTableHeadings(void);
  //  Build the initial page table
  virtual int buildTableHeader(ReportInternalData &data);
  // Build proteins table
  virtual int buildTableProteins(ReportInternalData &data);
  // Build protein detais (civerage) table
  virtual int buildTableProteinCoverage(ReportInternalData &data);
  // build contigs table - add contigs that map to a protein
  virtual int buildTableContigs(ReportInternalData &data);
  // build contigs table - add contigs that do not map to proteins
  virtual int buildTableContigsOrphan(ReportInternalData &data);
  // add a single contig to the table
  virtual int buildTableContig(ReportInternalData &data);
  // build clustar consensus table
  virtual int buildTableCluster(ReportInternalData &data);
  // build input spectra table
  virtual int buildTableInputSpectra(ReportInternalData &data);
  // build input spectra table when there is no cluster layer
  virtual int buildTableInputSpectra2(ReportInternalData &data);


  //////////////////////////////////////////////////////////////////////////////
  // Methods for data access

  // get input spectra original file and spectra indices from combined spectra index. Used when consensus cluster is not present
  void getinputSpectraFromCombinedSpectra(int &fileIndex, int &spectraIndex, int combinedSpectraIndex);
  // Get the cluster consensus index, given the file index and input spectra index
  int  getConsensusFromInputSpectra(int fileIndex, int inputSpectra);
  // get the file index and input spectra index pair list given the cluster consensus index
  list<pair<unsigned,unsigned> > *getInputSpectraFromConsensus(int consensus);
  // get the contig index given the cluster consensus index
  int  getContigFromConsensus(int consensus);
  // get the total number of clusters
  int getClusterCount(void);
  // get the total number of spectra
  int getSpectraCount(void);
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
  // Get a protein description given it's index
  string getProteinDescription(int protein);
  // Get a protein sequence given it's index
  string getProteinSequence(int protein);
  // replaces | char by space
  void cleanProteinName(string &proteinName);

  // get the model name
  string getModelName(Spectrum *spectrum);


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
  int getSequenceHomolog(ReportInternalData &data);
  // Generate Reference sequence, given the contig index
  int getSequenceReference(ReportInternalData &data);
  // translate a sequence in SequenceMapping structure into a string sequence
  string translateSequence(SequenceMapping &sequenceMapping);
  // translate a sequence form a string to the SequenceMapping structure
  void translateSequenceReverse(SequenceMapping &oSequenceMapping, string &iSequence);
  // Remove the () in sequences like MM(MM)MM, where M is any aa
  bool stringExcessParamsFromSequence(string &iSequence);


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
  int getSequence(int contigIdx, specnets::SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, specnets::Spectrum &contigSpectra, SequenceMapping &sequenceMapping, double &leadMass);
  // First step in common sequence generation
  int generateSequenceStep0( SequenceMapping &proteinData, int proteinIndex);
  // Second step in common sequence generation
  int generateSequenceStep1(specnets::SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, specnets::Spectrum &contigSpectra, ContigMatchData  &contigInfo, int contigIdx);
  // Third step in common sequence generation
  int generateSequenceStep2(int ProteinIndex, SequenceMapping &proteinData, ContigMatchData &source, SequenceMapping &m_sequenceMapping);
  // Third step in protein coverage generation
  int generateSequenceStep2B(int contigIdx, SequenceMapping &proteinData, ContigMatchData &source);

  // Generate statistics (B%, Y% and BY intensity %)
  void generateStatistics(specnets::Spectrum &spectrum, string &peptide, float &m_B, float &m_Y, float &m_BYint);
  //void generateStatistics(specnets::Spectrum &spectrum, string &peptide);


  // get input spectra based on GenoMS options
  list<pair<unsigned,unsigned> > *getInputSpectraFromConsensus(ReportInternalData &data);

  // build mass instervals
  void buildMassIntervals(vector<float> &masses, SequenceMapping &sm);


 public:


  //////////////////////////////////////////////////////////////////////////////
  // Spectrum data needed for several table fields.

  SpsFiles *spsFiles;


  // Constructors and destructor
  ReportTableGenerator();
  // Destructor deletes colTypes vector
  ~ReportTableGenerator();

  // set data files
  virtual void setSpsFiles(SpsFiles *f) {spsFiles = f;};

  // initialize data. Shoud be invoked after setting sps data files
  virtual void init(const ReportGeneratorData &reportGeneratorData);

  // Builds table from input data
  virtual int buildTables(ReportTableHeader *,ReportTableProtein *, ReportTableProteinCoverage *, ReportTableContig *, ReportTableClusterConsensus *, ReportTableInputSpectra *);

  void dump_contigData(ostream &sout, int contig, int star, int homolog, DiffData &diffData, SequenceMapping &iSequenceMapping);

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
