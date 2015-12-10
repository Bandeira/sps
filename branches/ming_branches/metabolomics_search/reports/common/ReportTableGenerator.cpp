////////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <algorithm>

#include "ReportTableGenerator.h"
#include "utils.h"
#include "Tokenizer.h"

#include "SpectrumAnnotStatistics.h"
#include "spectrum_scoring.h"

#include "Logger.h"
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {

using namespace std;
using namespace specnets;

////////////////////////////////////////////////////////////////////////////////
#define TABLE_SEP_L1  '|'
#define TABLE_SEP_L2  '@'
#define TABLE_SEP_L3  '&'
#define TABLE_SEP_L4  '!'

////////////////////////////////////////////////////////////////////////////////
// Report table base methods
////////////////////////////////////////////////////////////////////////////////
ReportTableGenerator::ReportTableGenerator() :
  m_consensusSpecSet(NULL),
  m_clusterData(NULL),
  m_fasta(NULL),
  m_abruijn(NULL),
  m_starSpectra(NULL),
  m_contigIndices(NULL),
  m_contigNames(NULL),
  m_input_index(NULL),
  m_contigsSpectra(NULL),
  m_homglue_ref_midx(NULL),
  m_homglue_ref_mp(NULL),
  m_sps_seqs(NULL),
  m_contigs_midx(NULL),
  m_contigs_mp(NULL),
  m_homglueMatches(NULL),
  m_homglue_matches_midx(NULL),
  m_homglue_matches_mp(NULL)
{
}
////////////////////////////////////////////////////////////////////////////////
ReportTableGenerator::~ReportTableGenerator()
{
  if(m_consensusSpecSet != NULL)      delete m_consensusSpecSet;
  if(m_clusterData != NULL)           delete m_clusterData;
  if(m_fasta != NULL)                 delete m_fasta;
  if(m_abruijn != NULL)               delete m_abruijn;
  if(m_starSpectra != NULL)           delete m_starSpectra;
  if(m_contigIndices != NULL)         delete m_contigIndices;
  if(m_contigNames != NULL)           delete m_contigNames;
  if(m_input_index != NULL)           delete m_input_index;
  if(m_contigsSpectra != NULL)        delete m_contigsSpectra;
  if(m_homglue_ref_midx != NULL)      delete m_homglue_ref_midx;
  if(m_homglue_ref_mp != NULL)        delete m_homglue_ref_mp;
  if(m_sps_seqs != NULL)              delete m_sps_seqs;
  if(m_contigs_midx != NULL)          delete m_contigs_midx;
  if(m_contigs_mp != NULL)            delete m_contigs_mp;
  if(m_homglueMatches != NULL)        delete m_homglueMatches;
  if(m_homglue_matches_midx != NULL)  delete m_homglue_matches_midx;
  if(m_homglue_matches_mp != NULL)    delete m_homglue_matches_mp;
}
////////////////////////////////////////////////////////////////////////////////
// File loading section.
////////////////////////////////////////////////////////////////////////////////
// Main load method. Checks for a table to load. If non existent, calls the loadData method, and then saves the table
int ReportTableGenerator::loadData(const ReportGeneratorData &reportGeneratorData)
{
  if(loadInputSpectraList(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameInputSpectraList);
  if(loadInputSpectraFiles(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameClusterMS);
  if(loadScanFiles(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameScanFiles);
  if(loadConsensusSpectraFile(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameConsensusSpectra);
  if(loadMsClusterData(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading clusterData.bin");

  if(loadProteinsFile(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameProteins);

  if(loadContigIndices(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameContigIndices);
  if(loadHomglueRefMidx(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueRefMidx);
  if(loadHomglueRefMp(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueRefMp);

  if(loadContigsMidxAll(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameContigsMidx);
  if(loadContigsMpAll(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameContigsMp);
  if(loadSpsSeqs(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameSpsSeqs);

  if(loadStarSpectra(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameStarSpectra);
  if(loadAbruijn(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameAbruijn);
  if(loadContigSpectra(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameContigSpectra);

  if(loadCspsMatchesMp(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueMatchesMp);
  if(loadCspsMatchesMidx(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueMatchesMidx);
  if(loadCspsSpectra(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameHomglueMatches);

  if(loadContigNames(reportGeneratorData) == ERROR)
    WARN_MSG("Error loading " << reportGeneratorData.filenameContigNames);

/*
  loadInputSpectraList(reportGeneratorData);
  loadInputSpectraFiles(reportGeneratorData);
  loadScanFiles(reportGeneratorData);
  loadConsensusSpectraFile(reportGeneratorData);
  loadMsClusterData(reportGeneratorData);
  loadProteinsFile(reportGeneratorData);
  loadContigIndices(reportGeneratorData);
  loadHomglueRefMidx(reportGeneratorData);
  loadHomglueRefMp(reportGeneratorData);
  loadContigsMidxAll(reportGeneratorData);
  loadContigsMpAll(reportGeneratorData);
  loadSpsSeqs(reportGeneratorData);
  loadStarSpectra(reportGeneratorData);
  loadAbruijn(reportGeneratorData);
  loadContigSpectra(reportGeneratorData);
  loadCspsMatchesMp(reportGeneratorData);
  loadCspsMatchesMidx(reportGeneratorData);
  loadCspsSpectra(reportGeneratorData);
*/

  // annotation model directory
  m_annotationModelDirectory = reportGeneratorData.annotationModelDir;
  // annotation model filename
  m_annotationModel = reportGeneratorData.annotationModel;

  // mass shift value
  m_massShift = reportGeneratorData.massShift;
  // peak mass tolerance
  m_peakMassTol = reportGeneratorData.peakMassTol;
  // parent mass tolerance
  m_parentMassTol = reportGeneratorData.parentMassTol;

  m_jobName = reportGeneratorData.job;
  m_userName = reportGeneratorData.user;

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// File name composition
string ReportTableGenerator::composeFileName(const string &projectDir, const string &fileName)
{
  // Compose output path
  string aux = projectDir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  // add filename
  aux += fileName;
  // return composed filename
  return aux;
}
////////////////////////////////////////////////////////////////////////////////
// read original input spectra filename list
int ReportTableGenerator::loadInputSpectraList(const ReportGeneratorData &reportGeneratorData)
{
  if(m_input_index) return OK;
  // create object
  m_input_index = new vector<string>();
  // get filename with full path
  string fn = reportGeneratorData.filenameInputSpectraList;
  if(!reportGeneratorData.absoluteFilenameInputSpectraList)
  fn = composeFileName(reportGeneratorData.projectDir, fn);
  // read the index file
  readFilesFromFile(fn, *m_input_index);
  // remove path from filenames
  for(int i = 0 ; i < m_input_index->size() ; i++)
    (*m_input_index)[i] = (*m_input_index)[i].substr((*m_input_index)[i].find_last_of("/\\")+1);
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// read pklbin spectra files
int ReportTableGenerator::loadInputSpectraFiles(const ReportGeneratorData &reportGeneratorData)
{
  // get filename with full path
  string fn = reportGeneratorData.filenameClusterMS;
  if(!reportGeneratorData.absoluteFilenameClusterMS)
  fn = composeFileName(reportGeneratorData.projectDir, fn);
  // read the index file
  readFilesFromFile(fn, m_inputSpectraPklbin);
  // read each specset file
  for(int i = 0 ; i < m_inputSpectraPklbin.size() ; i++) {
    //fn = composeFileName(reportGeneratorData.projectDir, m_inputSpectraPklbin[i].c_str());
    fn = m_inputSpectraPklbin[i];
    SpecSet currentSpecSet;
    if(currentSpecSet.LoadSpecSet_pklbin(fn.c_str()) == -1) {
      cout << "Couldn't open " << fn << endl;
      return ERROR;
    }
    // store read specset
    m_specSet.push_back(currentSpecSet);
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Load scan numbers files
int ReportTableGenerator::loadScanFiles(const ReportGeneratorData &reportGeneratorData)
{
  // get filename with full path
  string fn = reportGeneratorData.filenameScanFiles;
  if(!reportGeneratorData.absoluteFilenameScanFiles)
  fn = composeFileName(reportGeneratorData.projectDir, fn);
  // read the index file
  readFilesFromFile(fn, m_scanNumberFiles);
  // read each specset file
  for(int i = 0 ; i < m_scanNumberFiles.size() ; i++) {
    //fn = composeFileName(reportGeneratorData.projectDir, m_scanNumberFiles[i].c_str());
    fn = m_scanNumberFiles[i];
    // data holder
    vector<vector<int> > data;
    // load the data
    if(Load_binArray(fn.c_str(), data) < 0) {
      cout << "Couldn't open " << fn << endl;
      return ERROR;
    }
    // store read scan data
    m_specScan.push_back(data);
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Load cluster consensus spectra pklbin file
int ReportTableGenerator::loadConsensusSpectraFile(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_consensusSpecSet)
    return OK;
  // if not, create the object
  m_consensusSpecSet = new SpecSet();
  // get filename with full path
  string fn = reportGeneratorData.filenameConsensusSpectra;
  //if(!reportGeneratorData.absoluteFilenameConsensusSpectra)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // from <sps>/spectra/specs_ms.pklbin
  if(m_consensusSpecSet->LoadSpecSet_pklbin(fn.c_str()) == -1) {
    delete m_consensusSpecSet;
    m_consensusSpecSet = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// load MsCluster data.
int ReportTableGenerator::loadMsClusterData(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_clusterData)
    return OK;
  // if not, create the object
  m_clusterData = new ClusterData();
  // first try to load data in binary format
  return m_clusterData->load(reportGeneratorData.projectDir, true, false);
//  if(m_clusterData->loadData(reportGeneratorData.projectDir) == -1) {
//    // if cant find it, try to read the MsCluster output files
//    if(m_clusterData->loadMsClusterData(reportGeneratorData.projectDir) == -1) {
//      // if fail, there is an error. Possibly the wrong directory
//      return ERROR;
//    }
//    // If successful, save in binary format
//    if(m_clusterData->saveData(reportGeneratorData.projectDir) == -1) {
//      // if fail to save, there is an error. Possibly a permissions problem
//      return ERROR;
//    }
//  }
//  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// (4) Load Star Spectra
int ReportTableGenerator::loadStarSpectra(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_starSpectra)
    return OK;
  // if not, create the object
  m_starSpectra = new SpecSet();
  // get filename with full path
  string fn = m_fn_star = reportGeneratorData.filenameStarSpectra;
  //if(!reportGeneratorData.absoluteFilenameStarSpectra)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // load the data
  if(m_starSpectra->LoadSpecSet_pklbin(fn.c_str()) == -1) {
    delete m_starSpectra;
    m_starSpectra = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// (8) Load Abruijn graph
int ReportTableGenerator::loadAbruijn(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_abruijn)
    return OK;
  // if not, create the object
  m_abruijn = new abinfo_t();
  // get filename with full path
  string fn = m_fn_abruijn = reportGeneratorData.filenameAbruijn;
  //if(!reportGeneratorData.absoluteFilenameAbruijn)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // load the data
  if(Load_abinfo(fn.c_str(), *m_abruijn) == 0) {
    delete m_abruijn;
    m_abruijn = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Load proteins file (FASTA format).
int ReportTableGenerator::loadProteinsFile(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_fasta)
    return OK;
  // if not, create the object
  m_fasta = new DB_fasta();
  // get filename with full path
  string fn = reportGeneratorData.filenameProteins;
  //if(!reportGeneratorData.absoluteFilenameProteins)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // load the data
  if(m_fasta->Load(fn.c_str()) == 0) {
    delete m_fasta;
    m_fasta = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// (5) Load Spectra for contigs that matched a protein
int ReportTableGenerator::loadContigSpectra(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_contigsSpectra)
    return OK;
  // if not, create the object
  m_contigsSpectra = new SpecSet();
  // get filename with full path
  string fn = reportGeneratorData.filenameContigSpectra;
  //if(!reportGeneratorData.absoluteFilenameContigSpectra)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // load the data
  int ret = m_contigsSpectra->LoadSpecSet_pklbin(fn.c_str());
  if(ret == 0) {
    delete m_contigsSpectra;
    m_contigsSpectra = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// where to load contig indices to
int ReportTableGenerator::loadContigIndices(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_contigIndices)
    return OK;
  // if not, create the object
  m_contigIndices = new vector<vector<int> >();
  // get filename with full path
  string fn = reportGeneratorData.filenameContigIndices;
  //if(!reportGeneratorData.absoluteFilenameContigIndices)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // load the data
  if(Load_binArray(fn.c_str(), *m_contigIndices) < 0) {
    delete m_contigIndices;
    m_contigIndices = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// where to load contig indices to
// from <sps>/homology/homglue_ref_midx.pklbin
int ReportTableGenerator::loadHomglueRefMidx(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_homglue_ref_midx)
    return OK;
  // if not, create the object
  m_homglue_ref_midx = new SpecSet();
  // get filename with full path
  string fn = reportGeneratorData.filenameHomglueRefMidx;
  //if(!reportGeneratorData.absoluteFilenameHomglueRefMidx)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // from <sps>/spectra/specs_ms.pklbin
  if(m_homglue_ref_midx->LoadSpecSet_pklbin(fn.c_str()) == -1) {
    delete m_homglue_ref_midx;
    m_homglue_ref_midx = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// where to load contig indices to
// from <sps>/homology/homglue_ref_mp.bin
int ReportTableGenerator::loadHomglueRefMp(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_homglue_ref_mp)
    return OK;
  // if not, create the object
  m_homglue_ref_mp = new vector<vector<int> >();
  // get filename with full path
  string fn = reportGeneratorData.filenameHomglueRefMp;
  //if(!reportGeneratorData.absoluteFilenameHomglueRefMp)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // load the data
  if(Load_binArray(fn.c_str(), *m_homglue_ref_mp) < 0) {
    delete m_homglue_ref_mp;
    m_homglue_ref_mp = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
//
int ReportTableGenerator::loadContigsMidxAll(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_contigs_midx)
    return OK;
  // if not, create the object
  m_contigs_midx = new SpecSet();
  // get filename with full path
  string fn = reportGeneratorData.filenameContigsMidx;
  //if(!reportGeneratorData.absoluteFilenameContigsMidx)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // from <sps>/spectra/specs_ms.pklbin
  if(m_contigs_midx->LoadSpecSet_pklbin(fn.c_str()) == -1) {
    delete m_contigs_midx;
    m_contigs_midx = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// where to load contig indices to
// from <sps>/homology/contigs_mp_all.bin
int ReportTableGenerator::loadContigsMpAll(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_contigs_mp)
    return OK;
  // if not, create the object
  m_contigs_mp = new vector<vector<int> >();
  // get filename with full path
  string fn = reportGeneratorData.filenameContigsMp;
  //if(!reportGeneratorData.absoluteFilenameContigsMp)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // load the data
  if(Load_binArray(fn.c_str(), *m_contigs_mp) < 0) {
    delete m_contigs_mp;
    m_contigs_mp = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// where to load contig indices to
// from <sps>/assembly/sps_seqs.pklbin
int ReportTableGenerator::loadSpsSeqs(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_sps_seqs)
    return OK;
  // if not, create the object
  m_sps_seqs = new SpecSet();
  // get filename with full path
  string fn = m_fn_seqs = reportGeneratorData.filenameSpsSeqs;
  //if(!reportGeneratorData.absoluteFilenameSpsSeqs)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // from <sps>/spectra/specs_ms.pklbin
  if(m_sps_seqs->LoadSpecSet_pklbin(fn.c_str()) == -1) {
    delete m_sps_seqs;
    m_sps_seqs = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// (13) load <sps>/homology/homglue_matches_mp.pklbin
int ReportTableGenerator::loadCspsMatchesMp(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_homglue_matches_mp)
    return OK;
  // if not, create the object
  m_homglue_matches_mp = new vector<vector<int> >();
  // get filename with full path
  string fn = reportGeneratorData.filenameHomglueMatchesMp;
  //if(!reportGeneratorData.absoluteFilenameHomglueMatchesMp)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // load the data
  if(Load_binArray(fn.c_str(), *m_homglue_matches_mp) < 0) {
    delete m_homglue_matches_mp;
    m_homglue_matches_mp = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// (14) load <sps>/homology/homglue_matches_midx.bin
int ReportTableGenerator::loadCspsMatchesMidx(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_homglue_matches_midx)
    return OK;
  // if not, create the object
  m_homglue_matches_midx = new SpecSet();
  // get filename with full path
  string fn = reportGeneratorData.filenameHomglueMatchesMidx;
  //if(!reportGeneratorData.absoluteFilenameHomglueMatchesMidx)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // from <sps>/spectra/specs_ms.pklbin
  if(m_homglue_matches_midx->LoadSpecSet_pklbin(fn.c_str()) == -1) {
    delete m_homglue_matches_midx;
    m_homglue_matches_midx = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// (15) load <sps>/homology/homglue_matches.pklbin
int ReportTableGenerator::loadCspsSpectra(const ReportGeneratorData &reportGeneratorData)
{
  // Test if file already loaded
  if(m_homglueMatches)
    return OK;
  // if not, create the object
  m_homglueMatches = new SpecSet();
  // get filename with full path
  string fn = reportGeneratorData.filenameHomglueMatches;
  //if(!reportGeneratorData.absoluteFilenameHomglueMatches)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // from <sps>/spectra/specs_ms.pklbin
  if(m_homglueMatches->LoadSpecSet_pklbin(fn.c_str()) == -1) {
    delete m_homglueMatches;
    m_homglueMatches = NULL;
    cout << "Couldn't open " << fn << endl;
    return ERROR;
  }
  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// (--) load <sps>/homology/ref_sps_names.txt
int ReportTableGenerator::loadContigNames(const ReportGeneratorData &reportGeneratorData)
{
  if(m_contigNames) return OK;
  // create object
  m_contigNames = new vector<string>();
  // get filename with full path
  string fn = reportGeneratorData.filenameContigNames;
  //if(!reportGeneratorData.absoluteFilenameInputSpectraList)
  //fn = composeFileName(reportGeneratorData.projectDir, fn);
  // read the index file
  readFilesFromFile(fn, *m_contigNames);
  // return ok
  return 1;
}
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTables(ReportTableHeader           *tHeader,
                                      ReportTableProtein          *tProteins,
                                      ReportTableProteinCoverage  *tProteinCoverage,
                                      ReportTableContig           *tContig,
                                      ReportTableClusterConsensus *tClusterConsensus,
                                      ReportTableInputSpectra     *tInputSpectra)
{
  // Header table
  tableHeader           = tHeader;
  // Proteins table
  tableProteins         = tProteins;
  // Protein Coverage table
  tableProteinCoverage  = tProteinCoverage;
  //Contigs table
  tableContigs          = tContig;
  // CLuster consensus table
  tableClusterConsensus = tClusterConsensus;
  // Input spectra table
  tableInputSpectra     = tInputSpectra;

  // build tables
  buildTableProteins();
  // add orphan contigs to tables
  buildTableContigsOrphan();
  // buils table header
  buildTableHeader();


  // sort tables
  tableContigs->sortTable(0);
  tableClusterConsensus->sortTable(0);
  tableInputSpectra->sortTable(1);

  // populate IDs in input spectra table. This operation should be done after sort() to allow dictomic search when updating the table user sequence
  tableInputSpectra->populateIDs();
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableHeader(void)
{
  string links, texts;

  // sort elements for duplicate removal
  if(m_inputSpectraFiles.size())
    sort (m_inputSpectraFiles.begin(), m_inputSpectraFiles.end(), stringSortCompare);

  // duplicate find
  vector<string>::iterator it;

  // using default comparison
  it = unique (m_inputSpectraFiles.begin(), m_inputSpectraFiles.end(), stringUniqueCompare);

  // remove extra items
  m_inputSpectraFiles.resize( it - m_inputSpectraFiles.begin() );

  // cicle thru valid elements
  for(int i = 0 ; i < m_inputSpectraFiles.size() ; i++){

    if(i) {
      texts += TABLE_SEP_L1;
      links += TABLE_SEP_L1;
    }
    texts += m_inputSpectraFiles[i];
    links += parseInt(i);
  }


  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  // row data
  vector<string> dataRow;
  // Column 0: Job name
  dataRow.push_back(m_jobName);
  // Column 1: User name
  dataRow.push_back(m_userName);
  // Column 2: Status
  dataRow.push_back("100");
  // Column 3: Elapsed time
  dataRow.push_back("");
  // Column 4: Log
  dataRow.push_back("");
  // Column 5: Link to group by contig
  dataRow.push_back("contigs");
  // Column 6: Link to proteins
  dataRow.push_back("proteins");
  // Column 7: File name of cluster list
  dataRow.push_back("cluster.txt");
  // Column 8: File names
  dataRow.push_back(texts);
  // Column 9: File links
  dataRow.push_back(links);

  // store the row
  if(tableHeader)
    tableHeader->addDataRow(dataRow);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableProteins(void)
{
  if(!m_fasta) return ERROR;

   // cycle thru all proteins
  for(unsigned i = 0 ; i < m_fasta->IDs.size() ; i++ ) {

    //cout << "Building table proteins: " << i << endl;

    ////////////////////////////////////////////////////////////////////////////
    // Data needed to build this table and downstream tables

    // initialize consensus spectra counter for the protein
    m_consensusAcc = 0;
    // get protein name
    m_proteinName = getProteinName(i); // m_fasta->IDs[i];
    // get protein desc
    if(m_fasta)
      m_proteinDesc = m_fasta->getDesc(i);

    string sequence;
    int    proteinSize;
    int    covered;
    double percent;

    if(m_homglue_ref_midx && m_homglue_ref_mp && m_contigsSpectra)
      getProteinCoverage(i, *m_homglue_ref_midx, *m_homglue_ref_mp, *m_contigsSpectra, sequence, proteinSize, covered);

    percent = 100.0 * (double)covered / (double)proteinSize;

    ////////////////////////////////////////////////////////////////////////////
    // downstream table rows

    // build contig entries for this protein
    buildTableContigs(i);

    //cout << "Building table proteins: contigs done" << endl;

    // if the protein has not matched contigs. ignore it
    if(m_contigsMatchProtein.size() < 1)
      continue;

    // protein coverage
    buildTableProteinCoverage(i);


    ////////////////////////////////////////////////////////////////////////////
    // Table row build

    // row data
    vector<string> dataRow;
    // Column 0: index of protein  (1-based)
    dataRow.push_back(parseInt(i + 1));
    // Column 1: protein name
    dataRow.push_back(m_proteinName);
    // Column 2: protein description
    dataRow.push_back(m_proteinDesc);
    // Column 3: number of contigs
    dataRow.push_back(parseInt(m_contigsMatchProtein.size()));
    // Column 4: number of spectra
    dataRow.push_back(parseInt(m_consensusAcc));
    // Column 5: amino acids
    dataRow.push_back(parseInt(covered));
    // Column 6: coverage %
    dataRow.push_back(parseFloat(percent,1));
    // Column 7: protein sequence
    dataRow.push_back(sequence);

    // store the row
    if(tableProteins)
      tableProteins->addDataRow(dataRow);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Protein Coverage
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text          --> Protein ID
// cells[row][1] -> text          --> Protein name
// cells[row][2] -> text          --> Protein length (AAs)
// cells[row][3] -> text list     --> Protein sequence, separated by |
// cells[row][4] -> Contig data   --> CSPS Contigs, separated by |
// cells[row][5] -> Contig data   --> SPS Contigs, separated by |
//      Contig data: : items separated by &
//        0 -> Contig ID
//        1 -> Contig name
//        2 -> Contig start
//        3 -> Contig end
//        4 -> Contig items
//           Contig Item: items separated by @
//              0 -> Beginning
//              1 -> Span
//              0 -> Content
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableProteinCoverage(int protein)
{
  // where to hold protein structure data
  PdProteinInfo proteinData;

  // Process proteins file
  processProteinsFile(protein, proteinData);

  // Structure to hold CSPS contigs info -- intermidiate step
  std::map<int, std::vector<ContigMatchData> >  m_cspsContigInfo;

  // Structure to hold SPS contigs info -- intermidiate step
  std::map<int, std::vector<ContigMatchData> >  m_spsContigInfo;

  // Process csps files
  // CSPS_CONTIG_MATCHES_IDX  /  CSPS_CONTIG_MATCHES  /  CSPS_CONTIG_SPECTRA
  // /homology/homglue_matches_midx.pklbin  /  /homology/homglue_matches_mp.bin  /  /homology/homglue_matches.pklbin
  processContigFiles(*m_homglue_matches_midx, *m_homglue_matches_mp, *m_homglueMatches, m_cspsContigInfo);

  // Process sps files
  // SPS_CONTIG_MATCHES_IDX / SPS_CONTIG_MATCHES / SPS_CONTIG_SPECTRA
  // /homology/homglue_ref_midx.pklbin  /  /homology/homglue_ref_mp.bin  /  /spectra/contigs.pklbin
  processContigFiles(*m_homglue_ref_midx, *m_homglue_ref_mp, *m_contigsSpectra, m_spsContigInfo);


  // second step csps processing
  //DEBUG_MSG("Entering  processContigs(" << i << ", ...)");
  processContigs(protein, proteinData, proteinData.cspsDetails, m_cspsContigInfo);

  // second step sps processing
  //DEBUG_MSG("Entering  processContigs(" << i << ", ...)");
  processContigs(protein, proteinData, proteinData.spsDetails, m_spsContigInfo);

  // Replace sps contig numbers by their names
  //DEBUG_MSG("Entering  populateContigNames(" << i << ", ...)");
	populateContigNames(proteinData.spsDetails);

  // Generate table field output data
  //DEBUG_MSG("Entering  generateOutputFile2(" << i << ", ...)");
  generateCoverageOutput(proteinData);


  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  // row data
  vector<string> dataRow;

  // Column 0: protein index (1-based)
  dataRow.push_back(parseInt(protein + 1));
  // column 1: protein name
  dataRow.push_back(m_proteinName);
  // column 2: Reference sequence
  dataRow.push_back(parseInt(proteinData.proteinLength));
  // column 3: Protein sequence
  dataRow.push_back(proteinData.proteinSequenceEntryData);
  // column 4: csps contig data
  dataRow.push_back(proteinData.cspsEntryData);
  // column 5: sps contig data
  dataRow.push_back(proteinData.spsEntryData);

  // store the row
  if(tableProteinCoverage)
    tableProteinCoverage->addDataRow(dataRow);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Table Contigs
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text   --> contig index (1-based)
// cells[row][1] -> text   --> number of spectra
// cells[row][2] -> text   --> Reference sequence  --> generateSequence(i, return)
// cells[row][3] -> text   --> Homolog sequence   --> generateSequence(i, return)
// cells[row][4] -> text   --> DeNovo sequence --> generateSequence(i, return)
// cells[row][5] -> text   --> User sequence
// cells[row][6] -> text   --> Protein name
// cells[row][7] -> text   --> Protein description
// cells[row][8] -> text   --> spectra (coma separated indexes)
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableContigs(int protein)
{
  ////////////////////////////////////////////////////////////////////////////
  // Data needed to build this table and downstream tables

  // get homolog list for the protein
  m_contigsMatchProtein.clear();
  getHomologFromProtein(protein, m_contigsMatchProtein);

  // get all contigs for this protein
  m_contigs.clear();
  getAllContigFromContig(m_contigsMatchProtein, m_contigs);

  // cycle thru found contigs
  for(int i = 0 ; i < m_contigs.size() ; i++)
    buildTableContig(protein, i, m_contigsMatchProtein[i], m_contigs[i]);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableContigsOrphan(void)
{
  if(!m_sps_seqs)
    return ERROR;

  m_proteinName = "";
  m_proteinDesc = "";


  for(int i = 0 ; i < m_sps_seqs->size() ; i++)
    if((*m_sps_seqs)[i].size() != 0) {
      //cout << "looking for " << i << " ...";
      int homolog = getContigFromAllContig(i);
      if(homolog == -1) {
        //cout << "adding " << i << endl;
        buildTableContig(-1, i, -1, i);
      }// else
        //cout << "found" << endl;
    }
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableContig(int protein, int contig, int matchedContig, int allContigsContig)
{
  //cout << "Building table contigs: " << i << endl;

  ////////////////////////////////////////////////////////////////////////////
  // Sequences

  // get contig that maps to a protein
  // int homolog = getContigFromAllContig(m_contigs[contig]);
  int homolog = getContigFromAllContig(allContigsContig);
  // clear previous consensus cluster list
  m_consensus.clear();
  // get consensus spectra for this contig
  //getConsensusFromContig(m_contigs[contig], m_consensus);
  getConsensusFromContig(allContigsContig, m_consensus);
  // add concensus counter to accumulator
  m_consensusAcc += m_consensus.size();

  // offset masses
  double offsetHomolog, offsetReference;

  // build homolog sequence
  getSequenceHomolog(homolog, offsetHomolog);
  // build reference sequence
  getSequenceReference(homolog, offsetReference);

  // build deNovo sequence
  if(m_contigsSpectra)
    m_sequenceDeNovo = getSequenceDenovo( (*m_sps_seqs)[allContigsContig]);
    //m_sequenceDeNovo = getSequenceDenovo( (*m_sps_seqs)[contig]);
    //m_sequenceDeNovo = getSequenceDenovo( (*m_contigsSpectra)[matchedContig]);
  // build sequence deNovo structure
  translateSequenceReverse(sequenceMappingDeNovo, m_sequenceDeNovo);

  // translate homolog sequence
  m_sequenceHomolog = translateSequence(sequenceMappingHomolog);
  // translate reference sequence
  m_sequenceReference = translateSequence(sequenceMappingReference);

  // auxiliary variable to set reference sequence
  string auxRef;
  if(m_sequenceReference.compare(m_sequenceHomolog) && m_sequenceReference.length())
    auxRef = m_sequenceReference;


  //cout << " aa ; delta ; colspan ; startPosition ; covered" << endl;
  //for(int j = 0 ; j < sequenceMappingHomolog.processedAA.size() ; j++) {
  //  cout << sequenceMappingHomolog.processedAA[j].aa
  //       << " ; " << sequenceMappingHomolog.processedAA[j].delta
  //       << " ; " << sequenceMappingHomolog.processedAA[j].colspan
  //       << " ; " << sequenceMappingHomolog.processedAA[j].startPosition
  //       << " ; " << sequenceMappingHomolog.processedAA[j].covered
  //       << endl;
  //}

  //cout << "mass offset: reference: " << offsetReference << endl;
  //cout << "mass offset: homolog: " << offsetHomolog << endl;


  ////////////////////////////////////////////////////////////////////////////
  // Build mass intervals

  // get masses for the Reference sequences
  //getMassesCummulative((char*)m_sequenceReference.c_str(), m_referenceMasses, 0.0);

  // get masses for the Homolog sequences
  //getMassesCummulative((char*)m_sequenceHomolog.c_str(), m_homologMasses, 0.0);

  // get masses for the novo sequences
  //getMassesCummulative((char*)m_sequenceDeNovo.c_str(), m_deNovoMasses, 0.0);

  m_referenceMasses.push_back(sequenceMappingReference.startMass);
  for(int j = 0 ; j < sequenceMappingReference.processedAA.size() ; j++)
    m_referenceMasses.push_back(sequenceMappingReference.processedAA[j].massPoint);

  m_homologMasses.push_back(sequenceMappingHomolog.startMass);
  for(int j = 0 ; j < sequenceMappingHomolog.processedAA.size() ; j++)
    m_homologMasses.push_back(sequenceMappingHomolog.processedAA[j].massPoint);

  // process mass intervals
  string referenceMasses, homologMasses, deNovoMasses;

  processMassIntervals(m_referenceMasses, referenceMasses);
  processMassIntervals(m_homologMasses, homologMasses);
  processMassIntervals(m_deNovoMasses, deNovoMasses);


  ////////////////////////////////////////////////////////////////////////////
  // Get contig reversed state

  bool contigReversedFlag = getContigState(homolog);
  string contigReversedTableEntry;
  if(contigReversedFlag)
    contigReversedTableEntry = '1';

  ////////////////////////////////////////////////////////////////////////////
  // Data needed to build this table and downstream tables

  // build consensus list
  string consensusList;
  for(int j = 0 ; j < m_consensus.size() ; j++) {
    if(j != 0)
      consensusList += ',';
    consensusList += parseInt(m_consensus[j]);
  }

  string proteinToStore = (protein == -1 ? "-1" : parseInt(protein + 1));


  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  // row data
  vector<string> dataRow;

  // Column 0: contig index (1-based)
  //dataRow.push_back(parseInt(m_contigs[contig] + 1));
  dataRow.push_back(parseInt(allContigsContig + 1));
  // column 1: protein index (1-based)
  dataRow.push_back(proteinToStore);
  // column 2: number of spectra
  dataRow.push_back(parseInt(m_consensus.size()));
  // column 3: Reference sequence
  dataRow.push_back(auxRef);
  // column 4: Homolog sequence
  dataRow.push_back(m_sequenceHomolog);
  // column 5: DeNovo sequence
  dataRow.push_back(m_sequenceDeNovo);
  // column 6: User sequence
  dataRow.push_back("");
  // column 7: Protein name
  dataRow.push_back(m_proteinName);
  // column 8: Protein description
  dataRow.push_back(m_proteinDesc);
  // column 9: spectra (coma separated indexes)
  dataRow.push_back(consensusList);

  // data for contplot

  // column 10: Reference intervals
  dataRow.push_back(referenceMasses);
  // column 11: Homolog intervals
  dataRow.push_back(homologMasses);

  // column 12: Reference offset
  dataRow.push_back(parseFloat(offsetReference,2));
  // column 13: Homolog offset
  dataRow.push_back(parseFloat(offsetHomolog,2));

  // column 14: reverse flag
  dataRow.push_back(contigReversedTableEntry);

  // file names needed: 15 abruijn, 16 stars and 17 sps_seqs
  dataRow.push_back(m_fn_abruijn);
  dataRow.push_back(m_fn_star);
  dataRow.push_back(m_fn_seqs);

  // store the row
  if(tableContigs)
    tableContigs->addDataRow(dataRow);
  // Build the cluster pages
  buildTableCluster(allContigsContig, protein, homolog);
  //buildTableCluster(m_contigs[contig], protein, homolog);

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Consensus cluster table
//////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text   --> index in specset (1-based)
// cells[row][1] -> text   --> cluster index (used when filtered by cluster consensus)
// cells[row][2] -> text   --> contig ID
// cells[row][3] -> text   --> Reference sequence --> generateSequence()
// cells[row][4] -> text   --> Homolog sequence --> generateSequence()
// cells[row][5] -> text   --> DeNovo sequence --> generateSequence()
// cells[row][6] -> text   --> User sequence
// cells[row][7] -> text   --> mass value, from specs[i].parentMass
// cells[row][8] -> text   --> charge value, from specs[i].parentCharge
// cells[row][9] -> text   --> B%
// cells[row][10] -> text  --> Y&
// cells[row][11] -> text  --> BY Int %
// cells[row][12] -> text  --> spectra file name
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableCluster(int contig, int protein, int homolog)
{
  // cycle thru concensus spectra
  for(int i = 0 ; i < m_consensus.size() ; i++) {

    //cout << "Building table clusters: " << i << endl;

    // consensus index
    int consensus = m_consensus[i];

    Spectrum *spectrum = NULL;

    if(m_consensusSpecSet)
      spectrum = &((*m_consensusSpecSet)[consensus]);

    ////////////////////////////////////////////////////////////////////////////
    // Sequences

    // declare local reference holders
    string sequenceReference, sequenceHomolog, sequenceDeNovo;

    // build reference for cluster
    if(homolog != -1)
      propagateSequence(sequenceMappingReferenceCluster, contig, consensus, homolog, sequenceMappingReference);
    else
      sequenceMappingReferenceCluster.clear();

    // build homolog for cluster
    if(homolog != -1)
      propagateSequence(sequenceMappingHomologCluster, contig, consensus, homolog, sequenceMappingHomolog);
    else
      sequenceMappingHomologCluster.clear();

    // build deNovo for cluster
    propagateSequence(sequenceMappingDeNovoCluster, contig, consensus, homolog, sequenceMappingDeNovo);


    // translate reference sequence
    sequenceReference = translateSequence(sequenceMappingReferenceCluster);
    // translate homolog sequence
    sequenceHomolog = translateSequence(sequenceMappingHomologCluster);
    // translate deNovo sequence
    sequenceDeNovo = translateSequence(sequenceMappingDeNovoCluster);


    // auxiliary variable to set reference sequence
    string auxRef;
    bool auxRefCmp = sequenceReference.compare(sequenceHomolog) && sequenceReference.length();
    if(auxRefCmp)
      auxRef = sequenceReference;

    ////////////////////////////////////////////////////////////////////////////
    // Data needed to build this table

    float spectrumParentMass;
    int   spectrumParentCharge;

//    cout << " +++++++++++++++++++++++++++ " << sequenceHomolog << endl;
//    cout << endl << consensus << endl;
    // generat statistics data;
    if(homolog != -1) {
      if(auxRefCmp) {
        //cout << "sequenceReference" <<   sequenceReference << endl;
        if(m_consensusSpecSet)
        // if different, use reference sequence
          generateStatistics((*m_consensusSpecSet)[consensus], sequenceReference);
      } else {
        //cout << "sequenceHomolog" <<   sequenceHomolog << endl;
        // if equal, use homolog
        if(m_consensusSpecSet)
          generateStatistics((*m_consensusSpecSet)[consensus], sequenceHomolog);
      }
    } else
      generateStatistics((*m_consensusSpecSet)[consensus], sequenceDeNovo);


    if(spectrum)
      spectrumParentMass = spectrum->parentMass;

    if(spectrum)
      spectrumParentCharge = spectrum->parentCharge;


    ////////////////////////////////////////////////////////////////////////////
    // Table row build

    // row data
    vector<string> dataRow;

    // Column 0: consensus index (1-based)
    dataRow.push_back(parseInt(consensus + 1));
    // Column 1: contig (1-based)
    dataRow.push_back(parseInt(contig + 1));
    // Column 2: Protein (1-based)
    dataRow.push_back(parseInt(protein + 1));
    // Column 3: sequence - reference
    dataRow.push_back(auxRef);
    // Column 4: sequence - homolog
    dataRow.push_back(sequenceHomolog);
    // Column 5: sequence - de novo
    dataRow.push_back(sequenceDeNovo);
    // Column 6: sequence - user
    dataRow.push_back("");
    // Column 7: parent mass
    dataRow.push_back(parseFloat(spectrumParentMass, 2));
    // Column 8: charge
    dataRow.push_back(parseInt(spectrumParentCharge));
    // Column 9: B (%)
    dataRow.push_back(parseFloat(m_B, 2));
    // Column 10: Y (%)
    dataRow.push_back(parseFloat(m_Y, 2));
    // Column 11: BY Int (%)
    dataRow.push_back(parseFloat(m_BYint, 2));
    // Column 12: protein name
    dataRow.push_back(m_proteinName);

    // generate the input spectra entries for this cluster
    buildTableInputSpectra(contig, consensus, homolog);

    // store the row
    if(tableClusterConsensus)
      tableClusterConsensus->addDataRow(dataRow);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::buildTableInputSpectra(int contig, int consensus, int homolog)
{
  // sequence holders
  SequenceMapping sequenceMappingReferenceSpectra, sequenceMappingHomologSpectra, sequenceMappingDeNovoSpectra;

  ////////////////////////////////////////////////////////////////////////////
  // Table row build

  // get input spectra list
  list<pair<unsigned,unsigned> > *inputSpectra = getInputSpectraFromConsensus(consensus);

  //cout << endl << "Consensus " << consensus << " ---------------------" << endl;

  if(!inputSpectra)
    return ERROR;

  // build iterator
  list<pair<unsigned,unsigned> >::iterator it = inputSpectra->begin();

  // cycle thru input spectra
  for( ; it != inputSpectra->end() ; it++) {

    //cout << " ; " << it->first << " , " << it->second;

    // get spectrum
    Spectrum &spectrum = (m_specSet[it->first])[it->second];

    ////////////////////////////////////////////////////////////////////////////
    // Sequences

    // declare local reference holders
    string sequenceReference, sequenceHomolog, sequenceDeNovo;

    // translate reference sequence
    if(homolog != -1)
      sequenceReference = translateSequence(sequenceMappingReferenceCluster);
    // translate homolog sequence
    if(homolog != -1)
      sequenceHomolog = translateSequence(sequenceMappingHomologCluster);
   // translate deNovo sequence
    sequenceDeNovo = translateSequence(sequenceMappingDeNovoCluster);

    // auxiliary variable to set reference sequence
    string auxRef;
    bool auxRefCmp = sequenceReference.compare(sequenceHomolog) && sequenceReference.length();
    if(auxRefCmp)
      auxRef = sequenceReference;


    ////////////////////////////////////////////////////////////////////////////
    // Data needed to build this table

    // generat statistics data;
    if(homolog != -1) {
      if(auxRefCmp)
        // if different, use reference sequence
        generateStatistics(spectrum, sequenceReference);
      else
        // if equal, use homolog
        generateStatistics(spectrum, sequenceHomolog);
    } else
      generateStatistics(spectrum, sequenceDeNovo);


    ////////////////////////////////////////////////////////////////////////////
    // Get scan number

    // default is spectrum's scan number
    int scan = spectrum.scan;
    // get scan file for this input spectra file
    if(m_scanNumberFiles.size() > it->first) {
      // get the file index
      const vector<vector<int> > & dataFile = m_specScan[it->first];
      // get the scan number
      if(it->second < dataFile.size())
        if(dataFile[it->second].size())
          scan = dataFile[it->second][0];
    }


    ////////////////////////////////////////////////////////////////////////////
    // Table row build

    // row data
    vector<string> dataRow;

    // Column 0: table index (used for updates)
    dataRow.push_back("");
    // Column 1: spectrum index in specset (1-based)
    dataRow.push_back(parseInt(it->second + 1));
    // Column 2: scan number
    dataRow.push_back(parseInt(scan));
    // column 3: cluster index (1-based)
    dataRow.push_back(parseInt(consensus + 1));
    // column 4: Protein name
    dataRow.push_back(m_proteinName);
    // column 5: pklbin file index
    dataRow.push_back(parseInt(it->first));
    // column 6: pklbin file name
    dataRow.push_back(m_inputSpectraPklbin[it->first]);
    // column 7: Reference sequence
    dataRow.push_back(auxRef);
    // column 8: Homolog sequence
    dataRow.push_back(sequenceHomolog);
    // column 9: DeNovo sequence
    dataRow.push_back(sequenceDeNovo);
    // column 10: User sequence
    dataRow.push_back("");
    // column 11: mass value
    dataRow.push_back(parseFloat(spectrum.parentMass, 2));
    // column 12: charge value
    dataRow.push_back(parseInt(spectrum.parentCharge));
    // column 13: B%
    dataRow.push_back(parseFloat(m_B, 2));
    // column 14: Y&
    dataRow.push_back(parseFloat(m_Y, 2));
    // column 15: BY Int %
    dataRow.push_back(parseFloat(m_BYint, 2));
    // column 16: Original filename
    if(m_input_index)
      dataRow.push_back( (*m_input_index)[it->first] );

    // store the row
    if(tableInputSpectra)
      tableInputSpectra->addDataRow(dataRow);

    if(m_input_index)
      m_inputSpectraFiles.push_back( (*m_input_index)[it->first] );
  }
}
////////////////////////////////////////////////////////////////////////////////
// Dump the abruijn graph to the screen
void ReportTableGenerator::dump_abruijn(ostream &sout, abinfo_t *abruijn, bool web)
{
  if(!abruijn) return;

  sout << "----------------- abruijn graph----------------";
  if(web)
    sout << "<br>";
  sout << endl;

  // DUMP abruijn
  abinfo_t::iterator i0 = abruijn->begin();
  for(; i0 != abruijn->end() ; i0++) {
    sout << "Contig " << i0->first << endl;

    // output first vector in pair
    for(int i = 0 ; i < i0->second.first.first.size() ; i++) {
      sout << "     spectrum: " << i0->second.first.first[i] << "  --  flipped: " << i0->second.first.second[i];
      if(web)
        sout << "<br>";
      sout << endl;
    }

    sout << "     ------ ";
    if(web)
      sout << "<br>";
    sout << endl;

    // output second pair -> vector of pairs of vectors
    // output first vector in pair
    std::vector<std::pair< vector<int>, vector<double> > >::iterator i1 = i0->second.second.begin();
    for(; i1 != i0->second.second.end() ; i1++) {
      sout << "     ------ ";
      if(web)
        sout << "<br>";
      sout << endl;

      for(int i = 0 ; i < i1->first.size() ; i++) {
        sout << "     spectrum: " << i1->first[i] << "  --  peak mass: " << i1->second[i];
        if(web)
          sout << "<br>";
        sout << endl;
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
// Dump cluster data
void ReportTableGenerator::dump_clusterData(ostream &sout, ClusterData *cd, char *title, bool web)
{
  if(!cd) return;

  sout << "----------------- cluster data for " << title << " ----------------";
  if(web)
    sout << "<br>";
  sout << endl;

  cd->dump(sout, web);
}
////////////////////////////////////////////////////////////////////////////////
// Dump the bin_array to the screen
void ReportTableGenerator::dump_binArray(ostream &sout, vector<vector<int> > *a, char *title, bool web)
{
  if(!a) return;

  sout << "----------------- bin array graph for " << title << " ----------------";
  if(web)
    sout << "<br>";
  sout << endl;

  for(int i = 0 ; i < a->size() ; i++) {
    for(int j = 0 ; j < (*a)[i].size() ; j++) {
      if(j) sout << ", ";
      sout << (*a)[i][j];
    }
    if(web)
      sout << "<br>";
    sout << endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
// Dump the specset to the screen
void ReportTableGenerator::dump_specset(ostream &sout, specnets::SpecSet *set, char *title, bool web)
{
  if(!set) return;

  sout << "----------------- specset for " << title << " ----------------";
  if(web)
    sout << "<br>";
  sout << endl;

  for(int i = 0 ; i < set->size() ; i++) {
    sout << i << ": size: " << (*set)[i].size();
    for(int j = 0 ; j < (*set)[i].size() ; j++) {
      if(j) sout << ", ";
      sout << (*set)[i][j][0] << " ; " << (*set)[i][j][1];
    }
    if(web)
      sout << "<br>";
    sout << endl;
  }
}
////////////////////////////////////////////////////////////////////////////////
// Data mapping methods. Return
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::processMassIntervals(vector<float> &mass, string &intervals)
{
  float first = -1.0, second;

  for(int j = 0 ; j < mass.size() ; j++) {

    second = mass[j];

    if(first > second)
      break;

    if(j)
      intervals += TABLE_SEP_L1;

    intervals += parseFloat(second, 1);

    first = second;
  }
}
////////////////////////////////////////////////////////////////////////////////
// Data mapping methods. Return
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getConsensusFromInputSpectra(int fileIndex, int inputSpectra)
{
  if(!m_clusterData) return -1;

  return m_clusterData->getConsensusFromInputSpectra(fileIndex, inputSpectra);
}
////////////////////////////////////////////////////////////////////////////////
list<pair<unsigned,unsigned> > *ReportTableGenerator::getInputSpectraFromConsensus(int consensus)
{
  if(!m_clusterData) return NULL;

  return m_clusterData->getInputSpectraFromConsensus(consensus);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getContigFromConsensus(int consensus)
{
  if(!m_abruijn) return -1;
  if(consensus < 0) return -1;

  abinfo_t::iterator i0 = m_abruijn->begin();
  for( ; i0 != m_abruijn->end() ; i0++)
    // output first vector in pair
    for(int i = 0 ; i < i0->second.first.first.size() ; i++)
      if(i0->second.first.first[i] == consensus) {
        return i0->first;
      }

  return -1;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getConsensusFromContig(int contig, vector<int> &ret)
{
  if(!m_abruijn) return;

  if(contig < 0)
    return;

  abinfo_t::iterator i0 = m_abruijn->find(contig);
  if(i0 == m_abruijn->end())
    return;

  for(int i = 0 ; i < i0->second.first.first.size() ; i++)
    ret.push_back(i0->second.first.first[i]);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getContigFromInputSpectra(int fileIndex, int inputSpectra)
{
  // get consensus
  int consensus = getConsensusFromInputSpectra(fileIndex, inputSpectra);
  // get contig
  return getContigFromConsensus(consensus);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getAllContigFromContig(int contig)
{
  if(!m_contigIndices) return contig; // -1;

  if(contig < 0 || contig >= m_contigIndices->size())
    return -1;

  return (*m_contigIndices)[contig][0];
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getAllContigFromContig(vector<int> &contigs, vector<int> &allContigs)
{
  if(!m_contigIndices) {
    allContigs = contigs;
    return;
  }

  for(int i = 0 ; i < contigs.size() ; i++) {
    if((*m_contigIndices)[contigs[i]][0] >= 0)
      allContigs.push_back((*m_contigIndices)[contigs[i]][0]);
    }
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getContigFromAllContig(int contig)
{
  if(!m_contigIndices) return -1;

  for(int i = 0 ; i < m_contigIndices->size() ; i++)
    if((*m_contigIndices)[i][0] == contig)
      return i;
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getProteinFromHomolog(int homolog)
{
  if(!m_contigs_mp) return -1;
  // determine protein
  int protein = -1;
  if((homolog >= 0) && (homolog < m_contigs_mp->size()))
    protein = (*m_contigs_mp)[homolog][0];
  return protein;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getProteinFromReference(int reference)
{
  if(!m_homglue_ref_mp) return -1;
  // determine protein
  int protein = -1;
  if((reference >= 0) && (reference < m_homglue_ref_mp->size()))
    protein = (*m_homglue_ref_mp)[reference][0];
  return protein;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getHomologFromProtein(int protein, vector<int> &ret)
{
  // m_homglue_ref_mp --> m_contigs_mp
  if(!m_homglue_ref_mp) return;
  // exit empty if protein index invalid
  if(protein < 0) return;

  for(int i = 0 ; i < m_homglue_ref_mp->size() ; i++)  // m_contigs_mp
    if((*m_homglue_ref_mp)[i][0] == protein)
      ret.push_back(i);
}
////////////////////////////////////////////////////////////////////////////////
// Sequence generation methods
////////////////////////////////////////////////////////////////////////////////
// Step 0: get protein sequence into a data structure
int ReportTableGenerator::generateSequenceStep0(SequenceMapping &contigInfo, int proteinIndex)
{
  if(!m_fasta) return ERROR;

  string sequence = m_fasta->getSequence(proteinIndex);
  // cycle thru the sequence
  for(int j = 0 ; j < sequence.size() ; j++) {
    aaCell cell;
    // the aa sequence
    cell.aa = sequence[j];
    // just one
    cell.colspan = 0;
    // its position
    cell.startPosition = j;
    contigInfo.processedAA.push_back(cell);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Step 1: load MatchesIndex and contigSpectra into the data structure
int ReportTableGenerator::generateSequenceStep1(SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, SpecSet &contigSpectra, ContigMatchData &contigInfo, int contigIdx)
{
  // check boundaries
  if(contigIdx < 0 || contigIdx >= contigMatches.size())
    return ERROR;

  // get protein index
  m_proteinIdx = contigMatches[contigIdx][0];
  // get the flipped state
  m_flipped = contigMatches[contigIdx][1];


  // test if there is any data to store
  if(contigMatchesIndex[contigIdx].size() <= 0)
    return OK;

  // if this is the first sequence part, skip empty cells and state sequence beginning
  for(int j = 0 ; j < contigMatchesIndex[contigIdx].size() ; j++) {
    PairContigProteinMassIdx pairItem;
    pairItem.contigMassIdx = (int)contigMatchesIndex[contigIdx][j][0];
    pairItem.proteinMassIdx = (int)contigMatchesIndex[contigIdx][j][1];
    contigInfo.pair.push_back(pairItem);
  }

	// 3rd file: Contig spectra
	if(contigIdx < contigSpectra.size())
    if(contigSpectra[contigIdx].size() > 0) {

      // if this is the first sequence part, skip empty cells and state sequence beginning
   		for(int j = 0 ; j < contigSpectra[contigIdx].size() ; j++) {
        float mass = contigSpectra[contigIdx][j][0];
        contigInfo.contigMass.push_back(mass);
      }
    }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Step 2: Calculate alignment and generate sequence
int ReportTableGenerator::generateSequenceStep2(int contigIdx, SequenceMapping &proteinData, ContigMatchData &source, SequenceMapping &m_sequenceMapping)
{
  // get contig span and location
  int size = source.pair.size();
  if(size) {
    m_sequenceMapping.startPosition = source.pair[0].contigMassIdx;
    m_sequenceMapping.endPosition = source.pair[size-1].contigMassIdx;
  }

/*
  // we may need to add the leading mass
  float leadMass = source.contigMass[0];
  if(fabs(leadMass) > m_peakMassTol ) {
    // cell to hold one data item
    aaCell cell;
    // store mass value
    cell.delta = leadMass;
    // add the cell to the list
    m_sequenceMapping.processedAA.push_back(cell);
  } */

  // Loop thru contig size to generate info cells
  for(int j = 0 ; j < size-1 ; j++) {

    // get the corresponding contig indexes
    int contigIndex = source.pair[j].contigMassIdx;
    int contigIndexAfter = source.pair[j+1].contigMassIdx;
    // get the correspondif protein indexes
    int proteinIndex = source.pair[j].proteinMassIdx;
    int proteinIndexAfter = source.pair[j+1].proteinMassIdx;
    // Get the mass at this point
    float massCsps = source.contigMass[contigIndex];
    // and at the point after
    float massAfter = source.contigMass[contigIndexAfter];

    string sequence;

    // get AA sequence for interval
    //if( (contigIndex >= proteinData.startPosition ) && (contigIndexAfter < proteinData.endPosition) )
      for(int k = proteinIndex ; (k < proteinIndexAfter) && (k < proteinData.processedAA.size()) ; k++)
        sequence += proteinData.processedAA[k].aa;

    // get masses for sequence
    vector<float> masses;
    getMasses((char*)sequence.c_str(), masses);
    float mass = 0.0;
    for(int k = 0 ; k < masses.size() ; k++)
      mass += masses[k];

    // massAux used to calculate mass difference
    float massAux = (massAfter - massCsps) - mass;

    // cell to hold one data item
    aaCell cell;

    // save it
    cell.aa = sequence; //result;
    // store mass difference
    cell.delta = massAux;
    // calculate the colspan value for this cell
    cell.colspan = source.pair[j+1].proteinMassIdx - source.pair[j].proteinMassIdx - 1;
    // and the start position
    cell.startPosition = source.pair[j].contigMassIdx;
    // mass point at cell end
    cell.massPoint = massAfter;

    if(!j)
      m_sequenceMapping.startMass = massCsps;
    // add the cell to the list
    m_sequenceMapping.processedAA.push_back(cell);
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
// Step 2: Calculate protein sequence coverage
int ReportTableGenerator::generateSequenceStep2B(int contigIdx, SequenceMapping &proteinData, ContigMatchData &source)
{
  // get contig span and location
  int size = source.pair.size();
  int minimunIndex = source.pair[0].proteinMassIdx;
  int maximunIndex = source.pair[size-1].proteinMassIdx;

  // Loop thru contig size to generate info cells
  for(int j = 0 ; j < size-1 ; j++) {

    // get the corresponding contig indexes
    int contigIndex = source.pair[j].contigMassIdx;
    int contigIndexAfter = source.pair[j+1].contigMassIdx;
    // get the correspondif protein indexes
    int proteinIndex = source.pair[j].proteinMassIdx;
    int proteinIndexAfter = source.pair[j+1].proteinMassIdx;
    // Get the mass at this point
    float massCsps = source.contigMass[contigIndex];
    // and at the point after
    float massAfter = source.contigMass[contigIndexAfter];

    string sequence;

    // mark the protein's aa's covered by the csps contig
    //if( (contigIndex >= proteinData.startPosition ) && (contigIndexAfter < proteinData.endPosition) )
      for(int k = proteinIndex ; (k < proteinIndexAfter) && (k < proteinData.processedAA.size()) ; k++)
        proteinData.processedAA[k].covered = true;
  }

  return OK;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getProteinName(int protein)
{
  string result;

  if(!m_fasta) return result;

  if((protein >= 0) && (protein < m_fasta->IDs.size()))
    result = m_fasta->IDs[protein];
  cleanProteinName(result);
  return result;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getSequenceDenovo(Spectrum &s)
{
  AAJumps masses(1);
  string sequence;
  masses.getPeptideFromSpectrum(s,sequence);

  return sequence;
}
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::translateSequence(SequenceMapping &sequenceMapping)
{
  // return value holder
  string sequence;

  // cycle thru all sequence items
  for(int i = 0 ; i < sequenceMapping.processedAA.size() ; i++) {

    // test delta over peak mass tolerance
    if(fabs(sequenceMapping.processedAA[i].delta) > m_peakMassTol) {
      if(sequenceMapping.processedAA[i].aa.length()) {
        sequence += "(";
        sequence += sequenceMapping.processedAA[i].aa;
        sequence += ",";
        sequence += parseFloat(sequenceMapping.processedAA[i].delta, 2);
        sequence += ")";
      } else {
        sequence += "[";
        sequence += parseFloat(sequenceMapping.processedAA[i].delta, 2);
        sequence += "]";
      }
    } else {
      if(sequenceMapping.processedAA[i].aa.length() > 1)
        sequence += "(";
      sequence += sequenceMapping.processedAA[i].aa;
      if(sequenceMapping.processedAA[i].aa.length() > 1)
        sequence += ")";
    }
  }

  // return the sequence
  return sequence;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::translateSequenceReverse(SequenceMapping &oSequenceMapping, string &iSequence)
{
  string aa("ARNDCEQGHILKMFPSTWYV");
  string aux;

  oSequenceMapping.clear();

  // store the position after each tag
  size_t lastPosition = 0;
  int index = 0;

  while(lastPosition < iSequence.size()) {

    // store the processed sequence
    aaCell cell;
    cell.startPosition = index++;

    // in case of an AA
    if(aa.find_first_of(iSequence[lastPosition]) != string::npos) {

      cell.aa = iSequence[lastPosition];

    // in case of mass
      lastPosition++;

    } else if(iSequence[lastPosition] == '[') {

      // initialize walker
      size_t currentPosition = lastPosition + 1;
      // find end
      while(iSequence[currentPosition] != ']' && currentPosition < iSequence.size())
        currentPosition++;
      // get contents
      aux = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);
      // translate to a float
      cell.delta = getFloat(aux.c_str());
      lastPosition = currentPosition + 1;


    } else if(iSequence[lastPosition] == '(') {

      // initialize walker
      size_t currentPosition = lastPosition + 1;
      // find comma
      while(iSequence[currentPosition] != ',' && iSequence[currentPosition] != ')' && currentPosition < iSequence.size())
        currentPosition++;
      // get sub sequence contents
      aux = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);

      lastPosition = currentPosition + 1;

      if(iSequence[currentPosition] == ',') {
        // find end
        while(iSequence[currentPosition] != ')' && currentPosition < iSequence.size())
          currentPosition++;
        // get contents
        aux = iSequence.substr(lastPosition+1, currentPosition - lastPosition - 1);
        // translate to a float
        cell.delta = getFloat(aux.c_str());

        lastPosition = currentPosition + 1;
      }

    } else { // ERROR
    }

    // add the cell to the list
    oSequenceMapping.processedAA.push_back(cell);
  }

  oSequenceMapping.startPosition = 0;
  oSequenceMapping.endPosition = oSequenceMapping.processedAA.size() - 1;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSequenceReference(int reference, double &leadMass)
{
  // clear sequence mapping between protein and contig
  sequenceMappingReference.clear();
  // clear reference masses
  m_referenceMasses.clear();
  // test for needed files
  if( (!m_homglue_ref_midx) || (!m_homglue_ref_mp) || (!m_contigsSpectra) ) return -1;
  // generate the sequence structure
  return getSequence(reference, *m_homglue_ref_midx, *m_homglue_ref_mp, *m_contigsSpectra, sequenceMappingReference, leadMass);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSequenceHomolog(int homolog, double &leadMass)
{
  // clear sequence mapping between protein and contig
  sequenceMappingHomolog.clear();
  // clear homolog masses
  m_homologMasses.clear();
  // test for needed files
  if( (!m_contigs_midx) || (!m_contigs_mp) || (!m_contigsSpectra) ) return -1;
  // generate the sequence structure
  return getSequence(homolog, *m_contigs_midx, *m_contigs_mp, *m_contigsSpectra, sequenceMappingHomolog, leadMass);
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getProteinCoverage(int protein,
    SpecSet &contigMatchesIndex,
    std::vector<vector<int> > &contigMatches,
    SpecSet &contigSpectra,
    string &sequence, int &proteinSize, int &covered)
{
  // Protein sequence
  SequenceMapping proteinData;

  int ret;
  // Get the protein sequence
  ret = generateSequenceStep0(proteinData, protein);

  // get contigs associated with the protein
  vector<int> contigsIdx;
  for(int i = 0 ; i < contigMatches.size() ; i++)
    if(contigMatches[i][0] == protein)
      contigsIdx.push_back(i);

  // cycle thru all contigs
  for(int i = 0 ; i < contigsIdx.size() ; i++) {
    // Structure to hold contig info -- intermidiate step
    ContigMatchData  contigInfo;
    // Middle step
    ret = generateSequenceStep1(contigMatchesIndex, contigMatches, contigSpectra, contigInfo, contigsIdx[i]);
    // calculate sequence and
    ret = generateSequenceStep2B(contigsIdx[i], proteinData, contigInfo);
  }

  // build sequence for the protein

  // initial state is sequence unmaked
  bool state = false;
  // hold covered AAs count for percentage calculation
  covered = 0;

  // cycle thru all AAs
  for(int i = 0 ; i < proteinData.processedAA.size() ; i++) {

    if(proteinData.processedAA[i].covered)
      covered++;

    // if there is a state change, introduce the mark
    if(proteinData.processedAA[i].covered != state)
      sequence += TABLE_SEP_L1;

    // add the AA
    sequence += proteinData.processedAA[i].aa;

    // set the state as the current state
    state = proteinData.processedAA[i].covered;
  }

  // set the protein size
  proteinSize = proteinData.processedAA.size();

  return ret;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSequence(int contigIdx,
    SpecSet &contigMatchesIndex,
    std::vector<vector<int> > &contigMatches,
    SpecSet &contigSpectra,
    SequenceMapping &sequenceMapping,
    double &leadMass)
{
  // Structure to hold contig info -- intermidiate step
  ContigMatchData  contigInfo;
  // Protein sequence
  SequenceMapping proteinData;

  // determine protein
  int protein = -1;
  if(contigIdx < contigMatches.size())
    protein = contigMatches[contigIdx][0];

  if(protein < 0)
    return OK;

  int ret;
  // Get the protein sequence
  ret = generateSequenceStep0(proteinData, protein);
  // Middle step
  ret = generateSequenceStep1(contigMatchesIndex, contigMatches, contigSpectra, contigInfo, contigIdx);
  // calculate sequence and
  ret = generateSequenceStep2(contigIdx, proteinData, contigInfo, sequenceMapping);

  //cout << "Protein masses: ";
  //for(int i = 0 ; i < contigInfo.pair.size() ; i++)
  //  cout << contigInfo.pair[i].proteinMassIdx << " ; ";
  //cout << endl;

  //cout << "Masses: ";
  //for(int i = 0 ; i < contigInfo.contigMass.size() ; i++)
  //  cout << contigInfo.contigMass[i] << " ; ";
  //cout << endl;

  // get mass offset
  leadMass = contigInfo.contigMass[0];

  return ret;
}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::getContigState(unsigned contig)
{
  // test for needed files
  if(!m_contigs_mp) return false;

  if( contig < m_contigs_mp->size() && contig >= 0)
    return(*m_contigs_mp)[contig][2];
  return false;
}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::getAbruijnStarState(int contig, int star)
{
  // test for needed files
  if(!m_abruijn) return false;

  // contig in abruijn
  abContig_t &ab = (*m_abruijn)[contig];
  // contig data in abruijn
  abContigState_t &cd = ab.first;

  for(int i = 0 ; i < cd.first.size() ; i++) {
    if(cd.first[i] == star)
      return cd.second[i];
  }

  return false;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getReversedAbruijn(int contig, abContigData_t &nodes)
{
  // test for needed files
  if(!m_abruijn) return;

  // contig in abruijn
  abContig_t &ab = (*m_abruijn)[contig];
  // contig data in abruijn
  abContigData_t &cd = ab.second;

  // reverse the spectra - travel thru the contig nodes
  for(int i =  cd.size() - 1 ; i >= 0 ; i--) {
    // star IDs for node
    vector<int> IDs;
    // data for node
    vector<double> data;
    // cycle thru the data in the nodes
    for(int j = 0 ; j < cd[i].first.size() ; j++) {
      int id = cd[i].first[j];
      IDs.push_back(id);
      double parentMass = (*m_starSpectra)[id].parentMass;
      data.push_back(parentMass - AAJumps::massHion - cd[i].second[j]);
    }
    nodes.push_back(make_pair<vector<int>, vector<double> >(IDs, data));
  }

}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::getStarIndex(vector<pair<int, double> > &data, int star, double &value)
{
  for(int i = 0 ; i < data.size() ; i++)
    if(data[i].first == star) {
      value = data[i].second;
      return true;
    }
  return false;
}
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::getSequenceIndex(SequenceMapping &data, int index)
{
  for(int i = 0 ; i < data.processedAA.size() ; i++)
    if(data.processedAA[i].startPosition == index)
      return i;
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getSequenceBetween(aaCell &ret, SequenceMapping &data, int start, int end)
{
  ret.aa = "";
  ret.delta = 0.0;

  for(int i = start ; i < end && i < data.processedAA.size() ; i++) {
    ret.aa += data.processedAA[i].aa;
    ret.delta += data.processedAA[i].delta;
  }
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::propagateSequence(SequenceMapping &oSequenceMapping, int contig, int star, int homolog, SequenceMapping &iSequenceMapping)
{
  // test for needed files
  if( (!m_abruijn) || (!m_starSpectra) ) return;

  // contig in abruijn
  abContig_t &ab = (*m_abruijn)[contig];
  // contig data in abruijn
  abContigData_t *cd, nodes;
  // data holder
  DiffData diffData;

  // get abruijn nodes for contig in the right direction
  // check reversal state
  if( (homolog != -1) && getContigState(homolog)) {
    // of reversed, reverse the contig
    getReversedAbruijn(contig, nodes);
    // and use that one
    cd = &nodes;
  } else {
    // otherwise, use the original one
    cd = &(ab.second);
  }

  // gather start spectra data from abruijn
  for(int i = 0 ; i < cd->size() ; i++)
    for(int j = 0 ; j < (*cd)[i].first.size() ; j++)
      if((*cd)[i].first[j] == star)
        diffData.abruijnData.push_back(make_pair<int,double>(i, (*cd)[i].second[j]));

  // gather contig data from contig spectra
  for(int i = 0 ; i < (*m_sps_seqs)[contig].size() ; i++)
    diffData.contigData.push_back((*m_sps_seqs)[contig][i][0]);


  // clear sequence container
  oSequenceMapping.clear();

  int lastAaItem = iSequenceMapping.processedAA.size() +  iSequenceMapping.startPosition;
  // get first item
  //int start =  diffData.abruijnData[0].first;
  // variables used for data differences
  double s0 = 0.0, c0 = 0.0, c1, s1, sd;
  int proteinItemAtLastContigIndex = -1;
  bool firstItem = true;
  bool mark = false;

  //cout << contig << " ; " << homolog << " --> ";
  //for(int i = 0 ; i < iSequenceMapping.processedAA.size() ; i++) {
  //  cout << '(';
  //  cout << iSequenceMapping.processedAA[i].aa;
  //  cout << '|' << iSequenceMapping.processedAA[i].startPosition << ')';
  //}
  //cout << " ; ";

  // cycle thru all star data (masses)
  for(int i = 0 ; i < diffData.abruijnData.size() ; i++) {
    // get contig node index
    int i0 = diffData.abruijnData[i].first;
    // get star data (mass)
    s1 = diffData.abruijnData[i].second;
    // get contig data (mass)
    c1 = diffData.contigData[i0];
    // calc difference of differences
    sd = (s1-s0) - (c1-c0);
    // declare cell object
    aaCell cell;
    // make sure delta value is 0
    cell.delta = 0.0;
    // in case the current star node is:

    int proteinItemAtContigIndex = getSequenceIndex(iSequenceMapping, i0);

    // if no sequence item was found for this contig item
    if((proteinItemAtContigIndex == -1) && (proteinItemAtLastContigIndex == -1)) {
      //cout << "0(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";

      // ignore contig and star peaks with no AA correspondence
      continue;

    // the first:
    } else if(firstItem) { //if(i0 <= start) {
      //cout << "1(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";

      // just store the mass. This is the mass prefix
      // the aa string goes empty, the mass value is the interval mass
      cell.delta = s1;
      // say prefix item already done
      firstItem = false;

    // in a hole
    } else if(proteinItemAtContigIndex == -1 && proteinItemAtLastContigIndex != -1 && !mark) {
      //cout << "2(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";

      // get the sequence item
      cell = iSequenceMapping.processedAA[proteinItemAtLastContigIndex];
      // apply mass difference
      cell.delta += sd;

      mark = true;

      oSequenceMapping.processedAA.push_back(cell);

      s0 = s1;
      c0 = c1;

      continue;

    } else if(proteinItemAtContigIndex == -1) {
      //cout << "3(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";
      // mark as gapped and move to next
      mark = true;
      continue;

    // at the end of a hole
    } else if(mark) { //proteinItemAtLastContigIndex == -1) {
      //cout << "4(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";

      mark = false;
      // get the sequence item
      cell = iSequenceMapping.processedAA[proteinItemAtLastContigIndex];
      // apply mass difference
      cell.delta += sd;
      // get sequence between the 2 points from consensus
      getSequenceBetween(cell, iSequenceMapping, proteinItemAtLastContigIndex, proteinItemAtContigIndex);

    // inside sequence scope and star scope
    } else {
      //cout << "5(" << i0 << ';' << proteinItemAtLastContigIndex << ";" << proteinItemAtContigIndex << ") ";

      if(proteinItemAtContigIndex - proteinItemAtLastContigIndex > 1)
        getSequenceBetween(cell, iSequenceMapping, proteinItemAtLastContigIndex, proteinItemAtContigIndex);

      else {
        // get the sequence item
        cell = iSequenceMapping.processedAA[proteinItemAtLastContigIndex];
        // apply mass difference
        cell.delta += sd;
      }
    }

    // store the cell
    oSequenceMapping.processedAA.push_back(cell);

    // update lower mass intervalls items for next iteration as this interaction's upper masses
    proteinItemAtLastContigIndex = proteinItemAtContigIndex;
    s0 = s1;
    c0 = c1;
  }


  // suffix mass -- after star end, until contig end

  // get the last item in abruijn graph
  double ddc = s0; //diffData.abruijnData[diffData.abruijnData.size() - 1].second;
  // subtract from parent mass (along with massMH)
  double ddd = (*m_starSpectra)[star].parentMass - ddc - AAJumps::massMH;

  // if greather tham peak mass tolerance
  if(fabs(ddd) > m_parentMassTol) {
    // declare the cell
    aaCell cell2;
    // set the mass value
    cell2.delta = ddd;
    // store the cell
    oSequenceMapping.processedAA.push_back(cell2);
  }

  //cout << " --> ";
  //for(int i = 0 ; i < oSequenceMapping.processedAA.size() ; i++)
  //  cout << oSequenceMapping.processedAA[i].aa << ".";
  //cout << endl;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::dump_contigData(ostream &sout, int contig, int star, int homolog, DiffData &diffData, SequenceMapping &iSequenceMapping)
{
  // get the flipped state
  int flipped = getContigState(homolog);

  // calc diffs (star)
  double s0 = 0.0, s1;
  for(int i = 0 ; i < diffData.abruijnData.size() ; i++) {
    s1 = diffData.abruijnData[i].second;
    double sd = s1 - s0;
    diffData.abDiff.push_back(sd);
    s0 = s1;
  }

  // calc diffs (contig)
  double c0 = 0.0, c1;
  for(int i = 0 ; i < diffData.contigData.size() ; i++) {
    c1 = diffData.contigData[i];
    double cd = c1 - c0;
    diffData.contigDiff.push_back(cd);
    c0 = c1;
  }


  // diference of differences
  for(int i = 0 ; i < diffData.abDiff.size() ; i++) {
    int i1 = diffData.abruijnData[i].first;
    double df = diffData.abDiff[i] - diffData.contigDiff[i1];
    diffData.diffDiff.push_back(df);
  }

  sout << "---- " << contig << " -----" << star << " ------------" << endl;
  sout << "flipped: " << flipped << endl << endl;

  for(int i = 0 ; i < iSequenceMapping.processedAA.size() ; i++) {
    if(iSequenceMapping.processedAA[i].aa.length()) {
       if(fabs(iSequenceMapping.processedAA[i].delta) < m_peakMassTol)
         sout << iSequenceMapping.processedAA[i].aa;
       else
         sout << "(" << iSequenceMapping.processedAA[i].aa << "," << iSequenceMapping.processedAA[i].delta << ")";
    } else
      sout << "[" << iSequenceMapping.processedAA[i].delta << "]";
  }
  sout << endl << endl;

  sout << "abruijn: ";
  for(int i = 0 ; i < diffData.abruijnData.size() ; i++)
    sout << "(" << diffData.abruijnData[i].first << " ; " << diffData.abruijnData[i].second << ")  ";
  sout << endl << endl;

  sout << "contig: ";
  for(int i = 0 ; i < diffData.contigData.size() ; i++)
    sout << diffData.contigData[i] << " ; ";
  sout << endl << endl;

  sout << "abruijn diff: ";
  for(int i = 0 ; i < diffData.abDiff.size() ; i++)
    sout << diffData.abDiff[i] <<  " ; ";
  sout << endl << endl;

  sout << "config diff: ";
  for(int i = 0 ; i < diffData.contigDiff.size() ; i++)
    sout << diffData.contigDiff[i] << " ; ";
  sout << endl << endl;

  sout << "diff of diffs: ";
  for(int i = 0 ; i < diffData.abDiff.size() ; i++)
    sout << diffData.diffDiff[i] << " ; ";
  sout << endl << endl;
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getMassPrefix(int contig, int consensus, double &prefix, double &suffix)
{
  // test for needed files
  if(!m_abruijn) return;

  abinfo_t::iterator i0 = m_abruijn->find(contig);

  vector<double> masses;

  abContigData_t::iterator i1 = i0->second.second.begin();
  for(; i1 != i0->second.second.end() ; i1++)
    for(int i = 0 ; i < i1->first.size() ; i++)
      if(i1->first[i] == consensus)
        masses.push_back(i1->second[i]);

  // make sure its ordered
  sort(masses.begin(), masses.end());

  prefix = suffix = 0.0;
  if(!masses.size())
    return;

  // first is mass prefix
  prefix = masses[0];
  // last is suffix
  suffix = masses[masses.size() - 1];
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::generateStatistics(Spectrum &spectrum, string &peptide)
{
  psmPtr psm(new PeptideSpectrumMatch);
  //Set spectrum
  psm->m_spectrum = &spectrum;
  //Set PSM on spectrum (probably not necessary, but might make things easier
  spectrum.psmList.push_back(psm);

	MS2ScoringModel model;
	// Define annotation file location
	string annotationFile = m_annotationModelDirectory;
	annotationFile += '/';
	annotationFile += m_annotationModel;
	// Set annotatioln file
	model.LoadModel((char*)annotationFile.c_str());

	string allIons("all");

  float aux = m_massShift - AAJumps::massHion;

  // Annotate the spectrum
  psm->annotate(peptide, allIons, model, aux, aux, m_peakMassTol);

  // note that the 'max allowed jumps' parameter has to be set, but we're not using it for anything.
  AAJumps jumps(1);
  // data holders
  vector<float> prmMasses, srmMasses;
  float peptideMass;
  // Compute masses needed
  jumps.getPRMandSRMMasses(peptide, prmMasses, srmMasses, peptideMass);

  int size = spectrum.size();
  float spectrumBottom = spectrum[0][0];
  float spectrumTop = spectrum[size-1][0];

  //cout << "[" << spectrumBottom << " ; " << spectrumTop << "]" << endl;

  // clear statistics variables
  int b0 = 0;
  int y0 = 0;
  int bt = 0;
  int yt = 0;

  psmPtr currAnnotation = spectrum.psmList.front();

  for(int i = 0 ; i < prmMasses.size() ; i++) {
//    cout << i << " ; " << prmMasses[i];
    if((prmMasses[i] >= spectrumBottom - m_peakMassTol) && (prmMasses[i] <= spectrumTop + m_peakMassTol) ) {
//      cout << " ## ";
      bt++;
      for(register int j = 0 ; j <  size; j++)
        if(currAnnotation->m_peakAnnotations[j].first)
          if( (currAnnotation->m_peakAnnotations[j].second == i+1) && (
              (currAnnotation->m_peakAnnotations[j].first->name.compare("b")   == 0) ||
              (currAnnotation->m_peakAnnotations[j].first->name.compare("b++") == 0) )  ) {
            b0++;
//            cout << " ; " << currAnnotation->m_peakAnnotations[j].first->name;
            break;
        }
    }
//    cout << endl;
  }

  m_B = (float)b0 * 100.0 / (float)bt;


//  cout << " --------------------------- " << endl;
  for(int i = 0 ; i < srmMasses.size() ; i++) {
//    cout << i << " ; " << srmMasses[i];
    if((srmMasses[i] >= spectrumBottom - m_peakMassTol) && (srmMasses[i] <= spectrumTop + m_peakMassTol) ) {
//      cout << " ## ";
      yt++;
      for(register int j = 0 ; j <  size; j++)
        if(currAnnotation->m_peakAnnotations[j].first)
          if( (currAnnotation->m_peakAnnotations[j].second == i+1) && (
              (currAnnotation->m_peakAnnotations[j].first->name.compare("y")   == 0) ||
              (currAnnotation->m_peakAnnotations[j].first->name.compare("y++") == 0) )  ) {
            y0++;
//            cout << " ; " << currAnnotation->m_peakAnnotations[j].first->name;
            break;
          }
      }
//    cout << endl;
  }

  m_Y = (float)y0 * 100.0 / (float)yt;


  // get statistics
  SpectrumAnnotStatistics stats;

  m_BYint = stats.percentExplainedIntensity(*psm, allIons);
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::cleanProteinName(string &proteinName)
{
  for(int i = 0 ; i < proteinName.length() ; i++)
    if(proteinName[i] == '|')
      proteinName[i] = ' ';
}
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::processProteinsFile(int protein, PdProteinInfo &proteinData)
{
  // do we have a proteins' file?
  if(!m_fasta) return false;
  // is the protein index consistent with the proteins' file?
  if(m_fasta->size() <= protein) return false;

  // set protein name
  proteinData.proteinName = getProteinName(protein);
  // get the sequence
  string sequence = m_fasta->getSequence(protein);
  // cycle thru the sequence
  for(int j = 0 ; j < sequence.size() ; j++) {
    aaCell cell;
    // the aa sequence
    cell.aa = sequence[j];
    // just one
    cell.colspan = 0;
    // its position
    cell.startPosition = j;
    proteinData.proteinDetail.processedAA.push_back(cell);
  }
  // set protein size
  proteinData.proteinDetail.startPosition = 0;
  proteinData.proteinDetail.endPosition = proteinData.proteinDetail.processedAA.size();
  proteinData.proteinLength = proteinData.proteinDetail.endPosition - proteinData.proteinDetail.startPosition;
  // set internal protein index value
  proteinData.ProteinIndex = protein;

  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ReportTableGenerator::processContigFiles(SpecSet &contigMatchesIndex, std::vector<vector<int> > &contigMatches, SpecSet &contigSpectra,  std::map<int, std::vector<ContigMatchData> >  &contigs)
{
  // First file: homglue matches
	for(int i = 0 ; i < contigMatches.size() ; i++)
    // if equal to -1, then has no match
	  if(contigMatches[i][0] != -1) {
      // create strcuture for current contig
      ContigMatchData contigMatchData;
      // Store contig index
      contigMatchData.contigIndex = i;
	    // find curreint protein
	    std::map<int, std::vector<ContigMatchData> >::iterator it;
	    it = contigs.find(contigMatches[i][0]);
	    // if key exists, insert the new value
	    if(it != contigs.end()) {
	      it->second.push_back(contigMatchData);
	    } else {
	      // if it doesn't exist, create a new entry.
        // create the vector
        vector<ContigMatchData> aux;
        aux.push_back(contigMatchData);
       // Insert CSPS contig (index) info indexed by protein index
        contigs.insert(pair<int, vector<ContigMatchData> > (contigMatches[i][0], aux));
	    }
  	}

  // Second file: CSPS contig matches
	for(int i = 0 ; i < contigMatchesIndex.size() ; i++)
    // test if there is any data to store
	  if(contigMatchesIndex[i].size() > 0) {

	    // get the contig, by index (i)
      ContigMatchData *contig = getContigData(contigs, i);

      // if this is the first sequence part, skip empty cells and state sequence beginning
      if(contig)
    		for(int j = 0 ; j < contigMatchesIndex[i].size() ; j++) {
  		    PairContigProteinMassIdx pairItem;
  		    pairItem.contigMassIdx = (int)contigMatchesIndex[i][j][0];
  		    pairItem.proteinMassIdx = (int)contigMatchesIndex[i][j][1];
  		    contig->pair.push_back(pairItem);
  		  }
	  }

	// 3rd file: Contig spectra
	for(int i = 0 ; i < contigSpectra.size() ; i++)
	  if(contigSpectra[i].size() > 0) {

	    // get the contig, by index (i)
      ContigMatchData *contig = getContigData(contigs, i);

      // if this is the first sequence part, skip empty cells and state sequence beginning
      if(contig)
    		for(int j = 0 ; j < contigSpectra[i].size() ; j++) {
  		    float mass = contigSpectra[i][j][0];
  		    contig->contigMass.push_back(mass);
  		  }
	  }

  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::processContigs(int ProteinIndex, PdProteinInfo &proteinData, PdProteinDetail &target, std::map<int, std::vector<ContigMatchData> > &source)
{
  // Get csps contig information. If not found for this protein, NULL is returned.
  std::map<int, vector<ContigMatchData> >::iterator it2;
  it2 = source.find(proteinData.ProteinIndex);

  // Add CSPS contig information (if exists)
  if(it2 == source.end())
    return;

  // There may be several contigs for this protein. Must cycle thru them
  for(int i = 0 ; i < (*it2).second.size() ; i++) {

    // get csps contig span and location
    int size = (*it2).second[i].pair.size();
    int minimunIndex = (*it2).second[i].pair[0].proteinMassIdx;
    int maximunIndex = (*it2).second[i].pair[size-1].proteinMassIdx;

    // Define the data structure to hold it and set some initial data
    ProteinDetaisContigInfo contig;
    contig.startPosition = minimunIndex;
    contig.endPosition = maximunIndex;


    // Loop thru contig size to generate info cells
    for(int j = 0 ; j < size-1 ; j++) {

      // get the corresponding contig indexes
      int contigIndex = (*it2).second[i].pair[j].contigMassIdx;
      int contigIndexAfter = (*it2).second[i].pair[j+1].contigMassIdx;

      // get the correspondif protein indexes
      int proteinIndex = (*it2).second[i].pair[j].proteinMassIdx;
      int proteinIndexAfter = (*it2).second[i].pair[j+1].proteinMassIdx;
      // Get the mass at this point

      if(contigIndex >= (*it2).second[i].contigMass.size() || contigIndex >= (*it2).second[i].contigMass.size()) {
        DEBUG_MSG("Contig #" << j << ": Invalid mass index: " << contigIndex << " > " << (*it2).second[i].contigMass.size());
        return;
      }

      float massCsps = (*it2).second[i].contigMass[contigIndex];
      // and at the point after
      float massAfter = (*it2).second[i].contigMass[contigIndexAfter];

      string sequence;

      // get AA sequence for interval
      if( (contigIndex >= proteinData.proteinDetail.startPosition ) && (contigIndexAfter < proteinData.proteinDetail.endPosition) )
        for(int k = proteinIndex ; (k < proteinIndexAfter) && (k < proteinData.proteinDetail.processedAA.size()) ; k++)
          sequence += proteinData.proteinDetail.processedAA[k].aa;

      // get masses for sequence
      vector<float> masses;
      getMasses((char*)sequence.c_str(), masses);
      float mass = 0.0;
      for(int k = 0 ; k < masses.size() ; k++)
        mass += masses[k];

      string result = sequence;

      // massAux used to calculate mass difference
      float massAux = (massAfter - massCsps) - mass;

      if(fabs(massAux) > m_peakMassTol) {
        result = "(";
        result += sequence;
        result += ", ";
        result += parseFloat(massAux, 2);
        result += ")";
      }

      // cell to hold one data item
      aaCell cell;

      // save it
      cell.aa = result;
      // calculate the colspan value for this cell
      cell.colspan = (*it2).second[i].pair[j+1].proteinMassIdx - (*it2).second[i].pair[j].proteinMassIdx - 1;
      // and the start position
      cell.startPosition = (*it2).second[i].pair[j].proteinMassIdx;

      // add the cell to the list
      contig.processedAA.push_back(cell);
    }

    // set the contig name
    contig.name = parseInt((*it2).second[i].contigIndex);

    // set the base 1 index
    contig.base1Idx = (*it2).second[i].contigIndex + 1;

    // add the csps contig info to the csps contig list
    target.insert(pair<int, ProteinDetaisContigInfo > ((*it2).second[i].contigIndex, contig));
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::populateContigNames(PdProteinDetail &contig)
{
  PdProteinDetail::iterator it;
  for(it = contig.begin() ; it != contig.end() ; it++) {
    int contigID = it->first;
   	string name = getContigName(contigID);
    it->second.name = (name.length() > 0 ? name : parseInt(contigID + 1) );
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getContigName(int i)
{
  if(m_contigNames)
  	if( (i>= 0) && i < m_contigNames->size() )
	  	return (*m_contigNames)[i];
  return string("--");
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
ContigMatchData *ReportTableGenerator::getContigData( std::map<int, std::vector<ContigMatchData> > &contig, int contigIndex)
{
  std::map<int, vector<ContigMatchData> >::iterator it;
  // search for contig linked to given contig index
  for( it = contig.begin() ; it != contig.end() ; it++)
    // Cycle thru contig vector (for each protein)
    for(int i = 0 ; i < (*it).second.size() ; i++)
      if( (*it).second[i].contigIndex == contigIndex)
        // return pointer to struct
        return &((*it).second[i]);

  return NULL;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::generateCoverageOutput(PdProteinInfo &proteinData)
{
  // test for proteins with no contig mapped to them, and do not generate them
  if((proteinData.cspsDetails.size() == 0) && (proteinData.spsDetails.size() == 0))
    return 1;

  // Get the protein sequence. From this, we know the lenght
  ProteinDetaisContigInfo & proteinSequence = proteinData.proteinDetail;

  // Build a map key index. This is used to maintain the contig index order when outputing them under the protein sequence
  vector<int> spsID, cspsID;
  // get the csps contig indexes
  getOrder(proteinData.cspsDetails, cspsID);
  // get the sps contig indexes
  getOrder(proteinData.spsDetails, spsID);

  // generate protein sequence
  for(int j = 0 ;  j < proteinSequence.processedAA.size() ; j++) {
    if(j)
      proteinData.proteinSequenceEntryData += TABLE_SEP_L1;
    proteinData.proteinSequenceEntryData += proteinSequence.processedAA[j].aa;
  }

  // Add CSPS contig information (if exists)
  generateOutputContig(proteinData, cspsID, proteinData.cspsEntryData, proteinData.cspsDetails, false);

  // Add SPS contig information (if exists)
  generateOutputContig(proteinData, spsID, proteinData.spsEntryData, proteinData.spsDetails, true);

	return 1;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
int ReportTableGenerator::generateOutputContig(PdProteinInfo &proteinData, vector<int> &vectorID, std::string &page, PdProteinDetail &contig, bool link)
{
  // Add CSPS contig information (if exists)
  for(int j = 0 ; j < vectorID.size() ; j++)  {
    // get the contig sequence info
    ProteinDetaisContigInfo & contigSequence = contig[vectorID[j]];

    // get protein size
    int contigSize =  contigSequence.processedAA.size();

    // add contig separator
    if(j) page += TABLE_SEP_L1;

    // Write the contig id and link
    if(link) {
      // add ID
      page += getIntFromSeqName(contigSequence.name);
      page += TABLE_SEP_L2;
      // add name
      page += contigSequence.name;
      page += TABLE_SEP_L2;
    } else {
      // add id
      page += parseInt(contigSequence.base1Idx);
      page += TABLE_SEP_L2;
      // add name
      page += parseInt(contigSequence.base1Idx);
      page += TABLE_SEP_L2;
    }

    // start position
    page += parseInt(contigSequence.startPosition);
    page += TABLE_SEP_L2;
    // end position
    page += parseInt(contigSequence.endPosition);
    page += TABLE_SEP_L2;

    // cycle thru
    for(int k = 0 ; k < contigSize ; k++) {

      // add contig item separator
      if(k) page += TABLE_SEP_L3;

      // contig item start position
      page += parseInt(contigSequence.processedAA[k].startPosition);
      page += TABLE_SEP_L4;
      // contig item colspan
      page += parseInt(contigSequence.processedAA[k].colspan);
      page += TABLE_SEP_L4;
      // contig item content
      page += contigSequence.processedAA[k].aa;
      //page += TABLE_SEP_L4;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
void ReportTableGenerator::getOrder(PdProteinDetail &contig, vector<int> &order)
{
  // contig iterator
  PdProteinDetail::iterator it;
  // cycle thru all contigs
  for( it = contig.begin() ; it != contig.end() ; it++ ) {
    // get the contig sequence info
    ProteinDetaisContigInfo & contigSequence = it->second;
    // Check for empty contigs
    if(!contigSequence.processedAA.size()) continue;
    // add item
    order.push_back(it->first);
  }

  // Sort
  sort(order.begin(), order.end());
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
string ReportTableGenerator::getIntFromSeqName(string seq)
{
  vector<string> aux;
  stringSplit(seq, aux, ":");
  if(aux.size() > 1)
    return aux[1];
  return "";
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
