///////////////////////////////////////////////////////////////////////////////
#include <algorithm>

#include "ReportInterface.h"

#include "ParameterList.h"
#include "CommandLineParser.h"

#include "Logger.h"

#include "ReportTableHeader.h"
#include "ReportTableProtein.h"
#include "ReportTableProteinCoverage.h"
#include "ReportTableContig.h"
#include "ReportTableClusterConsensus.h"
#include "ReportTableInputSpectra.h"

#include "ReportHeader.h"
#include "ReportProtein.h"
#include "ReportInputSpectra.h"
#include "ReportContig.h"
#include "ReportCluster.h"

#include "ReportRendererHtml.h"
#include "ReportRendererHtmlDynamic.h"

#include "Timer.h"
#include "Tokenizer.h"

#include "copyright.h"

///////////////////////////////////////////////////////////////////////////////
using namespace specnets;
using namespace std;

namespace spsReports {
///////////////////////////////////////////////////////////////////////////////
ReportInterface::ReportInterface() :
   m_tableNameHeader("tableHeader.txt"),
   m_tableNameProtein("tableProtein.txt"),
   m_tableNameProteinCoverage("tableProteinCoverage.txt"),
   m_tableNameContig("tableContig.txt"),
   m_tableNameCluster("tableCluster.txt"),
   m_tableNameSpectra("tableSpectra.txt"),
   m_verbose(false)
{
}
///////////////////////////////////////////////////////////////////////////////
ReportInterface::~ReportInterface()
{
}
///////////////////////////////////////////////////////////////////////////////
string ReportInterface::composeFileName(const string &projectDir, const string &fileName)
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
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::buildDirectoryPath(const string &dir)
{
  if (! dir.empty())
    for (string::size_type i = 0, j = min(dir.find('/', i), dir.size()); i < dir.size(); i = j + 1, j = min(dir.find('/', i), dir.size()))
      if (! dir.substr(0, j).empty())
        if (mkdir(dir.substr(0, j).c_str(), ACCESSPERMS) == -1)
        //if (mkdir(dir.substr(0, j).c_str()) == -1)
          if (errno != EEXIST) {
            cerr << "error: cannot create directory '" << dir.substr(0, j) << "': " << strerror(errno) << endl;
            return ERROR;
          }
   return OK;
}
///////////////////////////////////////////////////////////////////////////////
void ReportInterface::dump_abruijn(ReportGeneratorData &data)
{
  ReportTableGenerator rtg;
  // load project data
  rtg.loadData(data);
  // generate the tables
  rtg.dump_abruijn(cout);
}
///////////////////////////////////////////////////////////////////////////////
void ReportInterface::dump_binArray(ReportGeneratorData &data)
{
  ReportTableGenerator rtg;
  // load project data
  rtg.loadData(data);
  // generate the tables
  rtg.dump(cout);
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::buildTables(ReportGeneratorData &data)
{
  Timer_c timer;

  // verbose output
  verboseOutput(cout, "Generating tables....", false);

  // build directory path for tables
  buildDirectoryPath(data.tablesDir);

  // set the table filenames
  ReportTableHeader           tHeader(data.tablesDir,           m_tableNameHeader);
  ReportTableProtein          tProtein(data.tablesDir,          m_tableNameProtein);
  ReportTableProteinCoverage  tProteinCoverage(data.tablesDir,  m_tableNameProteinCoverage);
  ReportTableContig           tContig(data.tablesDir,           m_tableNameContig);
  ReportTableClusterConsensus tCluster(data.tablesDir,          m_tableNameCluster);
  ReportTableInputSpectra     tInputSpectra(data.tablesDir,     m_tableNameSpectra);

  // table generator object
  ReportTableGenerator rtg;
  // load project data
  rtg.loadData(data);
  // generate the tables
  rtg.buildTables(&tHeader, &tProtein, &tProteinCoverage, &tContig, &tCluster, &tInputSpectra);

  // save the tables
  tHeader.saveTable();
  tProtein.saveTable();
  tProteinCoverage.saveTable();
  tContig.saveTable();
  tCluster.saveTable();
  tInputSpectra.saveTable();

  // verbose output
  verboseOutput(cout, timer.stop());

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::paginate(ReportGeneratorData &data, ReportRendererBase &rr, ReportBase &rep, string fnamePrefix, string pageName)
{
  Timer_c timer;

  // get the total number of elements, plus the vector of IDs
  vector<string> IDs, IDsEnd, fNames;
  // navigation bar name and suffix
  string fnameSuffix = ".html";
  int nPages = rep.getPageSplitStats(DEFAULT_ROWS_PER_TABLE, IDs, IDsEnd, fNames, fnamePrefix, fnameSuffix);
  // Declare filename holder and navigation bar holder variables
  string navBar;
  // build the navigation bar
  rr.buildNavigationBar(navBar, IDs, IDsEnd, fNames);
  // set the navigation bar
  rr.setNavigationBar(navBar);
  // cycle thru all pages
  for(int j = 0 ; j < nPages ; j++) {
    // verbose output
    stringstream aux; aux << "Generating " << pageName << " page: ";
    verboseOutput(cout, aux.str().c_str(), fNames[j].c_str(), "...", false);
    // add specified target directory
    string fn = composeFileName(data.outDir, fNames[j]);
    // open file to write to
    ofstream of(fn.c_str(), ios::out);
    // build the report page
    rr.renderReport(&rep, of, j * DEFAULT_ROWS_PER_TABLE, DEFAULT_ROWS_PER_TABLE);
    // close output file
    of.close();
    // verbose output (time)
    verboseOutput(cout, timer.restart());
  }
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::generateReportHtmlDynamic(ReportGeneratorData &data)
{
  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  //data.projectDir     = m_projectDir;

  vector<string>      aux;

  // static report renderer for contigs page
  ReportRendererHtmlDynamic rrHtmlDyn;
  // Specify project directory (needed for input datat)
  rrHtmlDyn.setProjectDir(data.projectDir);
  // Specify server cgi-bin location
  rrHtmlDyn.setServerLocation(data.server);
  // Specify cells per line
  rrHtmlDyn.setCellsPerLine(data.cellsPerLine);
  // specify fasta filename
  rrHtmlDyn.setFastaFilename(data.filenameProteins);

  // Specify project directory (needed to get data to generate reports dynamically)
  if(data.targetProjectDir.length())
    rrHtmlDyn.setProjectDirRel(data.targetProjectDir);
  else
    rrHtmlDyn.setProjectDirRel(data.projectDir);



  //////////////////////////////////////////////////////////////////////////////
  // header report generation

  // define report cluster object, and load tables
  ReportHeader rh(data.tablesDir, m_tableNameHeader);
  // load tables
  rh.load();
  // define file name for HTML report page
  string fn = "index.html";
  // add specified target directory
  fn = composeFileName(data.outDir, fn);
  // open file to write to
  ofstream ofh(fn.c_str(), ios::out);
  // build the report page
  rrHtmlDyn.renderReport(&rh, ofh);
  // close output file
  ofh.close();


  //////////////////////////////////////////////////////////////////////////////
  // End

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::generateReportHtml(ReportGeneratorData &data)
{
  // performance timers
  Timer_c timer, timer2;


  //////////////////////////////////////////////////////////////////////////////
  // data needed to generate tables

  //data.projectDir     = m_projectDir;

  vector<string>      aux;

  // static report renderer for contigs page
  ReportRendererHtml rrHtml;
  // Specify executables directory (needed for specplot module)
  rrHtml.setExeDir(data.exeDir);
  // Specify project directory (needed for input data and report output)
  rrHtml.setProjectDir(data.projectDir);

  // specify fasta filename
  rrHtml.setDisplayLevel(data.displayLevel);

  //////////////////////////////////////////////////////////////////////////////
  // proteins report generation

  // Proteins main page

  // verbose output
  verboseOutput(cout, "---- Generating protien pages ----");
  verboseOutput(cout, "Generating protiens page...", false);
  // define report protein list object, and load tables
  ReportProtein pl(data.tablesDir, m_tableNameProtein);
  // load table
  pl.load();
  // static report renderer for contigs page, and add specified target directory
  string fn = composeFileName(data.outDir, "proteins.html");
  // open file to write to
  ofstream of2(fn.c_str(), ios::out);
  // build the report page
  rrHtml.renderReport(&pl, of2);
  // close output file
  of2.close();

  // verbose output (time)
  verboseOutput(cout, timer.restart());


  // Child (specific) protein report pages

  // define specific protein report object
  ReportProtein p2(data.tablesDir, m_tableNameProtein, m_tableNameContig);
  // load table
  p2.load();
  // get the proteins column from the proteins table
  aux = p2.getTableColumn(0, TABLE_PROTEIN_FILTER_COL_PROTEIN);
  // cycle tru all proteins
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    p2.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "protein.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, rrHtml, p2, fnamePrefix, "protein");
  }


  //////////////////////////////////////////////////////////////////////////////
  // contigs report generation

  // Contigs main page

  // verbose output
  verboseOutput(cout, "---- Generating contig pages ----");
  // define report object for contig list
  ReportContig cl(data.tablesDir, m_tableNameContig);
  // load table
  cl.load();
  // generate pages with pagination
  paginate(data, rrHtml, cl, "contigs.", "contigs");


  // Child (specific) contig report pages

  rrHtml.clearNavigationBar();
  // define specific contig report object
  ReportContig rc(data.tablesDir, m_tableNameContig, m_tableNameCluster);
  // load table
  rc.load();
  // get the contig id column from the contigs table
  aux = rc.getTableColumn(0, TABLE_CONTIG_FILTER_COL_CONTIG);
  // cycle tru all contigs
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    rc.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "contig.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, rrHtml, rc, fnamePrefix, "contig");
  }


  //////////////////////////////////////////////////////////////////////////////
  // cluster report generation

  // Consensus cluster clusters page

  verboseOutput(cout, "---- Generating cluster pages ----");

  rrHtml.clearNavigationBar();
  // define report cluster object, and load tables
  ReportCluster clust(data.tablesDir, m_tableNameCluster, m_tableNameSpectra);
  // load tables
  clust.load();
  // get the cluster id column from the cluster table
  aux = clust.getTableColumn(0, TABLE_CLUSTER_FILTER_COL_CLUSTER);
  // cycle tru all clusters
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    clust.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "cluster.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, rrHtml, clust, fnamePrefix, "cluster");
  }


  //////////////////////////////////////////////////////////////////////////////
  // spectra report generation

  // verbose output
  verboseOutput(cout, "---- Generating input spectra pages ----");
  rrHtml.clearNavigationBar();
  // define report cluster object, and load tables
  ReportInputSpectra is(data.tablesDir, m_tableNameSpectra);
  // load tables
  is.load();
  // get the file index column from the input spectra table
  aux = is.getTableColumn(0, TABLE_SPECTRA_FILTER_COL_FILE);
  // sort elements for duplicate removal
  sort (aux.begin(), aux.end(), stringSortCompare);
  // duplicate find
  vector<string>::iterator it;
  // using default comparison
  it = unique (aux.begin(), aux.end(), stringUniqueCompare);
  // remove extra items
  aux.resize( it - aux.begin() );
  // cycle tru all input files
  for(int i = 0 ; i < aux.size() ; i++) {
    // set the filter for a specific contig in cluster table
    is.applyFilter(aux[i]);
    // set filename
    string fnamePrefix = "spectra.";  fnamePrefix +=  aux[i]; fnamePrefix += ".";
    // generate pages with pagination
    paginate(data, rrHtml, is, fnamePrefix, "spectra");
  }


  //////////////////////////////////////////////////////////////////////////////

  verboseOutput(cout, "---- Done! ----");
  // verbose output (time)
  verboseOutput(cout, timer2.stop());


  //////////////////////////////////////////////////////////////////////////////
  // header report generation

  rrHtml.clearNavigationBar();
  // define report cluster object, and load tables
  ReportHeader rh(data.tablesDir, m_tableNameHeader);
  // load tables
  rh.load();
  // define file name for HTML report page
  fn = "index.html";
  // add specified target directory
  fn = composeFileName(data.outDir, fn);
  // open file to write to
  ofstream ofh(fn.c_str(), ios::out);
  // build the report page
  rrHtml.renderReport(&rh, ofh);
  // close output file
  ofh.close();


  //////////////////////////////////////////////////////////////////////////////
  // End

  return OK;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::parseOptions(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------
  // Initialize directories used
  //--------------------------------------------------------------------------------------------
  // Get the execultable directory.
  string str = argv[0];
  size_t found;
  found = str.find_last_of("/\\");
  string exeDir = ".";
  if(found != string::npos)
    exeDir = str.substr(0, found);

  //--------------------------------------------------------------------------------------------
  // Parse the command line parameters
  //--------------------------------------------------------------------------------------------
  vector<CommandLineParser::Option> listOptions;

  //  c++ instruction                            cmd line option    parameter name    takes value?

  listOptions.push_back(CommandLineParser::Option("-help",            "help",             false));
  listOptions.push_back(CommandLineParser::Option("-version",         "VERSION",          false));
  listOptions.push_back(CommandLineParser::Option("-debug",           "debug",            false));
  listOptions.push_back(CommandLineParser::Option("-verbose",         "verbose",          false));

  // Data tables names and project location
  listOptions.push_back(CommandLineParser::Option("-exe-dir",         "EXE_DIR",          true));
  listOptions.push_back(CommandLineParser::Option("-project-dir",     "PROJECT_DIR",      true));
  listOptions.push_back(CommandLineParser::Option("-tables-dir",      "TABLES_DIR",       true));
  listOptions.push_back(CommandLineParser::Option("-table-header",    "TABLE_HEADER",     true));
  listOptions.push_back(CommandLineParser::Option("-table-protein",   "TABLE_PROTEIN",    true));
  listOptions.push_back(CommandLineParser::Option("-table-contig",    "TABLE_CONTIG",     true));
  listOptions.push_back(CommandLineParser::Option("-table-cluster",   "TABLE_CLUSTER",    true));
  listOptions.push_back(CommandLineParser::Option("-table-spectra",   "TABLE_SPECTRA",    true));

  listOptions.push_back(CommandLineParser::Option("-project-dir-server", "PROJECT_DIR_SERVER", true));


  // input data files

  listOptions.push_back(CommandLineParser::Option("-abruijn",           "FILE_COMP",        true));
  listOptions.push_back(CommandLineParser::Option("-sps-seqs",          "FILE_SEQS",        true));
  listOptions.push_back(CommandLineParser::Option("-contigs-mp",        "FILE_MP2",         true));
  listOptions.push_back(CommandLineParser::Option("-contigs-midx",      "FILE_MIDX2",       true));
  listOptions.push_back(CommandLineParser::Option("-homglue-ref-mp",    "FILE_REFMP",       true));
  listOptions.push_back(CommandLineParser::Option("-homglue-ref-midx",  "FILE_REFMIDX",     true));
  listOptions.push_back(CommandLineParser::Option("-stars",             "FILE_STARS",       true));
  listOptions.push_back(CommandLineParser::Option("-consensus-spectra", "FILE_MS",          true));
  listOptions.push_back(CommandLineParser::Option("-input-spectra-list","FILE_INDEX",       true));
  listOptions.push_back(CommandLineParser::Option("-cluster-ms",        "FILE_CLUSTERMS",   true));
  listOptions.push_back(CommandLineParser::Option("-cluster-scan",      "FILE_CLUSTERSCAN", true));
  listOptions.push_back(CommandLineParser::Option("-proteins",          "FILE_FASTA",       true));


  // variables
  listOptions.push_back(CommandLineParser::Option("-title",             "TITLE",            true));
  listOptions.push_back(CommandLineParser::Option("-shift-value",       "SHIFT_VALUE",      true));

  listOptions.push_back(CommandLineParser::Option("-job",               "JOB",              true));
  listOptions.push_back(CommandLineParser::Option("-user",              "USER",             true));
  listOptions.push_back(CommandLineParser::Option("-server",            "SERVER",           true));
  listOptions.push_back(CommandLineParser::Option("-cells-per-line",    "CELLS_PER_LINE",   true));

  listOptions.push_back(CommandLineParser::Option("-aa",                "AMINO_ACID_LIST",  true));
  listOptions.push_back(CommandLineParser::Option("-annotation-model",  "ANNOTATION_MODEL", true));


  listOptions.push_back(CommandLineParser::Option("-peakmasstol",       "TOLERANCE_PEAK",   true));
  listOptions.push_back(CommandLineParser::Option("-parentmasstol",     "TOLERANCE_PM",     true));

  listOptions.push_back(CommandLineParser::Option("-outdir",            "REPORT_DIR",       true));
  listOptions.push_back(CommandLineParser::Option("-prefix",            "PREFIX",           true));
  listOptions.push_back(CommandLineParser::Option("-outfile",           "OUTFILE",          true));
  listOptions.push_back(CommandLineParser::Option("-format",            "FORMAT",           true));
  listOptions.push_back(CommandLineParser::Option("-target",            "TARGET",           true));

  listOptions.push_back(CommandLineParser::Option("-no-msms-images",    "NO_MSMS_IMAGES",   true));


  // parameter file
  listOptions.push_back(CommandLineParser::Option("-p",                 "PARAMETER_FILE",   true));

  // Build tables command
  listOptions.push_back(CommandLineParser::Option("-build-tables",      "buildTables",      false));

  // Build report command
  listOptions.push_back(CommandLineParser::Option("-build-html",        "buildhtml",        false));

  // Build dynamic report command
  listOptions.push_back(CommandLineParser::Option("-build-html-dynamic","REPORT_DYNAMIC",   false));


  ////////////////////////////////////////////////////////////////////////////////
  // Execute the command line parser
  CommandLineParser clp(argc, argv, 0, listOptions);

  ////////////////////////////////////////////////////////////////////////////////
  // Checking for errors
  string parser_error;
  if (!clp.validate(parser_error)) {
    ERROR_MSG(parser_error);
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // The parameters' values
  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);


  // add EXE_DIR parameter
  commandLineParams.addIfDoesntExist("EXE_DIR", exeDir);


  ////////////////////////////////////////////////////////////////////////////////
  // Process the commands
  return processOptions(commandLineParams);

}
//////////////////////////////////////////////////////////////////////////////////
// Process options
//////////////////////////////////////////////////////////////////////////////////
int ReportInterface::processOptions(ParameterList &ip)
{
  ParameterList commandLineParams;

  ////////////////////////////////////////////////////////////////////////////////
  // Parameter file
  ////////////////////////////////////////////////////////////////////////////////

  if (ip.exists("PARAMETER_FILE")) {

    string parameterFilename = ip.getValue("PARAMETER_FILE");

    commandLineParams.readFromFile(parameterFilename);
    // Combine the command line parameters to the file ones
    //   Command line parameters take precedence (hence the overwrite flag set)
  }

  commandLineParams.addList(ip, true);

  ////////////////////////////////////////////////////////////////////////////////
  // help message control
  bool showHelp = true;

  // Data holder
  ReportGeneratorData data;

  ////////////////////////////////////////////////////////////////////////////////
  // "help" prints help and exits
  if (commandLineParams.exists("VERSION"))
    return version(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // the same for "version"
  if (commandLineParams.exists("help"))
    return help(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // sets the "verbose" flag
  if (commandLineParams.exists("verbose"))
    m_verbose = true;

  ////////////////////////////////////////////////////////////////////////////////
  // Boolean parameters

  if (commandLineParams.exists("SHIFT_VALUE"))
    ;


  ////////////////////////////////////////////////////////////////////////////////
  // File locations
  ////////////////////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////////////////////
  // File load section

  // Executables directory
  if (commandLineParams.exists("EXE_DIR"))
    data.exeDir = commandLineParams.getValue("EXE_DIR");

  // Project directory
  if (commandLineParams.exists("PROJECT_DIR")) {
    data.projectDir  = commandLineParams.getValue("PROJECT_DIR");
    data.outDir      = composeFileName(data.projectDir, data.outDir);
    data.tablesDir   = data.projectDir;
  }

  // Server cgi-bin location
  if (commandLineParams.exists("CELLS_PER_LINE"))
    data.cellsPerLine  = getInt(commandLineParams.getValue("CELLS_PER_LINE").c_str());

  // Server cgi-bin location
  if (commandLineParams.exists("SERVER"))
    data.server  = commandLineParams.getValue("SERVER");


  // Tables directory
  if (commandLineParams.exists("TABLES_DIR"))
     data.tablesDir = commandLineParams.getValue("TABLES_DIR");


  // table file names
  if (commandLineParams.exists("TABLE_HEADER"))
    m_tableNameHeader = commandLineParams.getValue("TABLE_HEADER");

  if (commandLineParams.exists("TABLE_PROTEIN"))
    m_tableNameProtein = commandLineParams.getValue("TABLE_PROTEIN");

  if (commandLineParams.exists("TABLE_CONTIG"))
    m_tableNameContig = commandLineParams.getValue("TABLE_CONTIG");

  if (commandLineParams.exists("TABLE_CLUSTER"))
    m_tableNameCluster = commandLineParams.getValue("TABLE_CLUSTER");

  if (commandLineParams.exists("TABLE_SPECTRA"))
    m_tableNameSpectra = commandLineParams.getValue("TABLE_SPECTRA");


  // Font file location
  if (commandLineParams.exists("FONT_LOCATION"))
    ;

  // Amino-acid masses file
  if (commandLineParams.exists("AMINO_ACID_LIST")) {
    string aux = commandLineParams.getValue("AMINO_ACID_LIST");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      // ReportInterfacetrum.setAminoAcidsFile(aux.substr(found + 1));
      // ReportInterfacetrum.setAminoAcidsFileDirectory(aux.substr(0, found));
    } else {
      // ReportInterfacetrum.setAminoAcidsFile(aux);
      // ReportInterfacetrum.setAminoAcidsFileDirectory(".");
    }
  }


  data.annotationModelDir = data.exeDir;


  // Annotation model
  if(commandLineParams.exists("ANNOTATION_MODEL")) {
    string aux = commandLineParams.getValue("ANNOTATION_MODEL");
    int found = aux.find_last_of("/\\");
    if(found != string::npos) {
      data.annotationModel = aux.substr(found + 1);
      data.annotationModelDir = aux.substr(0, found);
    } else {
      data.annotationModelDir = ".";
      data.annotationModel = aux;
    }
  }


  ////////////////////////////////////////////////////////////////////////////////
  //  Input data section - files needed to generate the tables

  // abruijn graph
  if(commandLineParams.exists("FILE_COMP")) {
    data.filenameAbruijn = commandLineParams.getValue("FILE_COMP");
    data.absoluteFilenameAbruijn = true;
  }

  if(commandLineParams.exists("FILE_SEQS")) {
    data.filenameSpsSeqs = commandLineParams.getValue("FILE_SEQS");
    data.absoluteFilenameSpsSeqs = true;
  }

  if(commandLineParams.exists("FILE_MP2")) {
    data.filenameContigsMp = commandLineParams.getValue("FILE_MP2");
    data.absoluteFilenameContigsMp = true;
  }

  if(commandLineParams.exists("FILE_MIDX2")) {
    data.filenameContigsMidx = commandLineParams.getValue("FILE_MIDX2");
    data.absoluteFilenameContigsMidx = true;
  }

  if(commandLineParams.exists("FILE_REFMP")) {
    data.filenameHomglueRefMp = commandLineParams.getValue("FILE_REFMP");
    data.absoluteFilenameHomglueRefMp = true;
  }

  if(commandLineParams.exists("FILE_REFMIDX")) {
    data.filenameHomglueRefMidx = commandLineParams.getValue("FILE_REFMIDX");
    data.absoluteFilenameHomglueRefMidx = true;
  }

  if(commandLineParams.exists("FILE_STARS")) {
    data.filenameStarSpectra = commandLineParams.getValue("FILE_STARS");
    data.absoluteFilenameStarSpectra = true;
  }

  if(commandLineParams.exists("FILE_MS")) {
    data.filenameConsensusSpectra = commandLineParams.getValue("FILE_MS");
    data.absoluteFilenameConsensusSpectra = true;
  }

  if(commandLineParams.exists("FILE_INDEX")) {
    data.filenameInputSpectraList = commandLineParams.getValue("FILE_INDEX");
    data.absoluteFilenameInputSpectraList = true;
  }

  if(commandLineParams.exists("FILE_CLUSTERMS")) {
    data.filenameClusterMS = commandLineParams.getValue("FILE_CLUSTERMS");
    data.absoluteFilenameClusterMS = true;
  }

  if(commandLineParams.exists("FILE_CLUSTERSCAN")) {
    data.filenameScanFiles = commandLineParams.getValue("FILE_CLUSTERSCAN");
    data.absoluteFilenameScanFiles = true;
  }

  if(commandLineParams.exists("FILE_FASTA")) {
    data.filenameProteins = commandLineParams.getValue("FILE_FASTA");
    data.absoluteFilenameProteins = true;
  }

  if(commandLineParams.exists("MATCHED_CONTIGS")) {
    data.filenameContigSpectra = commandLineParams.getValue("MATCHED_CONTIGS");
    data.absoluteFilenameContigSpectra = true;
  }

  if(commandLineParams.exists("MATCHED_CONTIGS_IDX")) {
    data.filenameContigIndices = commandLineParams.getValue("MATCHED_CONTIGS_IDX");
    data.absoluteFilenameContigIndices = true;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //  File save section

  // Output directory
  if (commandLineParams.exists("REPORT_DIR"))
    data.outDir = commandLineParams.getValue("REPORT_DIR");


  ////////////////////////////////////////////////////////////////////////////////
  //  Spectrum related parameters

  // spectrumInfo
  if (commandLineParams.exists("spectruminfo"))
    ;
    //ReportInterfacetrum.setSpectrumInfo(commandLineParams.getValue("spectruminfo"));

  // peak mass tolerance
  if (commandLineParams.exists("TOLERANCE_PEAK"))
    data.peakMassTol = commandLineParams.getValueFloat("TOLERANCE_PEAK");

  // parent mass tolerance
  if (commandLineParams.exists("TOLERANCE_PM"))
    data.peakMassTol = commandLineParams.getValueFloat("TOLERANCE_PM");



  ////////////////////////////////////////////////////////////////////////////////
  //  Report related parameters

  // Executables directory
  if (commandLineParams.exists("JOB"))
    data.job = commandLineParams.getValue("JOB");

  // Executables directory
  if (commandLineParams.exists("USER"))
    data.user = commandLineParams.getValue("USER");


  ////////////////////////////////////////////////////////////////////////////////
  //  Image related parameters

  // Title
  if (commandLineParams.exists("TITLE"))
    ;
    //ReportInterfacetrum.setTitle(commandLineParams.getValue("TITLE"));

  if (commandLineParams.exists("NO_MSMS_IMAGES"))
    data.displayLevel = 3;


  ////////////////////////////////////////////////////////////////////////////////
  //  Action commands

  // the debug flag
  if (commandLineParams.exists("debug")) {
    dump_abruijn(data);
    dump_binArray(data);
    return 0;
  }

  // Build text tables
  if (commandLineParams.exists("buildTables")) {
    showHelp = false;
    buildTables(data);
  }

  // test for both dynamic and static HTML reports -- no can do
  if (commandLineParams.exists("buildhtml") && commandLineParams.exists("REPORT_DYNAMIC")) {
    return error("static and dynamic report types cannot be specified simultaneously.");
  }

  // Build HTML reports
  if (commandLineParams.exists("buildhtml")) {
    showHelp = false;
    if(buildDirectoryPath(data.outDir) == OK)
      generateReportHtml(data);
  }

  // Build HTML reports
  if (commandLineParams.exists("REPORT_DYNAMIC")) {
    showHelp = false;

    data.targetProjectDir = commandLineParams.getValue("PROJECT_DIR_SERVER");
    if(buildDirectoryPath(data.outDir) == OK)
      generateReportHtmlDynamic(data);
  }


  if(showHelp)
    help(cout);


  // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: spsReports [OPTION]\n";
  sout << "Options:\n";
  sout << '\n';
  sout << "  --build-tables              Build report data tables\n";
  sout << "  --build-html                Build static html report pages\n";
  sout << "  --build-html-dynamic        Build dynamic html report entry page\n";
  sout << '\n';
  sout << "  --exe-dir PATH              Executables directory (defaults to spsReports exec dir)\n";
  sout << "  --project-dir PATH          Project directory (defaults to current directory)\n";
  sout << "  --tables-dir PATH           Table files directory (defaults to project directory)\n";
  sout << "  --project-dir-server        Directory for project after relocation. Used by dynamic reports. Defaults to --project-dir";
  sout << '\n';
  sout << "  --tables-header FILENAME    Header table filename\n";
  sout << "  --tables-protein FILENAME   Proteins table filename\n";
  sout << "  --tables-contig FILENAME    Contigs table filename\n";
  sout << "  --tables-cluster FILENAME   Cluster consensus table filename\n";
  sout << "  --tables-spectra FILENAME   Input spectra table filename\n";
  sout << '\n';
  sout << "  --aa FILE                   Amino acids file (txt format)\n";
  sout << "  --annotation-model FILE     Annotation model file (defaults to './model_cid.txt')\n";
  sout << '\n';
  sout << "  --peakmasstol NUMBER        Amino acid peak mass tolerance (defaults to 0.45)\n";
  sout << "  --parentmasstol NUMBER      Parent mass tolerance\n";
  sout << "  --shift-value               Specify and apply mass shift\n";
  sout << '\n';
  sout << "  --abruijn FILE              Abruijn graph\n";
  sout << "  --sps-seqs FILE             \n";
  sout << "  --contigs-mp FILE           \n";
  sout << "  --contigs-midx FILE         \n";
  sout << "  --homglue-ref-midx FILE     \n";
  sout << "  --stars FILE                Star spectra file\n";
  sout << "  --consensus-spectra FILE    Consensus spectra file\n";
  sout << "  --input-spectra-list FILE   Input spectra list file\n";
  sout << "  --cluster-ms FILE           \n";
  sout << "  --cluster-scan FILE         \n";
  sout << "  --proteins FILE             Fasta file containg protein information\n";
  sout << '\n';
  sout << "  --server URL                Path to server cgi-bin directory\n";
  sout << "  --user STRING               User name\n";
  sout << "  --job STRING                job name\n";
  sout << "  --title STRING              Specifies title\n";
  sout << "  --cells-per-line NUMBER     Specifies number of aa cells per line in protein coverage view\n";
  sout << "  --no-msms-images            Specifies that no MS/MS images are generated\n";
  sout << '\n';
  sout << "  --outdir PATH               Output directory (defaults to <--project-dir>/report)\n";
  sout << '\n';
  sout << "  --verbose                   Display progress information\n";
  sout << '\n';
  sout << "  --p FILE                    Read parameters from file\n";
  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
void ReportInterface::verboseOutput(ostream &os, const char *str, bool nl)
{
  if(!m_verbose) return;
  os << str;
  if(nl)
    os << endl;
  else
    os.flush();
}
///////////////////////////////////////////////////////////////////////////////
void ReportInterface::verboseOutput(ostream &os, const char *str, const char *val, bool nl)
{
  if(!m_verbose) return;
  os << str << val;
  if(nl)
    os << endl;
  else
    os.flush();
}
///////////////////////////////////////////////////////////////////////////////
void ReportInterface::verboseOutput(ostream &os, const char *str, const char *val, const char *term, bool nl)
{
  if(!m_verbose) return;
  os << str << val << term;
  if(nl)
    os << endl;
  else
    os.flush();
}
void ReportInterface::verboseOutput(ostream &os, const double &v, bool nl)
{
  if(!m_verbose) return;
  os << v;
  if(nl)
    os << endl;
  else
    os.flush();
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::version(ostream & sout)
{
  sout << "spsReports 3.0." << XSTR(SPS_VERSION) << '\n';
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ReportInterface::error(const string & a)
{
  cerr << "spsReports: invalid or missing option: " << a << endl;

  cerr << "Type 'spsReports --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
