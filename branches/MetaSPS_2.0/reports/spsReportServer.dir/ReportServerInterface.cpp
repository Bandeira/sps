///////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <fstream>

#include "ReportServerInterface.h"

#include "ParameterList.h"
#include "CommandLineParser.h"

#include "Logger.h"

#include "db_fasta.h"

#include "ReportTableHeader.h"
#include "ReportTableProtein.h"
#include "ReportTableProteinCoverage.h"
#include "ReportTableContig.h"
#include "ReportTableClusterConsensus.h"
#include "ReportTableInputSpectra.h"

#include "copyright.h"

///////////////////////////////////////////////////////////////////////////////
using namespace specnets;
using namespace std;

namespace spsReports {
///////////////////////////////////////////////////////////////////////////////
ReportServerInterface::ReportServerInterface() :
   m_tableNameHeader("tableHeader.txt"),
   m_tableNameProtein("tableProtein.txt"),
   m_tableNameProteinCoverage("tableProteinCoverage.txt"),
   m_tableNameContig("tableContig.txt"),
   m_tableNameCluster("tableCluster.txt"),
   m_tableNameSpectra("tableSpectra.txt")
{
}
///////////////////////////////////////////////////////////////////////////////
ReportServerInterface::~ReportServerInterface()
{
}
///////////////////////////////////////////////////////////////////////////////
string ReportServerInterface::composeFileName(const string &projectDir, const string &fileName)
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
int ReportServerInterface::parseOptions(int argc, char **argv)
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

  listOptions.push_back(CommandLineParser::Option("-help",                    "help",                   false));
  listOptions.push_back(CommandLineParser::Option("-version",                 "VERSION",                false));

  // Data tables names and project location
  listOptions.push_back(CommandLineParser::Option("-project-dir",             "PROJECT_DIR",            true));
  listOptions.push_back(CommandLineParser::Option("-tables-dir",              "TABLES_DIR",             true));
  listOptions.push_back(CommandLineParser::Option("-table-header",            "TABLE_HEADER",           true));
  listOptions.push_back(CommandLineParser::Option("-table-protein",           "TABLE_PROTEIN",          true));
  listOptions.push_back(CommandLineParser::Option("-table-protein-coverage",  "TABLE_PROTEIN_COVERAGE", true));
  listOptions.push_back(CommandLineParser::Option("-table-contig",            "TABLE_CONTIG",           true));
  listOptions.push_back(CommandLineParser::Option("-table-cluster",           "TABLE_CLUSTER",          true));
  listOptions.push_back(CommandLineParser::Option("-table-spectra",           "TABLE_SPECTRA",          true));

  // data access & location
  listOptions.push_back(CommandLineParser::Option("-table",                   "TABLE",                  true));
  listOptions.push_back(CommandLineParser::Option("-filter-field",            "FILTER_FIELD",           true));
  listOptions.push_back(CommandLineParser::Option("-filter-data",             "FILTER_DATA",            true));
  listOptions.push_back(CommandLineParser::Option("-update-field",            "UPDATE_FIELD",           true));
  listOptions.push_back(CommandLineParser::Option("-update-data",             "UPDATE_DATA",            true));
  listOptions.push_back(CommandLineParser::Option("-clear-data" ,             "CLEAR_DATA",            false));

  // Sequence to update in fasta file
  listOptions.push_back(CommandLineParser::Option("-sequence",                "SEQ",                    true));
  listOptions.push_back(CommandLineParser::Option("-description",             "DESC",                   true));
  listOptions.push_back(CommandLineParser::Option("-ID",                      "ID",                     true));
  listOptions.push_back(CommandLineParser::Option("-filename",                "FILENAME",               true));

  // parameter file
  listOptions.push_back(CommandLineParser::Option("-p",                       "PARAMETER_FILE",         true));

  // Get data command
  listOptions.push_back(CommandLineParser::Option("-get",                     "GET",                    false));

  // Update table command
  listOptions.push_back(CommandLineParser::Option("-update",                  "UPDATE",                 false));

  // Launch command
  listOptions.push_back(CommandLineParser::Option("-launch",                  "LAUNCH",                 true));

  // Status command
  listOptions.push_back(CommandLineParser::Option("-status",                  "STATUS",                 false));

  // Request ID. Used to fool "smart" browsers that cache AJAX requests
  listOptions.push_back(CommandLineParser::Option("-request-id",              "REQ_ID",           true));

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
int ReportServerInterface::processOptions(ParameterList &ip)
{

  ////////////////////////////////////////////////////////////////////////////////
  // help message control
  bool showHelp = true;
  bool missing  = false;
  bool bTable = false, bFilterField = false, bFilterData = false, bUpdateField = false, bUpdateData = false;

  ReportData data;

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
  // "help" prints help and exits
  if (commandLineParams.exists("VERSION"))
    return version(cout);

  ////////////////////////////////////////////////////////////////////////////////
  // the same for "version"
  if (commandLineParams.exists("help"))
    return help(cout);


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

  // Tables directory
  if (commandLineParams.exists("TABLES_DIR"))
     data.tablesDir = commandLineParams.getValue("TABLES_DIR");


  // table file names
  if (commandLineParams.exists("TABLE_HEADER"))
    m_tableNameHeader = commandLineParams.getValue("TABLE_HEADER");

  if (commandLineParams.exists("TABLE_PROTEIN"))
    m_tableNameProtein = commandLineParams.getValue("TABLE_PROTEIN");

  if (commandLineParams.exists("TABLE_PROTEIN_COVERAGE"))
    m_tableNameProteinCoverage = commandLineParams.getValue("TABLE_PROTEIN_COVERAGE");

  if (commandLineParams.exists("TABLE_CONTIG"))
    m_tableNameContig = commandLineParams.getValue("TABLE_CONTIG");

  if (commandLineParams.exists("TABLE_CLUSTER"))
    m_tableNameCluster = commandLineParams.getValue("TABLE_CLUSTER");

  if (commandLineParams.exists("TABLE_SPECTRA"))
    m_tableNameSpectra = commandLineParams.getValue("TABLE_SPECTRA");


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
  //  Data commands

  // get table ID
  if (commandLineParams.exists("TABLE")) {
    data.tableID = getInt(commandLineParams.getValue("TABLE").c_str());
    bTable = true;
  }

  if (commandLineParams.exists("FILTER_FIELD")) {
    data.filterField = getInt(commandLineParams.getValue("FILTER_FIELD").c_str());
    bFilterField = true;
  }

  if (commandLineParams.exists("FILTER_DATA")) {
    data.filterData = commandLineParams.getValue("FILTER_DATA");
    bFilterData = true;
  }

  if (commandLineParams.exists("UPDATE_FIELD")) {
    data.updateField = getInt(commandLineParams.getValue("UPDATE_FIELD").c_str());
    bUpdateField = true;
  }

  if (commandLineParams.exists("UPDATE_DATA")) {
    data.updateData = commandLineParams.getValue("UPDATE_DATA");
    bUpdateData = true;
  }

  if (commandLineParams.exists("CLEAR_DATA")) {
    data.updateData = "";
    bUpdateData = true;
  }


  ////////////////////////////////////////////////////////////////////////////////
  //  Action commands

  // Get table section
  if (commandLineParams.exists("GET")) {
    // disable help command
    showHelp = false;
    // test for all needed parameters present
    if(!bTable)
      return error("--table");
    if(!bFilterField && bFilterData)
      return error("--filter-data needs ---filter-field");
    if(bFilterField && !bFilterData)
      return error("--filter-field needs --filter-data");
    // call command
    getTableData(data);
  }

  // Update text tables
  if (commandLineParams.exists("UPDATE")) {
    // disable help command
    showHelp = false;
    // test for all needed parameters present
    if(!bTable)
      return error("--table");
    if(!bFilterField)
      return error("--filter-field");
    if(!bFilterData)
      return error("--filter-data");
    if(!bUpdateField)
      return error("--update-field");
    if(!bUpdateData)
      return error("--update-data");
    // call command
    updateTable(data);
  }

  // Get status file contents
  if (commandLineParams.exists("STATUS")) {
    // disable help command
    showHelp = false;

    // get status file
    if (!commandLineParams.exists("FILENAME")) {
      cout << "Status filename not provided." << endl;
      return 0;
    }

    // get the project directory
    string prjDir = ".";
    if (commandLineParams.exists("PROJECT_DIR"))
      prjDir  = commandLineParams.getValue("PROJECT_DIR");

    // compose filename
    string fn = composeFileName(prjDir, commandLineParams.getValue("FILENAME"));
    // open the file and read contents
    ifstream f;
    f.open(fn.c_str());
    string output;
    if (f.is_open()) {
      f >> output;
    }
    // close the file
    f.close();
    // output the contents
    cout << output;
  }

  // Lauch project from given point
  if (commandLineParams.exists("LAUNCH")) {
    // disable help command
    showHelp = false;

    // get the project directory
    string prjDir = ".";
    if (commandLineParams.exists("PROJECT_DIR"))
      prjDir  = commandLineParams.getValue("PROJECT_DIR");

    // get fasta file
    if (!commandLineParams.exists("FILENAME")) {
      cout << "Fasta filename not provided." << endl;
      return 0;
    }
    string fastaFilename = composeFileName(prjDir, commandLineParams.getValue("FILENAME"));

    // get sequence
    if (!commandLineParams.exists("ID")) {
      cout << "Protein ID not provided.";
      return 0;
    }
    string id = commandLineParams.getValue("ID");

    // get sequence
    if (!commandLineParams.exists("DESC")) {
      cout << "Protein description not provided.";
      return 0;
    }
    string desc = commandLineParams.getValue("DESC");

    // get sequence
    if (!commandLineParams.exists("SEQ")) {
      cout << "Protein sequence not provided.";
      return 0;
    }
    string seq = commandLineParams.getValue("SEQ");

    // add sequence to fasta
    specnets::DB_fasta  fasta;
    if(fasta.Load(fastaFilename.c_str()) == 0) {
      cout << "Error opening " << fastaFilename << endl;
      //return 0;
    }
    fasta.Add((char *)id.c_str(), (char *)desc.c_str(), (char *)seq.c_str());
    string saveFN = composeFileName(prjDir, "relaunch.protid.fasta");
    fasta.Save(saveFN.c_str());

    string toGo = "cd ";
    toGo += prjDir;
    toGo += " ; . ";

    // get program to launch name
    string prg = commandLineParams.getValue("LAUNCH");
    // compose location file name
    toGo += composeFileName(prjDir, prg);
    toGo += " &";

    // execute the program
    system(toGo.c_str());
  }

  // show help if no command specified
  if(showHelp)
    help(cout);

  // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int ReportServerInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: spsReports [OPTION]\n";
  sout << "Options:\n";
  sout << '\n';
  sout << "  --get                                Get data\n";
  sout << "  --update                             Update table cell\n";
  sout << "  --launch EXECUTABLE                  Launch program\n";
  sout << "  --status                             Get status from status file\n";
  sout << '\n';
  sout << "  --table ID                           Table ID in both get and update operations\n";
  sout << "  --filter-field ID                    Field index to apply the data filter\n";
  sout << "  --filter-data STRING                 Filter data to apply on the 'update field'\n";
  sout << "  --update-field ID                    Field index to be updated\n";
  sout << "  --update-data STRING                 Data to be placed in the update field\n";
  sout << '\n';
  sout << "  --exe-dir PATH                       Executables directory (defaults to spsReports exec dir)\n";
  sout << "  --project-dir PATH                   Project directory (defaults to current directory)\n";
  sout << "  --tables-dir PATH                    Table files directory (defaults to project directory)\n";
  sout << '\n';
  sout << "  --tables-header FILENAME             Header table filename\n";
  sout << "  --tables-protein FILENAME            Proteins table filename\n";
  sout << "  --tables-protein-coverage FILENAME   Proteins table filename\n";
  sout << "  --tables-contig FILENAME             Contigs table filename\n";
  sout << "  --tables-cluster FILENAME            Cluster consensus table filename\n";
  sout << "  --tables-spectra FILENAME            Input spectra table filename\n";
  sout << '\n';
  sout << "  --sequence AAs                       Protein sequence to add to the fasta file\n";
  sout << "  --description TEXT                   Description of protein to add\n";
  sout << "  --ID id                              ID of protein to add\n";
  sout << '\n';
  sout << "  --filename FILENAME                  Filename of file to open.\n";
  sout << '\n';
  sout << "  --p FILE                             Read parameters from file\n";
  sout << '\n';
  sout << "  --help                               Display this help and exit\n";
  sout << "  --version                            Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ReportServerInterface::version(ostream & sout)
{
  sout << "spsReportsServer 3.0." << XSTR(SPS_VERSION) << '\n';
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ReportServerInterface::error(const string & a)
{
  cerr << "spsReportsServer: invalid or missing option: " << a << endl;

  cerr << "Type 'spsReportsServer --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
ReportTableBase * ReportServerInterface::getTableObject(ReportData &data)
{
  ReportTableBase *tab = NULL;

  switch(data.tableID) {
  case 1:
    tab = new ReportTableHeader(data.tablesDir, m_tableNameHeader);
    break;
  case 2:
    tab = new ReportTableProtein(data.tablesDir, m_tableNameProtein);
    break;
  case 3:
    tab = new ReportTableContig(data.tablesDir, m_tableNameContig);
    break;
  case 4:
    tab = new ReportTableClusterConsensus(data.tablesDir, m_tableNameCluster);
    break;
  case 5:
    tab = new ReportTableInputSpectra(data.tablesDir, m_tableNameSpectra);
    break;
  case 6:
    tab = new ReportTableProteinCoverage(data.tablesDir, m_tableNameProteinCoverage);
    break;
  }
  return tab;
}
///////////////////////////////////////////////////////////////////////////////
int ReportServerInterface::getTableData(ReportData &data)
{
  ReportTableBase *tab = getTableObject(data);

  ReportTableData td;
  td.filterCol = data.filterField;
  td.filterText = data.filterData;
  td.updateCol = data.updateField;
  td.updateText = data.updateData;

  if(tab)
    tab->getData(cout, td, ';', '\n');
}
///////////////////////////////////////////////////////////////////////////////
int ReportServerInterface::updateTable(ReportData &data)
{
  ReportTableBase *tab = getTableObject(data);

  ReportTableData td;
  td.filterCol = data.filterField;
  td.filterText = data.filterData;
  td.updateCol = data.updateField;
  //td.updateText = data.updateData;

  // make sure no weird characteres are int the sequence
  for(int i = 0 ; i <  data.updateData.length() ; i++) {
    if(data.updateData[i] == '"')
      continue;
    td.updateText.push_back(data.updateData[i]);
  }


  if(tab)
    tab->update(td);
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
