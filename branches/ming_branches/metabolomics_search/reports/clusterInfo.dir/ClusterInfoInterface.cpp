///////////////////////////////////////////////////////////////////////////////
#include "ClusterInfoInterface.h"


#include "CommandLineParser.h"
#include "ParameterList.h"
#include "../utils.h"
#include "Logger.h"

#include "copyright.h"

#define DEFAULT_BIN_FILES_FILENAME  "bin_files.txt"
#define DEFAULT_INPUT_FILES_LIST  "spectra/input_index.txt"
#define DEFAULT_AUTFILE_NAME "ClusterInfo.csv"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
ClusterInfoInterface::ClusterInfoInterface() : outdir("."), outFileName(DEFAULT_AUTFILE_NAME),
  binFilesNames(DEFAULT_BIN_FILES_FILENAME),
  inputFilesNames(DEFAULT_INPUT_FILES_LIST),
  projectdir(".")
{
}

ClusterInfoInterface::~ClusterInfoInterface()
{
}

///////////////////////////////////////////////////////////////////////////////
void  ClusterInfoInterface::loadStringVector(string &inputFilesNames, vector<std::string> &inputFiles)
{
  std::vector<std::string> file;
  std::string line;
  file.clear();
  std::ifstream infile (inputFilesNames.c_str(), std::ios_base::in);
  while (getline(infile, line, '\n')) {
    inputFiles.push_back (line);
  }

  cout << "for: " << inputFilesNames << " got:" << endl;
  for(int i = 0 ; i < inputFiles.size() ; i++)
    cout << inputFiles[i] << endl;
}
///////////////////////////////////////////////////////////////////////////////
int ClusterInfoInterface::processOptions(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------
  // Initialize directories used
  //--------------------------------------------------------------------------------------------


  //--------------------------------------------------------------------------------------------
  // Parse the command line parameters
  //--------------------------------------------------------------------------------------------
  vector<CommandLineParser::Option> listOptions;

  //  c++ instruction                            cmd line option    parameter name    takes value?

  listOptions.push_back(CommandLineParser::Option("-help",            "help",             false));
  listOptions.push_back(CommandLineParser::Option("-version",         "VERSION",          false));


  listOptions.push_back(CommandLineParser::Option("-project-dir",     "PROJECT_DIR",      true));
  listOptions.push_back(CommandLineParser::Option("-inputspectra",    "INPUT_SPECTRA",    true));
  listOptions.push_back(CommandLineParser::Option("-inputbin",        "INPUT_BIN",        true));

  listOptions.push_back(CommandLineParser::Option("-outdir",          "OUTDIR",           true));
  listOptions.push_back(CommandLineParser::Option("-outfile",         "OUTFILE",          true));
  // parameter file
  listOptions.push_back(CommandLineParser::Option("-p",               "PARAMETER_FILE",   true));


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

  ////////////////////////////////////////////////////////////////////////////////
  // Parameter file
  ////////////////////////////////////////////////////////////////////////////////

  if (commandLineParams.exists("PARAMETER_FILE")) {

    string parameterFilename = commandLineParams.getValue("PARAMETER_FILE");

    ParameterList ip;
    ip.readFromFile(parameterFilename);
    // Combine the command line parameters to the file ones
    //   Command line parameters take precedence (hence the overwrite flag not set)
    commandLineParams.addList(ip, false);
  }

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

  // A file must be loaded
  bool fileLoaded = false;

  // input spectra file to load
  if (commandLineParams.exists("INPUT_SPECTRA")) {
    inputFilesNames = commandLineParams.getValue("INPUT_SPECTRA").c_str();
  }

  // input spectra file to load
  if (commandLineParams.exists("INPUT_BIN")) {
    binFilesNames = commandLineParams.getValue("INPUT_BIN").c_str();
  }

  // cluster binary file to load
  if (commandLineParams.exists("PROJECT_DIR")) {
    projectdir = commandLineParams.getValue("PROJECT_DIR").c_str();
  }

  string aux = projectdir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  aux += inputFilesNames;

  loadStringVector(aux, inputFiles);


  aux = projectdir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  aux += binFilesNames;

  loadStringVector(aux, scanFileNames);


  for(int i = 0 ; i < scanFileNames.size() ; i++) {
    vector<vector<int> > aux;
    Load_binArray(scanFileNames[i].c_str(), aux);
    scanNumbers.push_back(aux);
  }

  ////////////////////////////////////////////////////////////////////////////////
  //  File save section

  // Output directory
  if (commandLineParams.exists("OUTDIR"))
    outdir = commandLineParams.getValue("OUTDIR");

  // Output file name.
  if (commandLineParams.exists("OUTFILE"))
    outFileName = commandLineParams.getValue("OUTFILE");

  // compose output filename
  aux = outdir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  aux += outFileName;

  fileLoaded = clusterData.load(projectdir.c_str());

  //cout << "# of loaded clusters: " << clusterData.data.size() << endl;

  clusterData.writeCsv(&inputFiles, &scanNumbers, aux);

   // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int ClusterInfoInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: clusterinfo [OPTION]\n";
  sout << "Options:\n";
  sout << "  --project-dir DIRECTORY     Project directory\n";
  sout << "  --clusterbin FILE           Clustering binary data file\n";
  sout << "  --inputspectra FILE         Input spectra file list (txt format, defaults to spectra/bin_files.txt)\n";
  sout << '\n';
  sout << "  --outdir PATH               Output directory (defaults to .)\n";
  sout << "  --outfile FILE              Output filename (defaults to ClusterInfo.csv)\n";
  sout << '\n';
  sout << "  --p FILE                    Read parameters from file\n";
  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ClusterInfoInterface::version(ostream & sout)
{
  sout << "clusterinfo 1.0." << XSTR(SPS_VERSION) << '\n';
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int ClusterInfoInterface::error(const string & a)
{
  cerr << "clusterinfo: invalid or missing option: " << a << endl;

  cerr << "Type 'clusterinfo --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
