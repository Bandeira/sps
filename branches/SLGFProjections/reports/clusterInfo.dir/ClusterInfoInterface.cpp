///////////////////////////////////////////////////////////////////////////////
#include "ClusterInfoInterface.h"


#include "CommandLineParser.h"
#include "ParameterList.h"
#include "utils.h"
#include "Logger.h"

#include "copyright.h"

#define DEFAULT_BIN_FILES_FILENAME  "spectra/bin_files.txt"
#define DEFAULT_INPUT_FILES_LIST  "spectra/pklbin_files.txt"
#define DEFAULT_ORIGINAL_FILES_LIST  "spectra/input_index.txt"
#define DEFAULT_AUTFILE_NAME "ClusterInfo.csv"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
ClusterInfoInterface::ClusterInfoInterface() : outdir("."), outFileName(DEFAULT_AUTFILE_NAME),
  binFilesNames(DEFAULT_BIN_FILES_FILENAME),
  inputFilesNames(DEFAULT_INPUT_FILES_LIST),
  originalFilesNames(DEFAULT_ORIGINAL_FILES_LIST),
  projectdir(".")
{
}

ClusterInfoInterface::~ClusterInfoInterface()
{
}

///////////////////////////////////////////////////////////////////////////////
int  ClusterInfoInterface::loadStringVector(string &inputFilesNames, vector<std::string> &inputFiles)
{
  std::vector<std::string> file;
  std::string line;
  file.clear();
  std::ifstream infile (inputFilesNames.c_str(), std::ios_base::in | std::ios_base::binary);

  if(!infile) return 0;

  while (getline(infile, line, '\n')) {
    inputFiles.push_back (line);
  }

  cout << "for: " << inputFilesNames << " got:" << endl;
  for(int i = 0 ; i < inputFiles.size() ; i++)
    cout << inputFiles[i] << endl;

  return 1;
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
  listOptions.push_back(CommandLineParser::Option("-originalspectra", "ORIGINAL_SPECTRA", true));

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
  if (commandLineParams.exists("ORIGINAL_SPECTRA")) {
    originalFilesNames = commandLineParams.getValue("ORIGINAL_SPECTRA").c_str();
  }


  // input spectra file to load
  if (commandLineParams.exists("INPUT_BIN")) {
    binFilesNames = commandLineParams.getValue("INPUT_BIN").c_str();
  }

  // cluster binary file to load
  if (commandLineParams.exists("PROJECT_DIR")) {
    projectdir = commandLineParams.getValue("PROJECT_DIR").c_str();
  }


  // load input files names
  string aux = projectdir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  aux += inputFilesNames;

  if(!loadStringVector(aux, inputFiles)) {
    cerr << "Unable to load file: " << aux << endl;
    exit(0);
  }

  // load original files names
  aux = projectdir;
  if(aux[aux.length()-1] != '/')
    aux += '/';
  aux += originalFilesNames;

  if(!loadStringVector(aux, originalFiles)) {
    cerr << "Unable to load file: " << aux << endl;
    exit(0);
  }


  // read input data into memory structure
  for(int i = 0 ; i < inputFiles.size() ; i++) {
    SpecSet spec;
    //cout << inputFiles[i] << " loaded" << endl;
    if(spec.LoadSpecSet_pklbin(inputFiles[i].c_str()) <= 0) {
      cerr << "Unable to load file: " << inputFiles[i] << endl;
      return 0;
    }
    inputData.push_back(spec);
  }


  bool readFromPklbin = true;

  if(commandLineParams.exists("INPUT_BIN")) {

    aux = projectdir;
    if(aux[aux.length()-1] != '/')
      aux += '/';
    aux += binFilesNames;

    //cout << "chekcing for " << aux << endl;
    if(loadStringVector(aux, scanFileNames)) {
      //cout << "reading scan numbers from .bin" << endl;
      for(int i = 0 ; i < scanFileNames.size() ; i++) {
        vector<vector<int> > aux;
        Load_binArray(scanFileNames[i].c_str(), aux);
        //cout << aux.size() << " ; " << aux[0].size() << " ; " << aux[0][0] << " ; " << aux[0][1] << endl;
        scanNumbers.push_back(aux);
      }
      readFromPklbin = false;
    }
  }


  if(readFromPklbin) {
    //cout << "reading scan numbers from .pklbin" << endl;
    for(int i = 0 ; i < inputData.size() ; i++) {
      vector<vector<int> > aux;
      //cout << "Has " << spec.size() << " spectra" << endl;
      for(int j = 0 ; j < inputData[i].size() ; j++) {
        vector<int> aux2;
        aux2.push_back(inputData[i][j].scan);
        aux2.push_back(0);
        aux.push_back(aux2);
      }
      scanNumbers.push_back(aux);
    }
  }

/*
  if(readFromPklbin) {
    //cout << "reading scan numbers from .pklbin" << endl;
    for(int i = 0 ; i < inputFiles.size() ; i++) {
      vector<vector<int> > aux;
      SpecSet spec;
      //cout << inputFiles[i] << " loaded" << endl;
      if(spec.LoadSpecSet_pklbin(inputFiles[i].c_str()) <= 0) {
        cerr << "Unable to load file: " << inputFiles[i] << endl;
        return 0;
      }
      //cout << "Has " << spec.size() << " spectra" << endl;
      for(int j = 0 ; j < spec.size() ; j++) {
        vector<int> aux2;
        aux2.push_back(spec[j].scan);
        aux2.push_back(0);
        aux.push_back(aux2);
      }
      scanNumbers.push_back(aux);
    }
  } */

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

  cout << "# of loaded clusters: " << clusterData.data.size() << endl;

  clusterData.writeCsv(&originalFiles, &scanNumbers, aux, inputData);

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
  sout << "  --inputspectra FILE         Input spectra file list (txt format)\n";
  sout << "  --originalspectra FILE      Original spectra file list (txt format)\n";
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
