///////////////////////////////////////////////////////////////////////////////
#include "SpecServerInterface.h"

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "mzxml.h"

#include "Timer.h"

#include "copyright.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
SpecServerInterface::SpecServerInterface()
{
}

SpecServerInterface::~SpecServerInterface()
{
}
///////////////////////////////////////////////////////////////////////////////
int SpecServerInterface::processOptions(int argc, char **argv)
{
  //--------------------------------------------------------------------------------------------
  // Initialize directories used
  //--------------------------------------------------------------------------------------------
  // Get the execultable directory. Unless specified, renderer directory and font directory is expected
  string str = argv[0];

  //convert path
  size_t found;
  found = str.find_last_of("/\\");
  string aux;
  aux = '.';
  if(found != string::npos)
    aux = str.substr(0, found);
  string exePath = aux;

  exePath = replaceAll(exePath, "\\", "/");


  //--------------------------------------------------------------------------------------------
  // Parse the command line parameters
  //--------------------------------------------------------------------------------------------
  vector<CommandLineParser::Option> listOptions;

  //  c++ instruction                            cmd line option    parameter name    takes value?

  listOptions.push_back(CommandLineParser::Option("-help",            "help",             false));
  listOptions.push_back(CommandLineParser::Option("-version",         "VERSION",          false));
  listOptions.push_back(CommandLineParser::Option("-annot",           "ANNOTATIONS",      true));
//  listOptions.push_back(CommandLineParser::Option("-verbose",         "VERBOSE",          false));

  listOptions.push_back(CommandLineParser::Option("-pklbin",          "PKLBIN",           true));
  listOptions.push_back(CommandLineParser::Option("-mgf",             "MGF",              true));
  listOptions.push_back(CommandLineParser::Option("-pkl",             "PKL",              true));
  listOptions.push_back(CommandLineParser::Option("-mzxml",           "MZXML",            true));

  listOptions.push_back(CommandLineParser::Option("-spectrum",        "SPECTRUM",         false));
  listOptions.push_back(CommandLineParser::Option("-spectrumscan",    "SPECTRUMSCAN",     false));

  listOptions.push_back(CommandLineParser::Option("-request-id",      "REQ_ID",           true));

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

  return processOptions(commandLineParams);
}
////////////////////////////////////////////////////////////////////////////////
int SpecServerInterface::processOptions(ParameterList & ip)
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

  // pklbin file to load
  if (commandLineParams.exists("PKLBIN")) {
    string fileName = commandLineParams.getValue("PKLBIN").c_str();
    fileName = replaceAll(fileName, "\\", "/");
    fileLoaded = specSet.loadPklBin(fileName.c_str());
    if(!fileLoaded) {
      stringstream err; err << "Error loading pklbin file: " << fileName;
      return error(err.str());
    }
  }

  // mgf file to load
  if (commandLineParams.exists("MGF")) {
    string fileName = commandLineParams.getValue("MGF").c_str();
    fileName = replaceAll(fileName, "\\", "/");
    fileLoaded = specSet.LoadSpecSet_mgf(fileName.c_str());
    if(!fileLoaded) {
      stringstream err; err << "Error loading mgf file: " << fileName;
      return error(err.str());
    }
  }

  // pkl file to load
  if (commandLineParams.exists("PKL")) {
    string fileName = commandLineParams.getValue("PKL").c_str();
    fileName = replaceAll(fileName, "\\", "/");
    fileLoaded = specSet.LoadSpecSet_pkl(fileName.c_str());
    if(!fileLoaded) {
      stringstream err; err << "Error loading pkl file: " << fileName;
      return error(err.str());
    }
  }

  int ret;


  // mzxml file to load
  if (commandLineParams.exists("MZXML")) {
    // get the filename
    string mzxmlFileName = commandLineParams.getValue("MZXML");
    mzxmlFileName = replaceAll(mzxmlFileName, "\\", "/");
    // auxilizary vector needed
    vector<short> msLevel;
    // load method depends on the presence of spectrumscan specification
    int ret = LoadMzxml( (char * const)(mzxmlFileName.c_str()), specSet, & msLevel, 2);
    if(! ret) {
      stringstream err; err << "Error loading mzxml file.";
      return error(err.str());
    }
    fileLoaded = true;
  }

  // Check if a file was loaded
  if(!fileLoaded)
    return error("No file loaded.");


  ////////////////////////////////////////////////////////////////////////////////
  //  Spectrum related parameters

  // spectrum index used
  if (commandLineParams.exists("SPECTRUM")) {
    cout << specSet.size() << ';';
  }


  // Spectrum scan processing
  if(commandLineParams.exists("SPECTRUMSCAN")) {
    for(int i = 0 ; i < specSet.size() ; i++) {
      cout << specSet[i].scan << ';';
    }
  }

  cout << endl;
 
  // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int SpecServerInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: specplot [OPTION]\n";
  sout << "Options:\n";
  sout << "  --pklbin FILE               Spectrum data (pklbin format)\n";
  sout << "  --mgf FILE                  Spectrum data (mgf format)\n";
  sout << "  --pkl FILE                  Spectrum data (pkl format)\n";
  sout << "  --mzxml FILE                Spectrum data (mzxml format)\n";
  sout << '\n';
  sout << "  --spectrum                  Return number of spectra in file\n";
  sout << "  --spectrumscan              Return all Spectrum Scan numbers in file\n";
  sout << '\n';
  sout << "  --p FILE                    Read parameters from file\n";
  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecServerInterface::version(ostream & sout)
{
  sout << "SpecServer 1.0." << XSTR(SPS_VERSION) << '\n';
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int SpecServerInterface::error(const string & a)
{
  cerr << "SpecServer: invalid or missing option: " << a << endl;

  cerr << "Type 'SpecServer --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
