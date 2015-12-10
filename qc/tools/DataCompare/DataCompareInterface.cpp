///////////////////////////////////////////////////////////////////////////////
#include "DataCompareInterface.h"
#include "DataComparer.h"

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "Logger.h"
#include "mzxml.h"

#include "copyright.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
DataCompareInterface::DataCompareInterface()
{
}

DataCompareInterface::~DataCompareInterface()
{
}
///////////////////////////////////////////////////////////////////////////////
int DataCompareInterface::processOptions(int argc, char **argv)
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

#if  defined(__MINGW32__) || defined(__CYGWIN__)
  exePath = replaceAll(exePath, "\\", "/");
#endif

  //--------------------------------------------------------------------------------------------
  // Parse the command line parameters
  //--------------------------------------------------------------------------------------------
  vector<CommandLineParser::Option> listOptions;

  //  c++ instruction                            cmd line option    parameter name    takes value?

  listOptions.push_back(CommandLineParser::Option("-help",            "help",             false));
  listOptions.push_back(CommandLineParser::Option("-version",         "VERSION",          false));
//  listOptions.push_back(CommandLineParser::Option("-verbose",         "verbose",          false));

  listOptions.push_back(CommandLineParser::Option("-format",          "FORMAT",           true));
  listOptions.push_back(CommandLineParser::Option("-file1",           "FILE1",            true));
  listOptions.push_back(CommandLineParser::Option("-file2",           "FILE2",            true));
  listOptions.push_back(CommandLineParser::Option("-detail",          "DETAIL",           true));

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
//////////////////////////////////////////////////////////////////////////////////
int DataCompareInterface::processOptions(ParameterList &ip)
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
  // Boolean parameters


  ////////////////////////////////////////////////////////////////////////////////
  // File locations


  ////////////////////////////////////////////////////////////////////////////////
  // Object declaration
  
  DataCompare dataCompare;


  ////////////////////////////////////////////////////////////////////////////////
  // File load section


  // Detail level
  if (commandLineParams.exists("DETAIL")) {
    int aux = commandLineParams.getValueInt("DETAIL");
    dataCompare.setLevel(aux);
  }

  // Files format
  if (commandLineParams.exists("FORMAT")) {
    string aux = commandLineParams.getValue("FORMAT");
    dataCompare.setType(aux);
  } else {
    cerr << "Error: Files format not specified." << endl;
    return -1;
  }

  // File 1
  if (commandLineParams.exists("FILE1")) {
    string aux = commandLineParams.getValue("FILE1");
    int ret = dataCompare.addObject(aux, 0);
    if(!ret) {
      cerr << "Error: reading file: " << aux;
      return -1;
    }
  } else {
    cerr << "Error: missing input file file." << endl;
    return -1;
  }

  // File 2
  if (commandLineParams.exists("FILE2")) {
    string aux = commandLineParams.getValue("FILE2");
    int ret = dataCompare.addObject(aux, 1);
    if(!ret) {
      cerr << "Error: reading file: " << aux;
      return -1;
    }
  } else {
    cerr << "Error: missing input file file." << endl;
    return -1;
  }

  ////////////////////////////////////////////////////////////////////////////////
  //  File save section


  ////////////////////////////////////////////////////////////////////////////////
  //  Spectrum related parameters


  ////////////////////////////////////////////////////////////////////////////////
  //  Image related parameters

  //////////////////////////////////////////////////////////////////////////////
  // Other
  
  //////////////////////////////////////////////////////////////////////////////
  // Execute
  
  return dataCompare.exec();


  // return status ok
  return 0;
}
///////////////////////////////////////////////////////////////////////////////
int DataCompareInterface::help(ostream & sout)
{
  version(sout);

  sout << "Usage: dataCompare [OPTIONS]\n";
  sout << "Options:\n";

  sout << "  --format FORMAT             abruijn, bin (bin Array), pklbin, txt (report table file)\n";
  sout << "  --file1 FILENAME            1st file to compare\n";
  sout << "  --file2 FILENAME            2nd file to compare\n";
  sout << "  --detail VALUE              0 (default) to 4\n";
  sout << '\n';
  sout << "  --p FILE                    Read parameters from file\n";
  sout << '\n';
  sout << "  --help                      Display this help and exit\n";
  sout << "  --version                   Output version information and exit\n";
  sout << endl;

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int DataCompareInterface::version(ostream & sout)
{
  sout << "dataCompare 1.0." << XSTR(SPS_VERSION) << '\n';
  sout << COPYRIGHT1 << '\n';
  sout << COPYRIGHT2;
  sout << "\n\n";

  return 1;
}
///////////////////////////////////////////////////////////////////////////////
int DataCompareInterface::error(const string & a)
{
  cerr << "dataCompare: invalid or missing option: " << a << endl;

  cerr << "Type 'dataCompare --help' for more information." << endl << endl;

  return -1;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
