//
//  SpsModuleExecutor - stand alone execution driver for SpSModules
//
#include "CommandLineParser.h"
#include "ExecModuleFactoryLP.h"
#include "Logger.h"
#include "ParameterList.h"
#include "ParallelThreadedExecution.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"

#include "Logger.h"
#include "ExecFlr.h"

#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <unistd.h>
#include <time.h>

using namespace specnets;
using namespace std;

//-----------------------------------------------------------------------------
string getCurrentDirectory(void)
{
  string currentWorkindDir;
  currentWorkindDir = getcwd(NULL, 1024);
  return currentWorkindDir;

//  int dummy;
//  dummy = system("pwd > pwd.tmp");
//  ifstream ifs("pwd.tmp");
//  char buf[1024];
//  ifs.getline(buf, 1024);
//  return buf;
}
//-----------------------------------------------------------------------------

void addDefaultParameterValues(ParameterList &p)
{
  // Basic parameters
  p.addIfDoesntExist("TOLERANCE_PEAK", "0.4");

  // Grid parameters
  p.addIfDoesntExist("GRID_TYPE", "sge");
  p.addIfDoesntExist("GRID_NUMNODES", "0");
  p.addIfDoesntExist("GRID_NUMCPUS", "4");
  p.addIfDoesntExist("GRID_EXE_DIR", "../../../trunk/ExecFramework/FalseLocalizationRates");
  p.addIfDoesntExist("GRID_SGE_EXE_DIR", "");
  p.addIfDoesntExist("GRID_PARAMS", "-l h_vmem=1G");

  string currDir = getCurrentDirectory();
  p.addIfDoesntExist("PROJECT_DIR", currDir);
  p.addIfDoesntExist("RELATIVE_DIR", "1");
}

// -------------------------------------------------------------------------
string getCurrentTimeString(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  return asctime(timeinfo);
}
// -------------------------------------------------------------------------
bool performFlr(ParameterList & ip, bool gridExecutionFlag, bool resume)
{
  DEBUG_TRACE;
  ExecFlr moduleFlr(ip);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleFlr.validateParams(errorString);
  DEBUG_VAR(isValid);
  if (!isValid)
  {
    ERROR_MSG(errorString);
    return false;
  }

  isValid = moduleFlr.loadInputData();
  DEBUG_VAR(isValid);
  if (!isValid)
  {
    ERROR_MSG(errorString);
    return false;
  }

  bool returnStatus;

  if (!ip.exists("GRID_NUMNODES") or ip.getValueInt("GRID_NUMNODES") <= 0)
  {
    returnStatus = moduleFlr.invoke();
  }
  else
  {
    DEBUG_TRACE;
    int numNodes = ip.getValueInt("GRID_NUMNODES");

    string gridType = ip.getValue("GRID_TYPE");
    if (gridType == "pbs")
    {
      ParallelPbsExecution exec(&moduleFlr,
                                gridExecutionFlag,
                                !gridExecutionFlag,
                                resume);
      returnStatus = exec.invoke(numNodes);
    }
    else if (gridType == "sge")
    {
      ParallelSgeExecution exec(&moduleFlr,
                                gridExecutionFlag,
                                !gridExecutionFlag,
                                resume);
      returnStatus = exec.invoke(numNodes);
    }
  }
  // Test for return status
  DEBUG_VAR( returnStatus);
  if (!returnStatus)
  {
    ERROR_MSG("invoking moduleFlr exited in error.");
    return false;
  }

  DEBUG_TRACE;

  moduleFlr.saveOutputData();
  return true;
}
// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc < 2)
  {
    cerr << "Usage: main_flr parametersFileName [options]" << endl;
    return -1;
  }

  bool resumeFlag = false;
  bool gridExecutionFlag = false;

  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("t", "NUM_THREADS", true));
  listOptions.push_back(CommandLineParser::Option("lf", "LOG_FILE_NAME", true));
  listOptions.push_back(CommandLineParser::Option("ll", "LOG_LEVEL", true));

  // Command line parameters for ExecSpectraExtraction for CCMS Workflows
  listOptions.push_back(CommandLineParser::Option("ccms_output_library",
                                                  "OUTPUT_MGF",
                                                  true));
  listOptions.push_back(CommandLineParser::Option("ccms_input_library",
                                                  "EXISTING_LIBRARY_MGF",
                                                  true));
  listOptions.push_back(CommandLineParser::Option("ccms_input_spectradir",
                                                  "SPECTRA_DIR",
                                                  true));
  listOptions.push_back(CommandLineParser::Option("ccms_results_dir",
                                                  "RESULTS_DIR",
                                                  true));
  listOptions.push_back(CommandLineParser::Option("ccms_newresults_dir",
                                                  "NEWLIBRARYRESULTS_DIR",
                                                  true));
  listOptions.push_back(CommandLineParser::Option("ccms_input_spectradatafile",
                                                  "INPUTSPECTRA_DATAFILE",
                                                  true));
  listOptions.push_back(CommandLineParser::Option("ccms", "CCMS", true));

  CommandLineParser clp(argc, argv, 1, listOptions);
  string parser_error;
  if (!clp.validate(parser_error))
  {
    ERROR_MSG(parser_error);
    return -2;
  }

  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);

  if (commandLineParams.exists("RESUME_FLAG"))
  {
    resumeFlag = true;
  }
  DEBUG_VAR(resumeFlag);

  if (commandLineParams.exists("GRID_EXECUTION"))
  {
    gridExecutionFlag = true;
  }
  DEBUG_VAR(gridExecutionFlag);

  // Load the parameter file
  ParameterList ip;
  std::string paramsFileName = argv[1];
  std::string::size_type idx = paramsFileName.rfind(".");
  if (idx != std::string::npos && paramsFileName.substr(idx + 1) == "xml")
  {
    if (!ip.readFromProteosafeXMLFile(argv[1]))
    {
      ERROR_MSG("Can not read parameter file [" << argv[1] << "].");
      return -3;
    }
  }
  else
  {
    if (!ip.readFromFile(argv[1]))
    {
      ERROR_MSG("Can not read parameter file [" << argv[1] << "].");
      return -3;
    }
  }

  // Combine the command line parameters to the file ones
  //   Command line parameters take precedence (hence the overwrite flag set)
  ip.addList(commandLineParams, true);
  addDefaultParameterValues(ip);

  int logLevel = ip.getValueInt("LOG_LEVEL", 0);
  if (ip.exists("LOG_FILE_NAME"))
  {
    string logFileName = ip.getValue("LOG_FILE_NAME");
    Logger::setDefaultLogger(Logger::getLogger(logFileName, logLevel));
  }
  else
  {
    Logger::setDefaultLogger(Logger::getLogger(logLevel));
  }

  if (!ip.exists("EXE_DIR"))
  {
    // extract EXE_DIR from command line
    string exeDir(argv[0]);

    string mainFlrStr = "/main_flr";

    exeDir.erase(exeDir.length() - mainFlrStr.length(), mainFlrStr.length());

    mainFlrStr = "/ExecFramework/FalseLocalizationRates";

    if (exeDir.rfind(mainFlrStr) == exeDir.length() - mainFlrStr.length())
    {
      exeDir.erase(exeDir.length() - mainFlrStr.length(), mainFlrStr.length());
    }

    ip.setValue("EXE_DIR", exeDir);
  }

  if (ip.exists("EXE_DIR"))
  {
    string exeDir = ip.getValue("EXE_DIR");

    // if path begins with '~', exit program.
    if (exeDir[0] == '~')
    {
      cout
          << "EXE_DIR path begins with tilde (~). Paths beginning with tilde are not supported."
          << endl;
      exit(0);
    }

    // In case there is a "/" at the end of EXE_DIR.. remove it
    if (exeDir.length() > 2 && exeDir[exeDir.length() - 1] == '/')
    {
      exeDir = exeDir.substr(0, exeDir.length() - 1);
      ip.setValue("EXE_DIR", exeDir);
    }
  }

  if (ip.exists("GRID_EXE_DIR"))
  {
    // In case there is a "/" at the end of EXE_DIR.. remove it
    string gridExeDir = ip.getValue("GRID_EXE_DIR");
    if (gridExeDir.length() > 2 && gridExeDir[gridExeDir.length() - 1] == '/')
    {
      gridExeDir = gridExeDir.substr(0, gridExeDir.length() - 1);
      ip.setValue("GRID_EXE_DIR", gridExeDir);
    }
  }

  if (ip.exists("GRID_SGE_EXE_DIR"))
  {
    // In case there is a "/" at the end of GRID_SGE_EXE_DIR.. remove it
    string gridSgeExeDir = ip.getValue("GRID_SGE_EXE_DIR");
    if (gridSgeExeDir.length() > 2 && gridSgeExeDir[gridSgeExeDir.length() - 1]
        == '/')
    {
      gridSgeExeDir = gridSgeExeDir.substr(0, gridSgeExeDir.length() - 1);
      ip.setValue("GRID_SGE_EXE_DIR", gridSgeExeDir);
    }
  }

  if (!performFlr(ip, gridExecutionFlag, resumeFlag))
  {
    ERROR_MSG("Problem encountered during flr stage" );
    exit(-1);
  }

  DEBUG_VAR(ip.exists("GRID_NUMNODES"));

  if (!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0
      && (!resumeFlag && !gridExecutionFlag))
  {
    // If we are doing a grid execution (and are not actually on the grid)
    //    and we aren't resuming... then exit (we'll resume execution later)
    DEBUG_MSG("Files for grid execution have been saved.");
    DEBUG_MSG("Restart with -z option when grid execution has been completed.");
    exit(0);
  }

}

