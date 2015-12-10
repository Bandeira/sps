// Header Include
#include "ExecMsCluster.h"

// Module Includes
#include "Logger.h"

#include "spectrum.h"
#include "utils.h"

#include "ClusterData.h"
#include "ReportDefines.h"

#include "Specific.h"

// System Includes
#include <stdio.h>
#include <string.h>

using namespace specnets;
using namespace spsReports;
using namespace std;

// -------------------------------------------------------------------------
ExecMsCluster::ExecMsCluster(void) :
  m_combinedSpecSet(NULL)
{
  m_name = "ExecMsCluster";
  m_type = "ExecMsCluster";
  m_combinedSpecSet = new SpecSet();
  m_dataOwned = true;
}
// -------------------------------------------------------------------------
ExecMsCluster::ExecMsCluster(const ParameterList & params) :
  m_combinedSpecSet(NULL)
{
  m_params = params;
  m_name = "ExecMsCluster";
  m_type = "ExecMsCluster";
  m_combinedSpecSet = new SpecSet();
  m_dataOwned = true;
}

ExecMsCluster::ExecMsCluster(const ParameterList & inputParams,
                             const string& inputFiles,
                             SpecSet *specSet)
{
  m_combinedSpecSet = specSet;
  m_params = inputParams;
  m_name = "ExecMsCluster";
  m_type = "ExecMsCluster";
  m_dataOwned = false;
  m_combinedSpecSet->resize(0);
  m_inputFileList = inputFiles;
}

// -------------------------------------------------------------------------
ExecMsCluster::~ExecMsCluster(void)
{
  if (m_dataOwned)
  {
    delete m_combinedSpecSet;
  }
}
// -------------------------------------------------------------------------
ExecBase * ExecMsCluster::clone(const ParameterList & input_params) const
{
  return new ExecMsCluster(input_params);
}
// -------------------------------------------------------------------------
bool ExecMsCluster::invoke(void)
{
  vector<string> outMgfDirectoryListing;

  DEBUG_MSG("Entering  ExecMsCluster::invoke()");

  // check for clusterData file. If present, exit
  //ClusterData m_clusterData;
  //string projectDir = ".";
  //if (m_params.exists("PROJECT_DIR"))
  //	projectDir = m_params.getValue("PROJECT_DIR");
  //if(m_clusterData.check(projectDir))
  //  return true;


  // get working directory
  string lastWorkingDirectory = getcwd(NULL, 1024);

  // Get CLUSTER_MIN_SIZE's value
  int clusterMinSize = getInt(m_params.getValue("CLUSTER_MIN_SIZE").c_str());

  // Test for cluster #. If > 0 we need to ...
  if (clusterMinSize > 0)
  {

    formatMsClusterInput();

    // Create a couple of folders
    string makeDirectory = "mkdir ";
    // create TMP directory
    string aux = makeDirectory;
    aux += MS_CLUSTER_TEMP_DIR;
    DEBUG_MSG("system call: " << aux);
    //system(aux.c_str());
    mkdir_if_not_exist(aux.c_str());

    // create out directory
    aux = makeDirectory;
    aux += MS_CLUSTER_OUT_DIR;
    DEBUG_MSG("system call: " << aux);
    //system(aux.c_str());
    mkdir_if_not_exist(aux.c_str());

    // create out/mgf directory
    aux = makeDirectory;
    aux += MS_CLUSTER_MGF_DIR;
    DEBUG_MSG("system call: " << aux);
    //system(aux.c_str());
    mkdir_if_not_exist(aux.c_str());

    // create out/clust directory
    aux = makeDirectory;
    aux += MS_CLUSTER_DATA_SUBPATH;
    DEBUG_MSG("system call: " << aux);
    //system(aux.c_str());
    mkdir_if_not_exist(aux.c_str());

    // EXE_DIR
    std::string exeDir = m_params.getValue("EXE_DIR");
    // rtrim.
    rtrim(exeDir);

    // Call MsCluster
    int ret = callMsCluster(exeDir);
    if (ret != 0)
    {
      ERROR_MSG("Error executing MsCluster");
      return false;
    }

    // Remove tmp directory files
    aux = MS_CLUSTER_TEMP_DIR;
    removeFolder(aux);

    // out/mgf files info
    ret = chdir(MS_CLUSTER_MGF_DIR);
    // Read directory listing
    //outMgfDirectoryListing = directoryContents("./*.mgf", true);
    outMgfDirectoryListing = directoryContents(".", "", ".mgf", true);

    // Output mgf file list
    for (int i = 0; i < outMgfDirectoryListing.size(); i++)
      DEBUG_MSG(outMgfDirectoryListing[i]);

    // Move back to base directory
    aux = "../..";
    if (chdir(aux.c_str()) != 0)
    {
      ERROR_MSG("Failed to change directory");
      return false;
    }

  } // if (clusterMinSize > 0)

  // Merge input files
  if (chdir("./out/mgf") != 0)
  {
    ERROR_MSG("Failed to change directory");
    return false;
  }

  bool ret = mergeInputFiles(outMgfDirectoryListing);
  if (ret == false)
  {
    ERROR_MSG("Error merging input files.");
    return false;
  }

  string aux = "../..";
  if (chdir(aux.c_str()) != 0)
  {
    ERROR_MSG("Failed to change directory");
    return false;
  }

  if (!formatMsClusterOutput())
  {
    return false;
  }

  // Delete "clusters_*.bin" files

  vector<string> clusterDirectoryListing;
  // get directory listing for cluster_*.bin files
  //clusterDirectoryListing = directoryContents("clusters_*.bin", true);
  clusterDirectoryListing = directoryContents(".", "clusters_", ".bin", true);
  // delete aquired files
  for (int i = 0; i < clusterDirectoryListing.size(); i++)
  {
    if (remove(clusterDirectoryListing[i].c_str()) != 0)
    {
      ERROR_MSG("Failed to remove file");
      return false;
    }
  }

  // Return to previous woking directory
  if (chdir(lastWorkingDirectory.c_str()) != 0)
  {
    ERROR_MSG("Failed to change directory");
    return false;
  }

  if (!cleanupMsClusterInput())
  {
    return false;
  }

  DEBUG_MSG("Exiting  ExecMsCluster::invoke()");

  return true;
}

// -------------------------------------------------------------------------
bool ExecMsCluster::loadInputData(void)
{
}
// -------------------------------------------------------------------------
bool ExecMsCluster::saveOutputData(void)
{
  // Save concatenated pklbin file
  if (m_params.exists("OUTPUT_SPECTRA"))
  {

    if (!ExecMergeConvert::saveSpecsetMultiple(m_params.getValue("OUTPUT_SPECTRA"),
                                               m_combinedSpecSet))
    {
      return false;
    }
  }

  // generate the <sps>/clusterData.bin file
  ClusterData m_clusterData;
  string projectDir = ".";
  if (m_params.exists("PROJECT_DIR"))
    projectDir = m_params.getValue("PROJECT_DIR");
  int ret = m_clusterData.load(projectDir, true, true);

  return (ret == OK);
}
// -------------------------------------------------------------------------
bool ExecMsCluster::saveInputData(std::vector<std::string> & filenames)
{
}
// -------------------------------------------------------------------------
bool ExecMsCluster::loadOutputData(void)
{
}
// -------------------------------------------------------------------------
vector<ExecBase*> const & ExecMsCluster::split(int numSplit)
{
  m_subModules.resize(0);
  return m_subModules;
}
// -------------------------------------------------------------------------
bool ExecMsCluster::merge(void)
{
}
// -------------------------------------------------------------------------
bool ExecMsCluster::validateParams(std::string & error)
{
  m_isValid = false;

  VALIDATE_PARAM_EXIST("EXE_DIR");
  VALIDATE_PARAM_EXIST("CLUSTER_MIN_SIZE");
  VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
  VALIDATE_PARAM_EXIST("INPUT_SPECS_DIR");

  m_isValid = true;
  return true;
}
bool ExecMsCluster::formatMsClusterInput()
{

  if (m_inputFileList.length() == 0)
  {
    ERROR_MSG("Input cluster file list not specified!!");
    return false;
  }
  string lines("");
  string line;
  ifstream infile;
  m_msClustInputFileList = m_inputFileList;
  m_msClustInputFileList += "-format_prec.txt";
  infile.open(m_inputFileList.c_str(), ios::binary);
  if (infile.is_open())
  {
    while (infile.good())
    {
      getline(infile, line);
      if (line.length() == 0)
      {
        break;
      }
      lines += line;
      lines += "\n";
    }
  }
  else
  {
    ERROR_MSG("Could not read file " << m_inputFileList);
    return false;
  }
  infile.close();
  vector<string> inputFiles;
  splitText(lines.c_str(), inputFiles, "\n");
  string outputContents("");
  m_inputSpectra.resize(inputFiles.size());
  int idxUse = 0;
  for (int i = 0; i < inputFiles.size(); i++)
  {
    string nextFile(inputFiles[i]);

    if (nextFile.length() == 0)
    {
      continue;
    }
    SpecSet specSet;
    if (!ExecMergeConvert::loadSpecset(nextFile, &specSet))
    {
      return false;
    }
    m_inputSpectra[idxUse].resize(specSet.size());

    for (unsigned int j = 0; j < specSet.size(); j++)
    {
      m_inputSpectra[idxUse][j].resize(0);
      m_inputSpectra[idxUse][j].copyNP(specSet[j]);

      if (specSet[j].parentMass > 5000)
      {
        WARN_MSG("MSCluster hack for spectrum "<<j<<" in "<<inputFiles[i]);
        float z = specSet[j].parentCharge;
        if (z > 0)
        {
          WARN_MSG("converting to m/z with fake charge 0");
          specSet[j].parentMass = (specSet[j].parentMass + AAJumps::massHion
              * (z - 1)) / z;
          specSet[j].parentCharge = 1;
        }
        else
        {
          WARN_MSG("deleting spectrum");
          specSet[j].parentMass = 0;
          specSet[j].resize(0);
        }
      }
    }
    nextFile += "-format_prec.mgf";
    if (!ExecMergeConvert::saveSpecset(nextFile, &specSet))
    {
      return false;
    }
    outputContents += nextFile;
    outputContents += "\n";
    idxUse++;
  }

  ofstream myfileOut(m_msClustInputFileList.c_str(), ios::binary);
  if (myfileOut.is_open())
  {
    myfileOut << outputContents;
    myfileOut.close();
  }
  else
  {
    ERROR_MSG("Could not write to file " << m_msClustInputFileList);
    return false;
  }

  return true;

}

bool ExecMsCluster::cleanupMsClusterInput()
{
  if (m_msClustInputFileList.length() == 0)
  {
    ERROR_MSG("Input MS-Cluster file list not specified!!");
    return false;
  }
  string lines("");
  string line;
  ifstream infile;
  infile.open(m_msClustInputFileList.c_str(), ios::binary);
  if (infile.is_open())
  {
    while (infile.good())
    {
      getline(infile, line);
      if (line.length() == 0)
      {
        break;
      }
      if (remove(line.c_str()) != 0)
      {
        ERROR_MSG("Failed to remove \'" << line << "\'");
// This should not be a fatal error
//        return false;
      }
      else
      {
        DEBUG_MSG("Successfully removed \'" << line << "\'");
      }
    }
  }
  else
  {
    ERROR_MSG("Could not read file " << m_msClustInputFileList);
    return false;
  }
  infile.close();
  if (remove(m_msClustInputFileList.c_str()) != 0)
  {
    ERROR_MSG("Failed to remove \'" << m_msClustInputFileList << "\'");
// This should not be a fatal error
//    return false;
  }
  else
  {
    DEBUG_MSG("Successfully removed \'" << m_msClustInputFileList << "\'");
  }
  return true;
}

bool ExecMsCluster::formatMsClusterOutput()
{
  ClusterData clusterData;
  string projectDir = ".";
  if (m_params.exists("PROJECT_DIR"))
    projectDir = m_params.getValue("PROJECT_DIR");
  clusterData.load(projectDir, false, true);
  DEBUG_TRACE;
  //DEBUG_VAR(m_combinedSpecSet->size());
  //DEBUG_VAR(m_inputSpectra.size());
  for (map<unsigned, list<pair<unsigned, unsigned> > >::iterator clustIt =
      clusterData.data.begin(); clustIt != clusterData.data.end(); clustIt++)
  {
    unsigned int newIdx = clustIt->first;
    unsigned int oldIdx = clustIt->second.front().second;
    unsigned int oldFileIdx = clustIt->second.front().first;
    unsigned int oldScan = (*m_combinedSpecSet)[newIdx].scan;
    (*m_combinedSpecSet)[newIdx].copyNP(m_inputSpectra[oldFileIdx][oldIdx]);
    (*m_combinedSpecSet)[newIdx].scan = oldScan;
  }

  return true;
}
////////////////////////////////////////////////////////////////////////////////
// -------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////
int ExecMsCluster::callMsCluster(string & exeDir)
{
  // TOLERANCE_PM
  float pmTol = 0.0;
  if (m_params.exists("TOLERANCE_PM") != 0)
    pmTol = getFloat(m_params.getValue("TOLERANCE_PM").c_str());

  // MsCluster command
  std::string msClusterCommand;

  // Path
  msClusterCommand = exeDir;

  replace(msClusterCommand.begin(),
          msClusterCommand.end(),
          DIR_SEP_OTHER,
          DIR_SEP);

  // Executable name
  msClusterCommand += DIR_SEP;
  msClusterCommand += "MsCluster_bin";

  // --list
  msClusterCommand += " --list ";
  msClusterCommand += m_msClustInputFileList;

  string outDir = m_params.getValue("MSCLUSTER_OUT_DIR");
  string tmpDir = m_params.getValue("MSCLUSTER_TMP_DIR");

  // set tmp and out dirs
  msClusterCommand += " --tmp-dir ";//./spectra/tmp";
  msClusterCommand += MS_CLUSTER_TEMP_DIR;
  msClusterCommand += " --out-dir ";//./spectra/out";
  msClusterCommand += MS_CLUSTER_OUT_DIR;

  // --output-name
  msClusterCommand += " --output-name clusters";
  // --window
  msClusterCommand += " --window ";
  msClusterCommand += parseFloat(pmTol, 7);

  // --fragment-tolerance

  // --model-dir
  msClusterCommand += " --model-dir ";
  msClusterCommand += exeDir;
  msClusterCommand += "/resources/Models_mscluster";

  // --fragment-tolerance
  msClusterCommand += " --fragment-tolerance ";
  msClusterCommand += m_params.getValue("TOLERANCE_PEAK");

  // CLUSTER_MODEL
  if (m_params.exists("CLUSTER_MODEL") != 0)
  {
    msClusterCommand += " --model ";
    msClusterCommand += m_params.getValue("CLUSTER_MODEL");

    if (m_params.getValue("GUESS_CHARGE", "no") == "yes")
    {
      msClusterCommand += " --assign-charges ";
    }

    if (m_params.getValue("CORRECT_PM", "no") == "yes")
    {
      msClusterCommand += " --correct-pm ";
    }
  }

  // MIN_SPECTRUM_QUALITY
  if (m_params.exists("MIN_SPECTRUM_QUALITY") != 0)
  {
    msClusterCommand += " --sqs ";
    msClusterCommand += m_params.getValue("MIN_SPECTRUM_QUALITY");
  }

  // CLUSTER_PMTOL_PPM
  if (m_params.exists("CLUSTER_PMTOL_PPM") != 0)
  {
    msClusterCommand += " --precursor-ppm ";
    msClusterCommand += m_params.getValue("CLUSTER_PMTOL_PPM");
  }

  if (m_params.exists("MSCLUSTER_MIX_PROB") != 0)
  {
    msClusterCommand += " --mixture-prob ";
    msClusterCommand += m_params.getValue("MSCLUSTER_MIX_PROB");
  }

  msClusterCommand += " --memory-gb 0.35 ";

  // Force MsCluster_bin output into a file
  msClusterCommand += " > .";
  msClusterCommand += DIR_SEP;
  msClusterCommand += "log_mscluster.txt";

  //Send to sdtout DEBUG
  DEBUG_MSG("call msCluster: " << msClusterCommand);

  //string libPath = exeDir + "/libs";
  //addEnvironmentVariable(msClusterCommand, "LD_LIBRARY_PATH", libPath);

  return spsSystem(msClusterCommand.c_str());
}
////////////////////////////////////////////////////////////////////////////////
void ExecMsCluster::removeFolderFiles(string &dir, vector<string> &files)
{
  for (int i = 0; i < files.size(); i++)
  {
    string aux = dir;
    aux += "/";
    aux += files[i];
    remove(aux.c_str());
  }
}
////////////////////////////////////////////////////////////////////////////////
void ExecMsCluster::removeFolder(string &dir)
{
  // Variable to store directory listing
  vector<string> tmpDirectoryListing;
  // Auxiliary variables for directory composition
  string aux;
  aux = dir;
  aux += "/*";
  // get directory contents
  tmpDirectoryListing = directoryContents(aux.c_str(), "", "", false);
  // delete files in the directory
  removeFolderFiles(dir, tmpDirectoryListing);
  // now remove the directory
  remove(dir.c_str());
}
////////////////////////////////////////////////////////////////////////////////
// Merge-load the input files
bool ExecMsCluster::mergeInputFiles(vector<string> &outMgfDirectoryListing)
{
  m_combinedSpecSet->resize(0);
  // Merge all pklbin files into one
  for (int i = 0; i < outMgfDirectoryListing.size(); i++)
  {
    // Define a SpecSet for pklbin loading
    SpecSet specSet;
    // Load specta data and test for read errors
    if (specSet.LoadSpecSet_mgf(outMgfDirectoryListing[i].c_str()) <= 0)
      return false;

    // Merge the file into allSpecSet
    m_combinedSpecSet->insert(m_combinedSpecSet->end(),
                              specSet.begin(),
                              specSet.end());
    //for(int j = 0 ; j < specSet.specs.size() ; j++)
    //  m_combinedSpecSet->specs.push_back(specSet.specs[j]);

    // Delete the file
    //remove(outMgfDirectoryListing[i].c_str());
  }

  // Default scan numbers: [1,num clusters]
  for (unsigned int i = 0; i < m_combinedSpecSet->size(); i++)
  {

    (*m_combinedSpecSet)[i].scan = i + 1;
  }

  return true;
}
////////////////////////////////////////////////////////////////////////////////
