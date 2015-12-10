// Header Include
#include "ExecMsCluster.h"

// Module Includes
#include "Logger.h"

#include "spectrum.h"
#include "utils.h"

#include "ClusterData.h"


// System Includes
#include <stdio.h>
#include <string.h>

// Defines used by the module

#define WORKING_DIRECTORY         "spectra"
#define FILES_PKLBIN              "./specs_ms_*.pklbin"


using namespace specnets;
using namespace spsReports;
using namespace std;

// -------------------------------------------------------------------------
ExecMsCluster::ExecMsCluster(void) : m_combinedSpecSet(NULL), m_spectrumNames(NULL)
{
    m_name            = "ExecMsCluster";
    m_type            = "ExecMsCluster";
    m_spectrumNames   = new vector<string>();
    m_combinedSpecSet = new SpecSet();
    m_dataOwned       = true;
}
// -------------------------------------------------------------------------
ExecMsCluster::ExecMsCluster(const ParameterList & params) : m_combinedSpecSet(NULL), m_spectrumNames(NULL)
{
    m_params          = params;
    m_name            = "ExecMsCluster";
    m_type            = "ExecMsCluster";
    m_spectrumNames   = new vector<string>();
    m_combinedSpecSet = new SpecSet();
    m_dataOwned       = true;
}

ExecMsCluster::ExecMsCluster(const ParameterList & inputParams, SpecSet *specSet, vector<string> *spectrumNames)
{
    m_spectrumNames   = spectrumNames ;
    m_combinedSpecSet = specSet;
    m_params          = inputParams;
    m_name            = "ExecMsCluster";
    m_type            = "ExecMsCluster";
    m_dataOwned       = false;
    m_combinedSpecSet->resize(0);
}

// -------------------------------------------------------------------------
ExecMsCluster::~ExecMsCluster(void)
{
  if(m_dataOwned) {
    delete m_spectrumNames;
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

  // get working directory
  string lastWorkingDirectory = getcwd(NULL, 1024);

  // Get CLUSTER_MIN_SIZE's value
  int clusterMinSize = getInt(m_params.getValue("CLUSTER_MIN_SIZE").c_str());

  // Test for cluster #. If > 0 we need to ...
  if (clusterMinSize > 0) {

    // Create a couple of folders
    string makeDirectory = "mkdir ";
    // create TMP directory
    string aux = makeDirectory;
    aux += TMP_DIRECTORY;
    DEBUG_MSG("system call: " << aux);
    //system(aux.c_str());
    mkdir_if_not_exist(aux.c_str());
    
    // create out directory
    aux = makeDirectory;
    aux += OUT_DIRECTORY;
    DEBUG_MSG("system call: " << aux);
    //system(aux.c_str());
    mkdir_if_not_exist(aux.c_str());

    // create out/mgf directory
    aux = makeDirectory;
    aux += MGF_DIRECTORY;
    DEBUG_MSG("system call: " << aux);
    //system(aux.c_str());
    mkdir_if_not_exist(aux.c_str());

    // create out/clust directory
    aux = makeDirectory;
    aux += CLUST_DIRECTORY;
    DEBUG_MSG("system call: " << aux);
    //system(aux.c_str());
    mkdir_if_not_exist(aux.c_str());


    // EXE_DIR
    std::string exeDir = m_params.getValue("EXE_DIR");
    // rtrim.
    rtrim(exeDir);

    // Call MsCluster
    int ret = callMsCluster(exeDir);
    if(ret != 0) {
      ERROR_MSG("Error executing MsCluster");
      return false;
    }

    // Delete previously generated mgf files DEBUG
    readFilesFromFile("cluster_files.txt", mgfFiles);
    aux = ".";
    removeFolderFiles(aux, mgfFiles);

    // Remove tmp directory files
    aux = TMP_DIRECTORY;
    removeFolder(aux);

    // out/mgf files info
    ret = chdir(MGF_DIRECTORY);
    // Read directory listing
    //outMgfDirectoryListing = directoryContents("./*.mgf", true);
    outMgfDirectoryListing = directoryContents(".", "", ".mgf", true);

    // Output mgf file list
    for(int i = 0 ; i < outMgfDirectoryListing.size() ; i++)
      DEBUG_MSG(outMgfDirectoryListing[i]);


    // Move back to base directory
    aux = "../..";
    chdir(aux.c_str());

  } // if (clusterMinSize > 0)

  // Merge input files
  chdir("./out/mgf");

  bool ret = mergeInputFiles(outMgfDirectoryListing );
  if(ret == false) {
    ERROR_MSG("Error merging input files.");
    return false;
  }

  string aux = "../..";
  chdir(aux.c_str());


  // Delete "clusters_*.bin" files

  vector<string> clusterDirectoryListing;
  // get directory listing for cluster_*.bin files
  //clusterDirectoryListing = directoryContents("clusters_*.bin", true);
  clusterDirectoryListing = directoryContents(".", "clusters_", ".bin", true);
  // delete aquired files
  for(int i = 0 ; i < clusterDirectoryListing.size() ; i++)
    remove(clusterDirectoryListing[i].c_str());


  // Return to previous woking directory
  chdir(lastWorkingDirectory.c_str());

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
  if (m_params.exists("OUTPUT_SPECS")) {
    // get filename
    string outSpecsFilename = m_params.getValue("OUTPUT_SPECS");
    // save the specs file
    m_combinedSpecSet->SaveSpecSet_pklbin(outSpecsFilename.c_str());
  }
  
  // generate the <sps>/clusterData.bin file
  ClusterData m_clusterData;
  string projectDir = ".";
  if(m_params.exists("PROJECT_DIR"))
    projectDir = m_params.getValue("PROJECT_DIR");
  m_clusterData.load(projectDir, true, true);
  
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
  
  replace(msClusterCommand.begin(), msClusterCommand.end(), DIR_SEP_OTHER, DIR_SEP);

  // Executable name
  msClusterCommand += DIR_SEP;
  msClusterCommand += "MsCluster_bin";  
  
  // --list
  msClusterCommand += " --list ./spectra/cluster_files.txt";

  // set tmp and out dirs
  msClusterCommand += " --tmp-dir ./spectra/tmp";
  msClusterCommand += " --out-dir ./spectra/out";

  // --output-name
  msClusterCommand += " --output-name clusters";
  // --window
  msClusterCommand += " --window ";
  msClusterCommand += parseFloat(pmTol, 7);

  // --fragment-tolerance

  // --model-dir
  msClusterCommand += " --model-dir ";
  msClusterCommand += exeDir;
  msClusterCommand += "/Models_mscluster";

  // --fragment-tolerance
  msClusterCommand += " --fragment-tolerance ";
  msClusterCommand += m_params.getValue("TOLERANCE_PEAK");

  // CLUSTER_MODEL
  if (m_params.exists("CLUSTER_MODEL") != 0) {
    msClusterCommand += " --model ";
    msClusterCommand += m_params.getValue("CLUSTER_MODEL");

    if (m_params.getValue("GUESS_CHARGE","no") == "yes") {
      msClusterCommand += " --assign-charges ";
    }

    if (m_params.getValue("CORRECT_PM","no") == "yes") {
      msClusterCommand += " --correct-pm ";
    }
  }

  // MIN_SPECTRUM_QUALITY
  if (m_params.exists("MIN_SPECTRUM_QUALITY") != 0) {
    msClusterCommand += " --sqs ";
    msClusterCommand += m_params.getValue("MIN_SPECTRUM_QUALITY");
  }

  // CLUSTER_PMTOL_PPM
  if (m_params.exists("CLUSTER_PMTOL_PPM") != 0) {
    msClusterCommand += " --precursor-ppm ";
    msClusterCommand += m_params.getValue("CLUSTER_PMTOL_PPM");
  }

  if (m_params.exists("CLUSTER_MIX_PROB") != 0) {
    msClusterCommand += " --mixture-prob ";
    msClusterCommand += m_params.getValue("CLUSTER_MIX_PROB");
  }

  msClusterCommand += " --memory-gb 0.35 ";

  // Force MsCluster_bin output into a file
  msClusterCommand += " > .";
  msClusterCommand += DIR_SEP;
  msClusterCommand += "log_mscluster.txt";

  //Send to sdtout DEBUG
  DEBUG_MSG("call msCluster: " << msClusterCommand);

  // Do the system call to invoke MsCluster_bin. ret contains the return value:
  // -1 for error, other for MsCluster_bin return value
  return system(msClusterCommand.c_str());
}
////////////////////////////////////////////////////////////////////////////////
void ExecMsCluster::removeFolderFiles(string &dir, vector<string> &files)
{
  for(int i = 0 ; i < files.size() ; i++) {
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
  for(int i = 0 ; i < outMgfDirectoryListing.size() ; i++) {
    // Define a SpecSet for pklbin loading
    SpecSet specSet;
    // Load specta data and test for read errors
    if(specSet.LoadSpecSet_mgf(outMgfDirectoryListing[i].c_str()) <= 0)
      return false;
      
    // Merge the file into allSpecSet
    m_combinedSpecSet->insert( m_combinedSpecSet->end(), specSet.begin(), specSet.end() );
    //for(int j = 0 ; j < specSet.specs.size() ; j++)
    //  m_combinedSpecSet->specs.push_back(specSet.specs[j]);

    // Delete the file
    //remove(outMgfDirectoryListing[i].c_str());
  }

  // Default scan numbers: [1,num clusters]
  for(unsigned int i=0; i<m_combinedSpecSet->size(); i++)
    {
      
      (*m_combinedSpecSet)[i].scan = i+1;
	}

  return true;
}
////////////////////////////////////////////////////////////////////////////////
