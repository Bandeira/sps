//

// Module Includes
#include "CommandLineParser.h"
#include "Logger.h"
#include "ExecMsCluster.h"
#include "ExecGenoMS.h"
#include "ExecAlignment.h"
#include "ExecFilterAligns.h"
#include "ExecFilterPairs.h"
#include "ExecFilterStarPairs.h"
#include "ExecReportSpsplot.h"
#include "ExecReportProteinCoverage.h"
#include "ExecAssembly.h"
#include "ExecPrmScoring.h"
#include "ExecSpecProtAlign.h"
#include "ExecProtProtAlign.h"
#include "ExecTagSearch.h"
#include "ExecHomologyAssembly.h"
#include "ExecMainSpecnets.h"
#include "FileUtils.h"
#include "ParallelThreadedExecution.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"
#include "ParameterList.h"
#include "PeptideSpectrumMatchSet.h"

// Specnets Includes
#include "utils.h"  // stringSplit
#include "clusters.h"
#include "ClusterData.h"
#include "abruijn.h"
#include "copyright.h"


// System Includes
#include <fstream>
#include <iostream>
#include <string.h>
#include <string>
#include <unistd.h>
#include <algorithm>

using namespace specnets;
using namespace std;

enum Stage
{
  STAGE_BEGIN             = 0,
  STAGE_SCORING           = 1,
  STAGE_FILTERPAIRS       = 2,
  STAGE_FILTERALIGNS      = 3,
  STAGE_ALIGNMENT         = 4,
  STAGE_FILTERSTARPAIRS   = 5,
  STAGE_ASSEMBLY          = 6,
  STAGE_TAGSEARCH         = 7,
  STAGE_SPECNETS          = 8,
  STAGE_SPECPROTALIGN     = 9,
  STAGE_PROTPROTALIGN     = 10,
  STAGE_HOMOLOGYASSEMBLY  = 11,
  STAGE_GENOMS            = 12,
  STAGE_MERGE             = 13,
  STAGE_REPORT            = 14

};


// System Defines
#define TEST_VALID {          \
  DEBUG_VAR(isValid);         \
  if (!isValid) {             \
    ERROR_MSG(errorString);   \
    return false;             \
  }                           \
}

#define TEST_RETURN(__module, __var) {                                                    \
  if(__var.size() == 0) {                                                               \
    ERROR_MSG("invoking " << __module << ": empty return set \"" << #__var << "\".");     \
  /*  return false; */                                                                     \
  }                                                                                   \
}

#define TEST_RETURN_STATUS(__module) {                            \
  DEBUG_VAR( returnStatus);                                     \
  if (!returnStatus) {                                          \
    ERROR_MSG("invoking " << __module << " exited in error.");    \
    return false;                                              \
  }                                                             \
}

#define TEST_SAVE_OUPUT_DATA(__module) {                                      \
  DEBUG_VAR(returnStatus);                                                  \
  if (!returnStatus) {                                                      \
    ERROR_MSG("Saving output data for " << __module << " exited in error.");  \
    return false;                                                           \
  }                                                                         \
}

//-----------------------------------------------------------------------------
static string getCurrentDirectory(void)
{
  int dummy;
  dummy = system("pwd > pwd.tmp");
  ifstream ifs("pwd.tmp");
  char buf[1024];
  ifs.getline(buf, 1024);
  return buf;
}


//-----------------------------------------------------------------------------
string getProjPath(const ParameterList & pl, const string & addPath)
{
  string projDir = pl.getValue("PROJECT_DIR", "");
  bool relDir = pl.getValueBool("RELATIVE_DIR", true);
  return getPath(projDir, addPath, relDir);
}


//-----------------------------------------------------------------------------
string getCurrentTimeString(void)
{
  time_t rawtime;
  struct tm * timeinfo;
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  return asctime(timeinfo);
}

//-----------------------------------------------------------------------------
void addDefaultParameterValues(ParameterList &p)
{
	// Basic parameters
	p.addIfDoesntExist("TOLERANCE_PEAK",        "0.4");
	p.addIfDoesntExist("TOLERANCE_PM",          "1.5");
	p.addIfDoesntExist("RESOLUTION",            "0.1");

	// Preprocessing parameters
	p.addIfDoesntExist("CLUSTER_MIN_SIZE",      "1");
	p.addIfDoesntExist("CLUSTER_MODEL",         "LTQ_TRYP");
	p.addIfDoesntExist("PEPNOVO_MODEL",         "CID_IT_TRYP");
	p.addIfDoesntExist("MIN_SPECTRUM_QUALITY",  "0.15");
	p.addIfDoesntExist("CORRECT_PM",            "no");
	p.addIfDoesntExist("GUESS_CHARGE",          "no");

	// Alignment parameters
	p.addIfDoesntExist("AA_DIFF_COUNT",         "2");
	p.addIfDoesntExist("MIN_SHIFT",             "0");
	p.addIfDoesntExist("MIN_MOD_MASS",          "-100");
	p.addIfDoesntExist("MAX_MOD_MASS",          "100");
	p.addIfDoesntExist("MAX_NUM_MODS",          "1");
	p.addIfDoesntExist("MIN_RATIO",             "0.35");

	p.addIfDoesntExist("MAX_PVALUE",            "0.05");
	p.addIfDoesntExist("MIN_MATCHED_PEAKS",     "4");
	p.addIfDoesntExist("MAX_AA_JUMP",           "2");
	p.addIfDoesntExist("MIN_OVERLAP_AREA",      "0.45");
	p.addIfDoesntExist("PENALTY_PTM",           "-200");  // Set to "0" for SpecNets
	p.addIfDoesntExist("PENALTY_SAME_VERTEX",   "-1000000");
	p.addIfDoesntExist("FILTER_TRIGS",          "yes");  // Set to "no" for SpecNets
	p.addIfDoesntExist("PARTIAL_OVERLAPS",      "1");  // Set to "0" for SpecNets
	p.addIfDoesntExist("TAGS_FILTER",           "");
	p.addIfDoesntExist("TAGS_MATCH_FLANK",      "1");
	p.addIfDoesntExist("TAGS_MATCH_COUNT",      "2");

	// Comparative Shotgun Protein Sequencing (CSPS) parameters
	//p.addIfDoesntExist("CLUSTALW_EXE_DIR",      "");
	p.addIfDoesntExist("CLUSTALW_MINSCORE",     "250");
//	p.addIfDoesntExist("FORCE_REFERENCE",       "-1");

	// De novo sequencing parameters
	p.addIfDoesntExist("SPSPATH_MIN_NUM_PEAKS",       "5");
	p.addIfDoesntExist("SPSPATH_MIN_NUM_SPECS",       "2");
	p.addIfDoesntExist("SPS_MIN_EDGES_TO_COMPONENT",  "1");

	// tagsearch/matchma parameters
	p.addIfDoesntExist("TAG_LEN",                     "6");
	p.addIfDoesntExist("DOUBLE_AA_JUMPS",             "1");
	p.addIfDoesntExist("MATCH_TAG_FLANKING_MASSES",   "0");  // Set to 2 for SpecNets
	p.addIfDoesntExist("MAX_NUM_TAGS",                "0");
	p.addIfDoesntExist("MAX_NUM_MODS",                "2");
	p.addIfDoesntExist("MIN_MATCHED_PEAKS_DB",        "6");
	p.addIfDoesntExist("TAG_MATCH_TOP_SCORING_ONLY",  "1");

	// Networks parameters (pathproj)
//	p.addIfDoesntExist("MIN_PERC_EXPINT",   "0.01");
//	p.addIfDoesntExist("MIN_PERC_TP",       "0.01");

	// Grid parameters
	p.addIfDoesntExist("GRID_TYPE",         "sge");
	p.addIfDoesntExist("GRID_NUMNODES",     "0");
	p.addIfDoesntExist("GRID_NUMCPUS",      "1");
	p.addIfDoesntExist("GRID_EXE_DIR",      "");
	p.addIfDoesntExist("GRID_SGE_EXE_DIR",  "");
	p.addIfDoesntExist("GRID_PARAMS",       "-l h_vmem=1G");

  string currDir = getCurrentDirectory();
  p.addIfDoesntExist("PROJECT_DIR", currDir);
  p.addIfDoesntExist("RELATIVE_DIR",      "1");

	// Reporting parameters
	p.addIfDoesntExist("REPORT_DIR",        "report");
  p.addIfDoesntExist("SPS_PROJECTS",      "./sps_projects.txt");

}

//-----------------------------------------------------------------------------
bool loadInitialData(ParameterList & ip,
                     string convertCmd,
	                   string inputFileNames,
	                   string outputPrefix,
	                   SpecSet *ms2spectra = (SpecSet *)0)
{
  DEBUG_TRACE;
  //---------------------------------
  // Load spectra
  //---------------------------------

  // Get current working directory -- needed for report interactivity
  string currentWorkindDir;
  currentWorkindDir = getcwd(NULL, 1024);
  string aux; // used for path composition

  // List of generated data files
  vector<string> mfgFileList, pklbinFileList, binFileList;
  //
	vector<string> filenames, filenameTokens;
	string extension;
	ostringstream cmd, cmd2;
	stringSplit(inputFileNames,filenames, ";");

  string aux2 = getProjPath(ip, "./spectra/input_index.txt");
	ofstream input_index(aux2.c_str() ,ios_base::out);
	if(!input_index.is_open()) {
    ERROR_MSG("ERROR opening spectra/input_index.txt");
    return false;
	}

  DEBUG_TRACE;
	unsigned int baseIdx, totalSpectrumCount = 0;
	SpecSet *tmpSpecs = new SpecSet;
	for(unsigned int i=0; i<filenames.size(); i++) {
    DEBUG_VAR(filenames[i]);
		input_index << filenames[i] <<"\n";
		stringSplit(filenames[i],filenameTokens,".");
		// check for for invalid filename/extension
		if(filenameTokens.size() == 0)
		  continue;

		extension = filenameTokens[filenameTokens.size()-1];
		std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

		cmd.str("");
		cmd2.str("");
		if(extension == "pklbin") {
      DEBUG_VAR(extension);
		  // Load file (pklbin)
			if(tmpSpecs->LoadSpecSet_pklbin(filenames[i].c_str()) <= 0) {
        ERROR_MSG("Could not open file " << filenames[i]);
        delete tmpSpecs;
        return false;
      }
			totalSpectrumCount += tmpSpecs->size();
			// generate sequential name, 1 based index
			cmd << outputPrefix << i+1 << ".pklbin";
			// generate sequential name, 1 based index, for .bin file output
			cmd2 <<  outputPrefix << i+1 << ".bin";
			// save it with the name we want (pklbin)
			tmpSpecs->SaveSpecSet_pklbin(cmd.str().c_str(), cmd2.str().c_str());

      //=================================================
      // LARS - TEMPORARY MEASURE SO REPORTS KEEP WORKING
      //=================================================
			tmpSpecs->SaveSpecSet_bin(cmd2.str().c_str());
      //=================================================
      // END OF TEMPORARY MEASURE SO REPORTS KEEP WORKING
      //=================================================

			// add saved name to pklbin file list
			//aux = currentWorkindDir;
			//aux += '/';
			//aux += cmd.str();
			aux = cmd.str();
      DEBUG_VAR(aux);
			pklbinFileList.push_back(aux);
  		// add saved .bin file name to bin file list
  		//aux = currentWorkindDir;
    	//aux += '/';
    	//aux +=  cmd2.str();
			binFileList.push_back(aux);
		} else {
      DEBUG_VAR(extension);
		  // If it is an mgf file
			if(extension == "mgf") {
		  // Load file (mgf)
				if(tmpSpecs->LoadSpecSet_mgf(filenames[i].c_str()) <= 0) {
          ERROR_MSG("Could not open file " << filenames[i]);
          delete tmpSpecs;
          return false;
				}
				totalSpectrumCount += tmpSpecs->size();
  			// generate sequential name, 1 based index
				cmd << outputPrefix << i+1 << ".pklbin";
  			// generate sequential name, 1 based index, for .bin file output
  			cmd2 <<  outputPrefix << i+1 << ".bin";
	  		// save it with as pklbin, with the name we want
				tmpSpecs->SaveSpecSet_pklbin(cmd.str().c_str(), cmd2.str().c_str());

        //=================================================
        // LARS - TEMPORARY MEASURE SO REPORTS KEEP WORKING
        //=================================================
  			tmpSpecs->SaveSpecSet_bin(cmd2.str().c_str());
        //=================================================
        // END OF TEMPORARY MEASURE SO REPORTS KEEP WORKING
        //=================================================

				// add saved file name to pklbin file list
  			pklbinFileList.push_back(cmd.str());
				// add saved .bin file name to bin file list
    		//aux = currentWorkindDir;
      	//aux += '/';
      	//aux +=  cmd2.str();
      	aux =  cmd2.str();
  			binFileList.push_back(aux);
			} else {
				cmd.str("");
				if(extension=="pkl" or extension=="mzxml") {
					cmd << convertCmd << extension << " " << filenames[i];
				} else {
					ERROR_MSG("Could not parse " << filenames[i] <<" due to unknown extension file type " << extension);
					ms2spectra->resize(0);
					delete tmpSpecs;
					return false;
				}

				cmd << " " << outputPrefix << i+1;
				int dummy;
				dummy = system(cmd.str().c_str());
				if(dummy != 0) {
          ERROR_MSG("convert failed for file " << filenames[i]);
          DEBUG_MSG("convert command: " << cmd.str());
          delete tmpSpecs;
          return false;
				}
				cmd.str("");
				cmd << outputPrefix << i+1 << ".pklbin";
				totalSpectrumCount += tmpSpecs->LoadSpecSet_pklbin(cmd.str().c_str(), true);
  			pklbinFileList.push_back(cmd.str());
  			if(extension=="pkl" or extension=="mzxml") {
  				cmd.str("");
	  			cmd << outputPrefix << i+1 << ".bin";
    			//aux = currentWorkindDir;
    			//aux += '/';
    			//aux += cmd.str();
  	  		binFileList.push_back(aux);
  	  	}
			}
		}
		DEBUG_VAR(totalSpectrumCount);
	}
	if(ms2spectra) {
		ms2spectra->resize(totalSpectrumCount);
		DEBUG_VAR(ms2spectra->size());
		totalSpectrumCount = 0;
		for(unsigned int i=0; i<pklbinFileList.size(); i++) {
			DEBUG_VAR(pklbinFileList[i]);
			if(tmpSpecs->LoadSpecSet_pklbin(pklbinFileList[i].c_str()) <= 0) {
				ERROR_MSG("Could not open file " << pklbinFileList[i]);
				delete tmpSpecs;
				return false;
			}
			DEBUG_VAR(tmpSpecs->size());
			/*  DEBUG_VAR((*tmpSpecs)[0].annotation.size());
			DEBUG_VAR((*tmpSpecs)[0].annotation_peptide.size());
			DEBUG_VAR((*tmpSpecs)[0].ionTypes.size());
			 */
			for(unsigned int j=0; j<tmpSpecs->size(); j++) {
				(*ms2spectra)[totalSpectrumCount+j] = (*tmpSpecs)[j];
				(*ms2spectra)[totalSpectrumCount+j].scan = totalSpectrumCount+j+1;
			}
			totalSpectrumCount += tmpSpecs->size();
			DEBUG_VAR(totalSpectrumCount);

			//

		}
	}
	input_index.close();
	delete tmpSpecs;

  DEBUG_TRACE;

  // save pklbin file list
  aux = getProjPath(ip, "./pklbin_files.txt");
  if(!writeFileIndex(aux.c_str(), pklbinFileList)) {
    ERROR_MSG("Could not open \"pklbin_files.txt\".");
    return false;
  }

  // save bin file list
  if(binFileList.size()) {
    string aux = getProjPath(ip, "./bin_files.txt");
    if(!writeFileIndex(aux.c_str(), binFileList)) {
      ERROR_MSG("Could not open \"bin_files.txt\".");
      return false;
    }
  }

	// Convert files to MGF and save MGFs index file

  // Load files in directory and save them as mgf
  for(int i = 0 ; i < pklbinFileList.size() ; i++) {
    // Specset object - used to load specset from file
    SpecSet specSet;
    // Load specta data and test for read errors
    if(specSet.LoadSpecSet_pklbin(pklbinFileList[i].c_str()) <= 0) {
      ERROR_MSG("Could not open file " << pklbinFileList[i]);
      return 0;
    }

    // Hack to get around MSCluster crash while it's being fixed
    for(unsigned int j=0; j<specSet.size(); j++) {
    	if(specSet[j].parentMass > 5000) {
    		cerr<<"WARNING: MSCluster hack for spectrum "<<j<<" in "<<pklbinFileList[i].c_str()<<": ";
    		float z=specSet[j].parentCharge;
    		if(z > 0) {
    			cerr << "converting to m/z with fake charge 0\n"; cerr.flush();
    			specSet[j].parentMass = (specSet[j].parentMass + AAJumps::massHion*(z-1) ) / z;
    			specSet[j].parentCharge = 1;
    		} else {
    			cerr << "deleting spectrum\n"; cerr.flush();
    			specSet[j].parentMass = 0;
    			specSet[j].resize(0);
    		}
    	}
    }

    // Get dot before file name extension position
    int lastDotPosition = pklbinFileList[i].string::find_last_of(".");
    // Get file name without extension
    string mgfFileName = pklbinFileList[i].substr(0, lastDotPosition + 1);
    // add extension
    mgfFileName += "mgf";
    // Save the file in mgf format
    specSet.SaveSpecSet_mgf(mgfFileName.c_str());
    // add mgf file name to list for later use (delete the file)
    mfgFileList.push_back(mgfFileName);
    // output operation DEBUG
    DEBUG_MSG("Convert pklbin to mgf: " << pklbinFileList[i] << " --> " << mgfFileName);
  }

  // save mgf file list
  aux = getProjPath(ip, "./spectra/cluster_files.txt");
  if(!writeFileIndex(aux.c_str(), mfgFileList)) {
    ERROR_MSG("Could not open \"cluster_files.txt\".");
    return false;
  }


  DEBUG_TRACE;
  return true;
}

//-----------------------------------------------------------------------------
bool performMsCluster(ParameterList & ip,
		                  SpecSet &clustSpectra,
		                  vector<string> &spectrumNames)
{
  DEBUG_TRACE;

  ParameterList msclusterParams;
  msclusterParams.addIfExists(ip,"CLUSTER_MIN_SIZE");
  msclusterParams.addIfExists(ip,"EXE_DIR");
  msclusterParams.addIfExists(ip,"TOLERANCE_PEAK");
  msclusterParams.addIfExists(ip,"TOLERANCE_PM");
  msclusterParams.addIfExists(ip,"CLUSTER_MODEL");
  msclusterParams.addIfExists(ip,"MIN_SPECTRUM_QUALITY");
  msclusterParams.addIfExists(ip,"INPUT_SPECS_MS");

  msclusterParams.addIfExists(ip, "PROJECT_DIR");

  msclusterParams.setValue("INPUT_SPECS_DIR",     "spectra");
  msclusterParams.setValue("OUTPUT_SPECS",        "spectra/specs_ms.pklbin");

  DEBUG_TRACE;
  stringstream aux;
  msclusterParams.print(aux);
  DEBUG_MSG(aux.str());
  DEBUG_TRACE;
  ExecMsCluster moduleMsCluster(msclusterParams, &clustSpectra, &spectrumNames);

  string errorString;
  bool isValid = moduleMsCluster.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  bool returnStatus = moduleMsCluster.invoke();
  // Test for return status
  TEST_RETURN_STATUS("moduleMsCluster");

  DEBUG_TRACE;

  moduleMsCluster.saveOutputData();

  return true;
}

//-----------------------------------------------------------------------------
bool performScoring(ParameterList & ip,
                    SpecSet & inputSpectra,
                    SpecSet & outputSpectra)
{
  DEBUG_TRACE;

  ParameterList scoringParams;
  scoringParams.addIfExists(ip,"EXE_DIR");
  scoringParams.addIfExists(ip,"TOLERANCE_PEAK");
  scoringParams.addIfExists(ip,"TOLERANCE_PM");
  scoringParams.addIfExists(ip,"CLUSTER_MODEL");
  scoringParams.addIfExists(ip,"MIN_SPECTRUM_QUALITY");
  scoringParams.addIfExists(ip,"GUESS_CHARGE");
  scoringParams.addIfExists(ip,"CORRECT_PM");
  scoringParams.addIfExists(ip,"PEPNOVO_MODEL");

  scoringParams.addIfExists(ip, "PROJECT_DIR");

  ExecPrmScoring moduleScoring(scoringParams,
                               &inputSpectra,
                               &outputSpectra);

  string errorString;
  bool isValid = moduleScoring.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  bool returnStatus = moduleScoring.invoke();
  // Test for return status
  TEST_RETURN_STATUS("moduleScoring");

  DEBUG_TRACE;

  moduleScoring.saveOutputData();

  return true;
}
//-----------------------------------------------------------------------------
bool performGenoMS(ParameterList & ip)
{
  DEBUG_TRACE;

  ParameterList genoMSParams;

  genoMSParams.addIfExists(ip,"EXE_DIR");
  genoMSParams.addIfExists(ip,"TOLERANCE_PEAK");
  genoMSParams.addIfExists(ip,"TOLERANCE_PM");
  genoMSParams.addIfExists(ip,"PEAK_PENALTY");

  genoMSParams.addIfExists(ip,"ALL2ALL_SIMILARITY");
  genoMSParams.addIfExists(ip,"HMM_LATE_ADD");
  genoMSParams.addIfExists(ip,"FDR_CUTOFF");
  genoMSParams.addIfExists(ip,"MUTATION_MODE");
  genoMSParams.addIfExists(ip,"DBCOMBINED");
  genoMSParams.addIfExists(ip,"DBROOTNAME");
  genoMSParams.addIfExists(ip,"GENOMESEQ");
  genoMSParams.addIfExists(ip,"DIGEST");
  genoMSParams.addIfExists(ip,"TEMPLATECONSTRAINTFILE");
  genoMSParams.addIfExists(ip,"FIXEDMOD");
  genoMSParams.addIfExists(ip,"PROJECT_DIR");
  genoMSParams.addIfExists(ip,"RUN_DBSEARCH");


  genoMSParams.setValue("PRMS",     "spectra/specs_scored.prms");
  genoMSParams.setValue("SPECTRA",     "spectra/specs_ms_z.mgf");
  genoMSParams.setValue("OUTPUT_FILE","genoMS.out");
  genoMSParams.setValue("GENERATE_REPORTS","1");
  genoMSParams.setValue("LOG_FILE","genoMS.log");

  DEBUG_TRACE;
  stringstream aux;
  genoMSParams.print(aux);
  DEBUG_MSG(aux.str());
  DEBUG_TRACE;
  ExecGenoMS moduleGenoMS(genoMSParams);

  string errorString;
  bool isValid = moduleGenoMS.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  bool returnStatus = moduleGenoMS.invoke();
  // Test for return status
  TEST_RETURN_STATUS("moduleGenoMS");

  DEBUG_TRACE;

  return true;
}

bool performMergeOfCSPSAndGenoMS(ParameterList & ip)
{
  ///////////////////////////////////////////////
  //First we duplicate the input spectrum file
  ///////////////////////////////////////////////
  bool debug = false;
  DEBUG_TRACE;

  unsigned int cspsContigCount;
  unsigned int genomsContigCount;

  spsReports::ClusterData clusterData;


  SpecSet inputSpectra;
  SpecSet other;
  unsigned int specCount; // total (unique spectra)

  if(inputSpectra.LoadSpecSet_pklbin("./spectra/specs_ms.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/specs_ms.pklbin");
      return false;
    }

  DEBUG_MSG("First we have " << inputSpectra.size() << " inputSpectra!");
  if(other.LoadSpecSet_pklbin("./spectra/specs_ms.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/specs_ms.pklbin");
      return false;
    }



  inputSpectra.appendSpecSet(other);

  specCount = inputSpectra.size();

  DEBUG_MSG("Now we have " << specCount << " inputSpectra!");
  if(inputSpectra.SaveSpecSet_pklbin("./spectra/specs_ms.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving doubled spectra to ./spectra/specs_ms.pklbin");
      return false;
    }
  if(debug)
    {
      inputSpectra.SaveSpecSet_mgf("./spectra/specs_ms.mgf");
    }

  ///////////////////////////////////////////////
  //Also duplicate the clusterData.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Doubling clusterData.bin");
  clusterData.loadDataDouble(".");
  clusterData.saveData("..");


  ///////////////////////////////////////////////
  //Merge the star spectra
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging star spectra, csps + genoms");
    if(inputSpectra.LoadSpecSet_pklbin("./spectra/csps.stars.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/csps.stars.pklbin");
      return false;
    }

  if(other.LoadSpecSet_pklbin("./spectra/genoms.stars.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/genoms.stars.pklbin");
      return false;
    }

  inputSpectra.appendSpecSet(other);
  if(inputSpectra.SaveSpecSet_pklbin("./spectra/stars.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving doubled star spectra to ./spectra/stars.pklbin");
      return false;
    }

  DEBUG_MSG("Made " << inputSpectra.size() << " star spectra!");
    if(debug)
    {
      inputSpectra.SaveSpecSet_mgf("./spectra/stars.mgf");
    }
  ///////////////////////////////////////////////
  //Merge the contig spectra
  //////////////////////////////////////////////
  DEBUG_MSG("Merging contig spectra");
  if(inputSpectra.LoadSpecSet_pklbin("./assembly/csps.sps_seqs.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./assembly/csps.sps_seqs.pklbin");
      return false;
    }
  cspsContigCount = inputSpectra.size();

  if(other.LoadSpecSet_pklbin("./assembly/genoms.sps_seqs.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./assembly/genoms.sps_seqs.pklbin");
      return false;
    }
  genomsContigCount = other.size();
  inputSpectra.appendSpecSet(other);
  if(inputSpectra.SaveSpecSet_pklbin("./assembly/sps_seqs.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving merged contig spectra to ./assembly/sps_seqs.pklbin");
      return false;
    }
  DEBUG_MSG("Made " << inputSpectra.size() << " contig spectra!");
  DEBUG_MSG("CSPS has contigs 0-"<< (cspsContigCount-1));
  DEBUG_MSG("GenoMS has contigs "<<cspsContigCount << "-" << (genomsContigCount + cspsContigCount-1));

  if(debug)
    {
      inputSpectra.SaveSpecSet_mgf("./assembly/sps_seqs.mgf");
    }
  ///////////////////////////////////////////////
  //Merge abruin infos
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging abruijn info");
  abinfo_t info1;
  abinfo_t info2;
  string abinfoOutFile = "./assembly/component_info.bin";

  Load_abinfo("./assembly/csps.component_info.bin",info1);
  Load_abinfo("./assembly/genoms.component_info.bin",info2);

  if(!merge_abinfo(info1,info2,abinfoOutFile.c_str(),specCount/2))
    {
      ERROR_MSG("Error merging abruijn info");
      return false;
    }
  if(debug)
    {
      abinfo_t info;
      Load_abinfo("./assembly/component_info.bin",info);
      dumpAbInfo("./assembly/component_info.txt",info);
    }
  ///////////////////////////////////////////////
  //merge contig_mp.bin
  ///////////////////////////////////////////////

  DEBUG_MSG("Merging contig_mp_all.bin");
  vector<vector<int> > * m_contigs_mp_sub1;
  vector<vector<int> > * m_contigs_mp_sub2;

  m_contigs_mp_sub1 = new vector<vector<int> >();
  m_contigs_mp_sub2 = new vector<vector<int> >();

  DEBUG_MSG("Initialized vectors..");
  // load the data
  if(Load_binArray<int,vector>("./homology/csps.contigs_mp.bin", *m_contigs_mp_sub1) < 0) {

    return false;
  }

  if(Load_binArray<int,vector>("./homology/genoms.contigs_mp.bin", *m_contigs_mp_sub2) < 0) {

    return false;
  }
  DEBUG_MSG("Loaded contig_mp_all vectors");
  if(Save_doubleBinArray<int>("./homology/contigs_mp.bin",*m_contigs_mp_sub1,*m_contigs_mp_sub2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./homology/contigs_mp.bin");
      return false;

    }

DEBUG_MSG("CSPS contig_mp.bin = " << (*m_contigs_mp_sub1).size() << "x" << (*m_contigs_mp_sub1)[0].size());
  DEBUG_MSG("GenoMS contig_mp.bin = " << (*m_contigs_mp_sub2).size() << "x" << (*m_contigs_mp_sub2)[0].size());
  if(debug)
    {

      cout << "DEBUG: homology/contigs_mp.bin" << endl;

      Load_binArray<int,vector>("./homology/contigs_mp.bin",*m_contigs_mp_sub1);
      for(int i = 0; i < (*m_contigs_mp_sub1).size(); ++i)
	{
	  for(int j = 0; j < (*m_contigs_mp_sub1)[0].size(); ++j)
	    {
	      cout << " " << (*m_contigs_mp_sub1)[i][j];
	    }
	  cout << endl;
	}


    }
///////////////////////////////////////////////
  //Merge contigs_midx.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contigs_midx.pklbin");
  if(inputSpectra.LoadSpecSet_pklbin("./homology/csps.contigs_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.contigs_midx.pklbin");
      return false;
    }

  if(other.LoadSpecSet_pklbin("./homology/genoms.contigs_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.contigs_midx.pklbin");
      return false;
    }

  inputSpectra.appendSpecSet(other);
  if(inputSpectra.SaveSpecSet_pklbin("./homology/contigs_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/contigs_midx.pklbin");
      return false;
    }
  if(debug)
    {
      inputSpectra.SaveSpecSet_mgf("./homology/contigs_midx.mgf");
    }

  ///////////////////////////////////////////////
  //merge contig_mp_all.bin
  ///////////////////////////////////////////////

  DEBUG_MSG("Merging contig_mp_all.bin");
  vector<vector<int> > * m_contigs_mp1;
  vector<vector<int> > * m_contigs_mp2;

  m_contigs_mp1 = new vector<vector<int> >();
  m_contigs_mp2 = new vector<vector<int> >();

  DEBUG_MSG("Initialized vectors..");
  // load the data
  if(Load_binArray<int,vector>("./homology/csps.contigs_mp_all.bin", *m_contigs_mp1) < 0) {

    return false;
  }

  if(Load_binArray<int,vector>("./homology/genoms.contigs_mp_all.bin", *m_contigs_mp2) < 0) {

    return false;
  }
  DEBUG_MSG("Loaded contig_mp_all vectors");
  if(Save_doubleBinArray<int>("./homology/contigs_mp_all.bin",*m_contigs_mp1,*m_contigs_mp2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./homology/contigs_mp_all.bin");
      return false;

    }

DEBUG_MSG("CSPS contig_mp_all.bin = " << (*m_contigs_mp1).size() << "x" << (*m_contigs_mp1)[0].size());
  DEBUG_MSG("GenoMS contig_mp_all.bin = " << (*m_contigs_mp2).size() << "x" << (*m_contigs_mp2)[0].size());
  if(debug)
    {

      cout << "DEBUG: homology/contigs_mp_all.bin" << endl;

      Load_binArray<int,vector>("./homology/contigs_mp_all.bin",*m_contigs_mp1);
      for(int i = 0; i < (*m_contigs_mp1).size(); ++i)
	{
	  for(int j = 0; j < (*m_contigs_mp1)[0].size(); ++j)
	    {
	      cout << " " << (*m_contigs_mp1)[i][j];
	    }
	  cout << endl;
	}


    }
  ///////////////////////////////////////////////
  //Merge contigs_midx_all.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contigs_midx_all.pklbin");
  if(inputSpectra.LoadSpecSet_pklbin("./homology/csps.contigs_midx_all.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.contigs_midx_all.pklbin");
      return false;
    }

  if(other.LoadSpecSet_pklbin("./homology/genoms.contigs_midx_all.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.contigs_midx_all.pklbin");
      return false;
    }

  inputSpectra.appendSpecSet(other);
  if(inputSpectra.SaveSpecSet_pklbin("./homology/contigs_midx_all.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/contigs_midx_all.pklbin");
      return false;
    }
  if(debug)
    {
      inputSpectra.SaveSpecSet_mgf("./homology/contigs_midx_all.mgf");
    }
  ///////////////////////////////////////////////
  //Merge contigs.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contigs.pklbin");
  if(inputSpectra.LoadSpecSet_pklbin("./spectra/csps.contigs.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/csps.contigs.pklbin");
      return false;
    }
  unsigned int cspsMappedContigCount = inputSpectra.size();
  if(other.LoadSpecSet_pklbin("./spectra/genoms.contigs.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./spectra/genoms.contigs.pklbin");
      return false;
    }

  inputSpectra.appendSpecSet(other);
  if(inputSpectra.SaveSpecSet_pklbin("./spectra/contigs.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig spectra to ./spectra/contigs.pklbin");
      return false;
    }
  DEBUG_MSG("Made " << inputSpectra.size() << " mapped contig spectra!");
  DEBUG_MSG("CSPS has contigs 0-"<< (cspsMappedContigCount-1));
  DEBUG_MSG("GenoMS has contigs "<<cspsMappedContigCount << "-" << (inputSpectra.size()-1));
  if(debug)
    {
      inputSpectra.SaveSpecSet_mgf("./spectra/contigs.mgf");
    }
  ///////////////////////////////////////////////
  //merge contig_indices.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging contig_indices.bin");

  vector<vector<int> > * m_contigIndex1;
  vector<vector<int> > * m_contigIndex2;

  m_contigIndex1 = new vector<vector<int> >();
  m_contigIndex2 = new vector<vector<int> >();

  // load the data
  if(Load_binArray<int,vector>("./spectra/csps.contigs_indices.bin", *m_contigIndex1) < 0) {

    ERROR_MSG("Error loading binary array ./spectra/csps.contigs_indices.bin");
    return false;
  }

  if(Load_binArray<int,vector>("./spectra/genoms.contigs_indices.bin", *m_contigIndex2) < 0) {

    ERROR_MSG("Error loading binary array ./spectra/genoms.contigs_indices.bin");
    return false;
  }
  DEBUG_MSG("CSPS contig_indices.bin = " << (*m_contigIndex1).size() << "x" << (*m_contigIndex1)[0].size());
  DEBUG_MSG("GenoMS contig_indices.bin = " << (*m_contigIndex2).size() << "x" << (*m_contigIndex2)[0].size());


  //Update component indices from genoms result
  for(int i = 0; i < m_contigIndex2->size(); ++i)
    {
      (*m_contigIndex2)[i][0] += cspsContigCount;

    }

  if(Save_doubleBinArray<int>("./spectra/contigs_indices.bin",*m_contigIndex1,*m_contigIndex2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./spectra/contigs_indices.bin");
      return false;

    }

  if(debug)
    {

      cout << "DEBUG: spectra/contigs_indices.bin" << endl;

      Load_binArray<int,vector>("./spectra/contigs_indices.bin",*m_contigIndex1);
      for(int i = 0; i < (*m_contigIndex1).size(); ++i)
	{
	  for(int j = 0; j < (*m_contigIndex1)[0].size(); ++j)
	    {
	      cout << " " << (*m_contigIndex1)[i][j];
	    }
	  cout << endl;
	}


    }


  ///////////////////////////////////////////////
  //merge homglue_ref_mp.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_ref_mp.bin");
  vector<vector<int> > * m_homglueRefMp1;
  vector<vector<int> > * m_homglueRefMp2;

  m_homglueRefMp1 = new vector<vector<int> >();
  m_homglueRefMp2 = new vector<vector<int> >();

  // load the data
  if(Load_binArray<int,vector>("./homology/csps.homglue_ref_mp.bin", *m_homglueRefMp1) < 0) {

    return false;
  }

  if(Load_binArray<int,vector>("./homology/genoms.homglue_ref_mp.bin", *m_homglueRefMp2) < 0) {

    return false;
  }

    DEBUG_MSG("CSPS homglue_ref_mp.bin = " << (*m_homglueRefMp1).size() << "x" << (*m_homglueRefMp1)[0].size());
  DEBUG_MSG("GenoMS homglue_ref_mp.bin = " << (*m_homglueRefMp2).size() << "x" << (*m_homglueRefMp2)[0].size());

  if(Save_doubleBinArray("./homology/homglue_ref_mp.bin",*m_homglueRefMp1,*m_homglueRefMp2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./homology/homglue_ref_mp.bin");
      return false;

    }

  if(debug)
    {

      cout << "DEBUG: homology/homglue_ref_mp.bin" << endl;

      Load_binArray<int,vector>("./homology/homglue_ref_mp.bin",*m_homglueRefMp1);
      for(int i = 0; i < (*m_homglueRefMp1).size(); ++i)
	{
	  for(int j = 0; j < (*m_homglueRefMp1)[0].size(); ++j)
	    {
	      cout << " " << (*m_homglueRefMp1)[i][j];
	    }
	  cout << endl;
	}


    }

  ///////////////////////////////////////////////
  //Merge homglue_ref_midx.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_ref_midx.pklbin");
  if(inputSpectra.LoadSpecSet_pklbin("./homology/csps.homglue_ref_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_ref_midx.pklbin");
      return false;
    }

  if(other.LoadSpecSet_pklbin("./homology/genoms.homglue_ref_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_ref_midx.pklbin");
      return false;
    }

  inputSpectra.appendSpecSet(other);
  if(inputSpectra.SaveSpecSet_pklbin("./homology/homglue_ref_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_ref_midx.pklbin");
      return false;
    }

  if(debug)
    {
      DEBUG_MSG("homglue_ref_midx.pklbin has" << inputSpectra.size() << " spectra");
      inputSpectra.SaveSpecSet_mgf("./homology/homglue_ref_midx.mgf");
    }

  ///////////////////////////////////////////////
  //Merge homglue_matches.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_matches.pklbin");
  if(inputSpectra.LoadSpecSet_pklbin("./homology/csps.homglue_matches.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_matches.pklbin");
      return false;
    }

  if(other.LoadSpecSet_pklbin("./homology/genoms.homglue_matches.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_matches.pklbin");
      return false;
    }

  inputSpectra.appendSpecSet(other);
  if(inputSpectra.SaveSpecSet_pklbin("./homology/homglue_matches.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_matches.pklbin");
      return false;
    }
  if(debug)
    {
      DEBUG_MSG("homglue_matches.pklbin has" << inputSpectra.size() << " spectra");
      inputSpectra.SaveSpecSet_mgf("./homology/homglue_matches.mgf");
    }

  ///////////////////////////////////////////////
  //merge homglue_matches_mp.bin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_matches_mp.bin");
  vector<vector<int> > * m_homglueMatchMp1;
  vector<vector<int> > * m_homglueMatchMp2;

  m_homglueMatchMp1 = new vector<vector<int> >();
  m_homglueMatchMp2 = new vector<vector<int> >();

  // load the data
  if(Load_binArray<int,vector>("./homology/csps.homglue_matches_mp.bin", *m_homglueMatchMp1) < 0) {

    return false;
  }

  if(Load_binArray<int,vector>("./homology/genoms.homglue_matches_mp.bin", *m_homglueMatchMp2) < 0) {

    return false;
  }

    DEBUG_MSG("CSPS homglue_matches_mp.bin = " << (*m_homglueMatchMp1).size() << "x" << (*m_homglueMatchMp1)[0].size());
  DEBUG_MSG("GenoMS homglue_matches_mp.bin = " << (*m_homglueMatchMp2).size() << "x" << (*m_homglueMatchMp2)[0].size());

  if(Save_doubleBinArray("./homology/homglue_matches_mp.bin",*m_homglueMatchMp1,*m_homglueMatchMp2) < 0)
    {
      ERROR_MSG("Error merging arrays for ./homology/homglue_matches_mp.bin");
      return false;

    }

  ///////////////////////////////////////////////
  //Merge homglue_matches_midx.pklbin
  ///////////////////////////////////////////////
  DEBUG_MSG("Merging homglue_matches.pklbin");
  if(inputSpectra.LoadSpecSet_pklbin("./homology/csps.homglue_matches_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/csps.homglue_matches_midx.pklbin");
      return false;
    }

  if(other.LoadSpecSet_pklbin("./homology/genoms.homglue_matches_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem loading spectra from ./homology/genoms.homglue_matches_midx.pklbin");
      return false;
    }

  inputSpectra.appendSpecSet(other);
  if(inputSpectra.SaveSpecSet_pklbin("./homology/homglue_matches_midx.pklbin") <= 0)
    {
      ERROR_MSG("Problem saving contig index spectra to ./homology/homglue_matches_midx.pklbin");
      return false;
    }

  if(debug)
    {
      DEBUG_MSG("homglue_matches_midx.pklbin has" << inputSpectra.size() << " spectra");
      inputSpectra.SaveSpecSet_mgf("./homology/homglue_matches_midx.mgf");
    }

  ///////////////////////////////////////////////
  //merge ref_sps_names.txt
  ///////////////////////////////////////////////

  DEBUG_MSG("Merging ref_sps_names.txt");
  vector<string> genoMScontigNames;
  int i,j;
  int contigNum = cspsContigCount + 1;
  //The GenoMS contig numbers need to be updated
  readFilesFromFile("./homology/genoms.ref_sps_names.txt",genoMScontigNames);


  FILE * fp = fopen("./homology/genoms.ref_sps_names.txt","w");
  //We assume the names are GenoMS:i, where i is the contig number
  for(i = 0; i < genoMScontigNames.size(); ++i)
    {
      DEBUG_MSG("Orig genoMS contig Name: " << genoMScontigNames[i]);
      string newName = "";
      for(j = 0; j < genoMScontigNames[i].size(); ++j)
	{
	  if (genoMScontigNames[i][j] == ':')
	    break;
	  else
	    newName += genoMScontigNames[i][j];
	}
      newName += ":";
      newName += intToString(contigNum);
      newName += "\n";
      genoMScontigNames[i] = newName;
      DEBUG_MSG("New genoMS contig Name: " << genoMScontigNames[i]);

      fwrite((void *)(genoMScontigNames[i].c_str()),genoMScontigNames[i].size(),1,fp);
      contigNum += 1;

    }

  fclose(fp);


  if(!concatenateFiles("./homology/csps.ref_sps_names.txt","./homology/genoms.ref_sps_names.txt","./homology/ref_sps_names.txt"))
    return false;



  return true;

}

bool generateRelaunchScript(ParameterList & ip)
{
  string relauncherScriptName = "relauncher.sh";
  FILE* myfile;
  myfile = fopen(relauncherScriptName.c_str(),"r");
  if(myfile)
    {
      fclose(myfile);
      return true;
    }



  string dbFileName = ip.getValue("PROJECT_DIR",".") + "/relaunch.protid.fasta";
  string paramsFile = "relaunch.params";
  string exeDir = ip.getValue("EXE_DIR",".");



  string logLevel = ip.getValue("LOG_LEVEL","0");
  if(strlen(logLevel.c_str()) == 0)
    logLevel = "0";
  //Create the command

  //Change to the project directory
  string relauncherScript = "cd " + ip.getValue("PROJECT_DIR", ".") + "\n";

  relauncherScript += "echo \"Running\" > status.txt\n";

  //If we were currently running in merge mode, then we need to move the old CSPS files back to their old places
  // if(ip.getValueInt("MERGE_FLAG",0) == 1)
  myfile = fopen("assembly/csps.sps_seqs.pklbin","r");
  if(myfile)
    {
      fclose(myfile);
      relauncherScript += "mv assembly/csps.sps_seqs.pklbin assembly/sps_seqs.pklbin\n";
      relauncherScript += "mv assembly/csps.component_info.bin assembly/component_info.bin\n";
      relauncherScript += "mv homology/csps.contigs_midx_all.pklbin homology/contigs_midx_all.pklbin\n";
      relauncherScript += "mv homology/csps.contigs_midx.pklbin homology/contigs_midx.pklbin\n";
      relauncherScript += "mv homology/csps.contigs_mp_all.bin homology/contigs_mp_all.bin\n";
      relauncherScript += "mv homology/csps.contigs_mp.bin homology/contigs_mp.bin\n";
      relauncherScript += "mv homology/csps.homglue_matches_midx.pklbin homology/homglue_matches_midx.pklbin\n";
      relauncherScript += "mv homology/csps.homglue_matches_mp.bin homology/homglue_matches_mp.bin\n";
      relauncherScript += "mv homology/csps.homglue_matches.pklbin homology/homglue_matches.pklbin\n";
      relauncherScript += "mv homology/csps.homglue_ref_midx.pklbin homology/homglue_ref_midx.pklbin\n";
      relauncherScript += "mv homology/csps.homglue_ref_mp.bin homology/homglue_ref_mp.bin\n";
      relauncherScript += "mv homology/csps.ref_sps_names.txt homology/ref_sps_names.txt\n";
      relauncherScript += "mv spectra/csps.contigs_indices.bin spectra/contigs_indices.bin\n";
      relauncherScript += "mv spectra/csps.contigs.pklbin spectra/contigs.pklbin\n";
      relauncherScript += "mv spectra/csps.stars.pklbin spectra/stars.pklbin\n";
    }


  //execute the command
  relauncherScript += exeDir + "/main_specnets " + paramsFile;
  relauncherScript += " -i tagsearch";
  relauncherScript += " -ll " + logLevel;
  relauncherScript += " -lf " + ip.getValue("LOG_FILE_NAME", "relauncher_log.txt");
  int val = ip.getValueInt("GRID_EXECUTION");
  if(val != 0)
    relauncherScript += " -g";
  //val = ip.getValueInt("GENOMS_FLAG");
  //if(val != 0)
  //  relauncherScript += " -q";
  //val = ip.getValueInt("MERGE_FLAG");
  //if(val != 0)
  //  relauncherScript += " -m";
  relauncherScript += "\n";
  //Write the relauncher.sh
  FILE * f = fopen(relauncherScriptName.c_str(),"w");

  if(f == NULL)
    {
      return false;
    }

  val = fwrite(relauncherScript.c_str(),sizeof(char),strlen(relauncherScript.c_str()),f);
  if(val != strlen(relauncherScript.c_str()))
    {
      ERROR_MSG("Problem encountered writing relauncher file");
      ERROR_MSG(val);
      ERROR_MSG(strlen(relauncherScript.c_str()));
      return false;
    }
  fclose(f);

    //Write the updated params file
  ip.setValue("FASTA_DATABASE", dbFileName);
  ip.setValue("GENOMS_FLAG","0");
  ip.setValue("MERGE_FLAG","0");
  ip.writeToFile(getProjPath(ip,".") + "/" + paramsFile);


  DEBUG_MSG("Wrote relauncher.sh");
}


//-----------------------------------------------------------------------------
bool performFilterPairs(ParameterList & ip,
                        SpecSet & inputSpectra,
                        SpecSet & inputSpectraMS2,
                        SpectrumPairSet & filteredPairs,
                        vector<TwoValues<float> > & ratios,
                        vector<TwoValues<float> > & means,
                        vector<float> & varTerms,
                        list<vector<float> > & alignStats,
                        vector<vector<float> > & specStats,
                        std::vector<unsigned int> & idxKept,
                        std::vector<TwoValues<float> > & pvalues,
                        bool gridExecutionFlag,
                        bool resume)
{
  ParameterList filterpairsParams;
  filterpairsParams.addIfExists(ip,"TOLERANCE_PEAK");
  filterpairsParams.addIfExists(ip,"TOLERANCE_PM");
  filterpairsParams.addIfExists(ip,"PARTIAL_OVERLAPS");
  filterpairsParams.addIfExists(ip,"AA_DIFF_COUNT");
  filterpairsParams.addIfExists(ip,"MIN_OVERLAP_AREA");
  filterpairsParams.addIfExists(ip,"MIN_RATIO");
  filterpairsParams.addIfExists(ip,"MIN_SHIFT");
  filterpairsParams.addIfExists(ip,"MAX_SHIFT");
  filterpairsParams.addIfExists(ip,"PAIRS_MATCH_MODE");
  filterpairsParams.addIfExists(ip,"PAIRS_MIN_COSINE");

  filterpairsParams.addIfExists(ip,"GRID_TYPE");
  filterpairsParams.addIfExists(ip,"GRID_NUMNODES");
  filterpairsParams.addIfExists(ip,"GRID_NUMCPUS");
  filterpairsParams.addIfExists(ip,"GRID_EXE_DIR");
  filterpairsParams.addIfExists(ip,"GRID_SGE_EXE_DIR");
  filterpairsParams.addIfExists(ip,"GRID_PARAMS");

  filterpairsParams.addIfExists(ip,"MAX_PVALUE");
  filterpairsParams.addIfExists(ip,"RESOLUTION");
  filterpairsParams.addIfExists(ip,"TAGS_MATCH_FLANK");
  filterpairsParams.addIfExists(ip,"TAGS_MATCH_COUNT");

  filterpairsParams.addIfExists(ip, "PROJECT_DIR");

  filterpairsParams.addIfDoesntExist("USE_MIN_DIST_57","1");
  filterpairsParams.addIfDoesntExist("SPEC_TYPE_MSMS","0");

  filterpairsParams.setValue("MIN_NUM_MATCHED_PEAKS", ip.getValue("MIN_MATCHED_PEAKS"));
  filterpairsParams.setValue("MAX_SHIFT", ip.getValue("MAX_MOD_MASS"));

  filterpairsParams.setValue("INPUT_SPECS_PKLBIN", getProjPath(ip, "./spectra/specs_scored.pklbin") );
  filterpairsParams.setValue("GRID_DATA_DIR",      getProjPath(ip, "./aligns") );
  filterpairsParams.setValue("OUTPUT_ALIGNS",      getProjPath(ip, "./aligns/pairs_raw.bin") );
  filterpairsParams.setValue("OUTPUT_MEANS",       getProjPath(ip, "./aligns/means.bin") );
  filterpairsParams.setValue("OUTPUT_VARIANCE",    getProjPath(ip, "./aligns/vars.bin") );
  filterpairsParams.setValue("OUTPUT_RATIOS",      getProjPath(ip, "./aligns/ratios.bin") );

  filterpairsParams.writeToFile("debug_filterpairs.params");

  DEBUG_TRACE;

  stringstream aux;
  filterpairsParams.print(aux);
  DEBUG_MSG(aux.str());

  DEBUG_TRACE;

  if(filterpairsParams.getValue("PAIRS_MATCH_MODE","") == "cosine") {
  	for(unsigned int i=0; i<inputSpectraMS2.size(); i++) {
  		for(unsigned int j=0; j<inputSpectraMS2[i].size(); j++)
  			inputSpectraMS2[i][j][1] = sqrt(inputSpectraMS2[i][j][1]);
  		inputSpectraMS2[i].normalize2();
  	}
  }

  ExecFilterPairs moduleFilterPairs(filterpairsParams,
                                    &inputSpectra,
                                    &inputSpectraMS2,
                                    &filteredPairs,
                                    &ratios,
                                    &means,
                                    &varTerms,
                                    &alignStats,
                                    &specStats,
                                    &idxKept,
                                    &pvalues);

  string errorString;
  bool isValid = moduleFilterPairs.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  bool returnStatus;
  if(!ip.exists("GRID_NUMNODES") or ip.getValueInt("GRID_NUMNODES")<=0)
  {
    returnStatus = moduleFilterPairs.invoke();
  }
  else
  {
  	DEBUG_TRACE;
    int numNodes = ip.getValueInt("GRID_NUMNODES");

    string gridType = ip.getValue("GRID_TYPE");
    if (gridType == "pbs") {
    	ParallelPbsExecution exec(&moduleFilterPairs, gridExecutionFlag, !gridExecutionFlag, resume);
    	returnStatus = exec.invoke(numNodes);
    } else if (gridType == "sge") {
     	ParallelSgeExecution exec(&moduleFilterPairs, gridExecutionFlag, !gridExecutionFlag, resume);
    	returnStatus = exec.invoke(numNodes);
    }
  }


  // Test for return status
  TEST_RETURN_STATUS("moduleFilterPairs");

  DEBUG_TRACE;

  moduleFilterPairs.saveOutputData();

  return true;
}

//-----------------------------------------------------------------------------
bool performFilterAligns(ParameterList & ip,
                         SpectrumPairSet & filteredPairs,
                         vector<TwoValues<float> > & ratios,
                         vector<TwoValues<float> > & means,
                         vector<float> & varTerms,
                         std::vector<unsigned int> & idxKept,
                         std::vector<TwoValues<float> > & pvalues)
{
  ParameterList filteralignsParams;
  filteralignsParams.addIfExists(ip,"TOLERANCE_PM");
  filteralignsParams.addIfExists(ip,"MIN_RATIO");
  filteralignsParams.addIfExists(ip,"MAX_PVALUE");
  filteralignsParams.addIfExists(ip,"FILTER_TRIGS");

  filteralignsParams.setValue("INPUT_ALIGNS",   getProjPath(ip, "./aligns/pairs_raw.bin") );
  filteralignsParams.setValue("INPUT_MEANS",    getProjPath(ip, "./aligns/means.bin") );
  filteralignsParams.setValue("INPUT_VARIANCE", getProjPath(ip, "./aligns/vars.bin") );
  filteralignsParams.setValue("INPUT_RATIOS",   getProjPath(ip, "./aligns/ratios.bin") );
  filteralignsParams.setValue("OUTPUT_ALIGNS",  getProjPath(ip, "./aligns/pairs.bin") );

  filteralignsParams.writeToFile("debug_filteraligns.params");

  DEBUG_TRACE;

  stringstream aux;
  filteralignsParams.print(aux);
  DEBUG_MSG(aux.str());

  DEBUG_TRACE;

  ExecFilterAligns moduleFilterAligns(filteralignsParams,
                                      &filteredPairs,
                                      &ratios,
                                      &means,
                                      &varTerms,
                                      &idxKept,
                                      &pvalues);

  string errorString;
  bool isValid = moduleFilterAligns.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

  bool returnStatus = moduleFilterAligns.invoke();

  // Test for return status
  TEST_RETURN_STATUS("moduleFilterAligns");

  DEBUG_TRACE;

  // Save output data
  moduleFilterAligns.saveOutputData();

  return true;
}

//-----------------------------------------------------------------------------
bool performAlignment(ParameterList & ip,
                      AAJumps & jumps,
                      SpecSet &inputSpectra,
                      SpectrumPairSet &inputPairs,
                      SpecSet &pairAlignments,
                      SpecSet &starSpectraOnly,
                      SpecSet &starSpectra,
                      vector<unsigned int> &alignedSpectra)
{
  ParameterList alignParams;
  alignParams.addIfExists(ip,"TOLERANCE_PEAK");
  alignParams.addIfExists(ip,"TOLERANCE_PM");
  alignParams.addIfExists(ip,"RESOLUTION");
  alignParams.addIfExists(ip,"PARTIAL_OVERLAPS");
  alignParams.addIfExists(ip,"MAX_AA_JUMP");
  alignParams.addIfExists(ip,"PENALTY_PTM");
  alignParams.addIfExists(ip,"PENALTY_SAME_VERTEX");

  alignParams.addIfExists(ip, "PROJECT_DIR");

  alignParams.setValue("PENALTY_PTM_PEAKS",   "-1.0");
  alignParams.setValue("OUTPUT_STARS",       getProjPath(ip, "./spectra/stars_only.pklbin") );
  alignParams.setValue("OUTPUT_STARS_INDEX", getProjPath(ip, "./spectra/stars_indices.bin") );
  alignParams.setValue("OUTPUT_STARS_ALL",   getProjPath(ip, "./spectra/stars.pklbin") );
  if (ip.exists("AMINO_ACID_MASSES"))
  {
    alignParams.setValue("AMINO_ACID_MASSES", "../"
        + ip.getValue("AMINO_ACID_MASSES"));
  }

  alignParams.writeToFile("debug_align.params");

  DEBUG_TRACE;
  ExecAlignment moduleAlignment(alignParams,
                                &inputSpectra,
                                &inputPairs,
                                &pairAlignments,
                                &starSpectraOnly,
                                &starSpectra,
                                &alignedSpectra);

  string errorString;
  bool isValid = moduleAlignment.validateParams(errorString);
  TEST_VALID;

  DEBUG_TRACE;

#if 1
  bool returnStatus = moduleAlignment.invoke();
#else
  DEBUG_TRACE;
  ParallelThreadedExecution exec(&moduleAlignment);
  bool returnStatus = exec.invoke(1, 1);
#endif

  // Test for return status
  TEST_RETURN_STATUS("ExecAlignment");

  DEBUG_TRACE;
  // Save output data
  moduleAlignment.saveOutputData();

  return true;
}

//-----------------------------------------------------------------------------
bool performFilterStarPairs(ParameterList & ip,
                            SpectrumPairSet &inputPairs,
                            SpecSet &starSpectra,
                            vector<vector<float> > &ratios,
                            SpecSet &matchedPeaks)
{
  ParameterList filterstarpairsParams;
  filterstarpairsParams.addIfExists(ip,"TOLERANCE_PEAK");
  filterstarpairsParams.addIfExists(ip,"TOLERANCE_PM");
  filterstarpairsParams.addIfExists(ip,"MAX_MOD_MASS");
  filterstarpairsParams.addIfExists(ip,"PARTIAL_OVERLAPS");
  filterstarpairsParams.addIfExists(ip,"MIN_RATIO");
  filterstarpairsParams.addIfExists(ip,"MIN_MATCHED_PEAKS");
  filterstarpairsParams.addIfExists(ip,"MAX_AA_JUMP");
  filterstarpairsParams.addIfExists(ip,"PENALTY_PTM");
  filterstarpairsParams.addIfExists(ip,"PENALTY_SAME_VERTEX");

  filterstarpairsParams.addIfExists(ip, "PROJECT_DIR");

  filterstarpairsParams.setValue("SPEC_TYPE_MSMS", "0");
  filterstarpairsParams.setValue("OUTPUT_ALIGNS", getProjPath(ip, "aligns/pairs_stars.bin") );

  filterstarpairsParams.writeToFile(getProjPath(ip, "debug_filterstarpairs.params"));

  DEBUG_TRACE;
  ExecFilterStarPairs moduleFilterStarPairs(filterstarpairsParams,
                                            &inputPairs,
                                            &starSpectra,
                                            &ratios,
                                            &matchedPeaks);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleFilterStarPairs.validateParams(errorString);
  TEST_VALID;

#if 1
  bool returnStatus = moduleFilterStarPairs.invoke();
#else
  ParallelSgeExecution exec(&moduleFilterStarPairs);
  //ParallelThreadedExecution exec(&moduleFilterStarPairs);
  bool returnStatus = exec.invoke(1, 1);
#endif

  // Test for return status
  TEST_RETURN_STATUS("ExecFilterStarPairs");

  // Save output data
  returnStatus = moduleFilterStarPairs.saveOutputData();
//  TEST_SAVE_OUPUT_DATA("ExecFilterStarPairs");

  return true;
}

//-----------------------------------------------------------------------------

bool performAssembly(ParameterList & ip,
                     SpecSet & starSpectra,
                     SpectrumPairSet & starPairs,
                     Clusters & contigShifts,
                     abinfo_t & contigAbinfo)
{
  ParameterList assemblyParams;
  assemblyParams.addIfExists(ip,"TOLERANCE_PEAK");
  assemblyParams.addIfExists(ip,"TOLERANCE_PM");
  assemblyParams.addIfExists(ip,"MAX_MOD_MASS");
  assemblyParams.addIfExists(ip,"MAX_AA_JUMP");
  assemblyParams.addIfExists(ip,"PENALTY_PTM");
  assemblyParams.addIfExists(ip,"MIN_MATCHED_PEAKS");
  assemblyParams.addIfExists(ip,"PENALTY_SAME_VERTEX");
  assemblyParams.addIfExists(ip,"PENALTY_PTM");

  assemblyParams.addIfExists(ip, "PROJECT_DIR");

  assemblyParams.setValue("NO_SEQUENCING", "false");
  assemblyParams.setValue("ADD_ENDPOINTS", "true");
  assemblyParams.setValue("EDGE_SCORE_TYPE", "1");
  assemblyParams.setValue("SPEC_TYPE_MSMS", "0");
  assemblyParams.setValue("GRAPH_TYPE", "2");
  assemblyParams.setValue("OUTPUT_COMPLETE_ABRUIJN", "");
  assemblyParams.setValue("PATH_MIN_PEAKS", ip.getValue("SPSPATH_MIN_NUM_PEAKS"));
  assemblyParams.setValue("PATH_MIN_SPECS", ip.getValue("SPSPATH_MIN_NUM_SPECS"));
  assemblyParams.setValue("MIN_EDGES_TO_COMPONENT", ip.getValue("SPS_MIN_EDGES_TO_COMPONENT"));

  assemblyParams.setValue("OUTPUT_CLUSTERS", getProjPath(ip, "assembly/path_spectra_as_cluster.txt") );
  assemblyParams.setValue("OUTPUT_SPECS",    getProjPath(ip, "assembly/sps_seqs.pklbin") );
  assemblyParams.setValue("OUTPUT_ABINFO",   getProjPath(ip, "assembly/component_info.bin") );

  assemblyParams.writeToFile(getProjPath(ip, "debug_assembly.params"));

  DEBUG_TRACE;
  ExecAssembly moduleAssembly(assemblyParams,
                               &starSpectra,
                               &starPairs,
                               &contigShifts,
                               &contigAbinfo);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleAssembly.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleAssembly.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecAssembly");

  returnStatus = moduleAssembly.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecAssembly");

  return true;
}

//-----------------------------------------------------------------------------
bool performTagsearch(ParameterList & ip,
	                     SpecSet & spectra,
	                     DB_fasta &db,
	                     vector<unsigned int> *specsToSearch)
{
	ParameterList tagsearchParams;
	tagsearchParams.addIfExists(ip,"RESOLUTION");
	tagsearchParams.addIfExists(ip,"MAX_MOD_MASS");
	tagsearchParams.addIfExists(ip,"TOLERANCE_PEAK");
	tagsearchParams.addIfExists(ip,"TOLERANCE_PM");

  tagsearchParams.addIfExists(ip, "PROJECT_DIR");

	tagsearchParams.addIfExists(ip,"MAX_PARSIMONY");
	tagsearchParams.addIfExists(ip,"TAG_LEN");
	tagsearchParams.addIfDoesntExist("TAG_LEN",                   "6");
	tagsearchParams.addIfExists(ip,"DOUBLE_AA_JUMPS");
	tagsearchParams.addIfDoesntExist("DOUBLE_AA_JUMPS",           "1");
	tagsearchParams.addIfExists(ip,"MATCH_TAG_FLANKING_MASSES");
	tagsearchParams.addIfDoesntExist("MATCH_TAG_FLANKING_MASSES", "0");
	tagsearchParams.addIfExists(ip,"TAG_MATCH_TOP_SCORING_ONLY");

  tagsearchParams.setValue("OUTPUT_PSM",                getProjPath(ip, "assembly/tagsearchpsm.txt") );
	tagsearchParams.setValue("OUTPUT_PEPTIDES",           getProjPath(ip, "spectra/tagsearch.txt") );
	tagsearchParams.writeToFile("debug_tagsearch.params");

  DEBUG_TRACE;
  ExecTagSearch moduleTagsearch(tagsearchParams,
                               &spectra,
                               &db,
                               specsToSearch);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleTagsearch.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleTagsearch.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecTagSearch");

  // Saving output data
  returnStatus = moduleTagsearch.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecTagSearch");

	return true;
}

//-----------------------------------------------------------------------------

bool performSpecProtAlign(ParameterList & ip,
  	                      Clusters & contigs,
  	                      DB_fasta & db,
  	                      SpecSet * matchedSpectra,
  	                      vector<unsigned int> & matchedSpectraIndices)
{
	ParameterList specProtAlignParams;
	specProtAlignParams.addIfExists(ip,"RESOLUTION");
	specProtAlignParams.addIfExists(ip,"MAX_MOD_MASS");
	specProtAlignParams.addIfExists(ip,"TOLERANCE_PEAK");
	specProtAlignParams.addIfExists(ip,"TOLERANCE_PM");

  specProtAlignParams.addIfExists(ip, "PROJECT_DIR");

	if(ip.exists("MIN_NUM_MATCH_PEAKS_DB"))
		specProtAlignParams.setValue("MIN_NUM_MATCH_PEAKS",ip.getValue("MIN_NUM_MATCH_PEAKS_DB"));
	else
		specProtAlignParams.addIfDoesntExist("MIN_NUM_MATCH_PEAKS", "7");
	specProtAlignParams.addIfExists(ip,"MAX_NUM_MODS");
	specProtAlignParams.addIfDoesntExist("MAX_NUM_MODS",          "2");
	specProtAlignParams.addIfExists(ip,"MIN_MOD_MASS");
	specProtAlignParams.addIfDoesntExist("MIN_MOD_MASS",          "-100");
	specProtAlignParams.addIfExists(ip,"MAX_MOD_MASS");
	specProtAlignParams.addIfDoesntExist("MAX_MOD_MASS",          "100");
	specProtAlignParams.addIfExists(ip,"MAX_PARSIMONY");

	specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS",         getProjPath(ip, "spectra/contigs.pklbin") );
	specProtAlignParams.setValue("OUTPUT_MATCHED_SPECS_IDX",     getProjPath(ip, "spectra/contigs_indices.bin") );
	specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX",     getProjPath(ip, "homology/contigs_midx.pklbin") );
	specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS",         getProjPath(ip, "homology/contigs_mp.bin") );
	specProtAlignParams.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", getProjPath(ip, "homology/contigs_midx_all.pklbin") );
	specProtAlignParams.setValue("OUTPUT_MATCHED_PROTS_ALL",     getProjPath(ip, "homology/contigs_mp_all.bin") );
	specProtAlignParams.setValue("OUTPUT_PSM",                   getProjPath(ip, "homology/contigs_psm.txt") );
	specProtAlignParams.setValue("OUTPUT_REF_SPS_NAMES",         getProjPath(ip, "homology/ref_sps_names.txt") );

	specProtAlignParams.setValue("ENFORCE_ENDPEAKS",  "0");
	specProtAlignParams.setValue("SPEC_TYPE_MSMS",    "0");
	specProtAlignParams.setValue("MIN_RATIO",         "0.4");
	specProtAlignParams.writeToFile(getProjPath(ip, "debug_specprotalign.params") );

  DEBUG_TRACE;
  ExecSpecProtAlign moduleSpecProtAlign(specProtAlignParams,
                               &contigs,
                               &db,
                               matchedSpectra,
                               &matchedSpectraIndices);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleSpecProtAlign.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleSpecProtAlign.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecSpecProtAlign");

  returnStatus = moduleSpecProtAlign.saveOutputData();
  DEBUG_VAR(returnStatus);
  TEST_SAVE_OUPUT_DATA("ExecSpecProtAlign");

	return true;
}


//-----------------------------------------------------------------------------

bool performProtProtAlign(ParameterList & ip,
                          unsigned int refProtIdx,
                          set<unsigned int> dbIndexes,
                          DB_fasta & db)
{
  DEBUG_TRACE;
	ParameterList protProtAlignParams;

	protProtAlignParams.addIfExists(ip, "PROJECT_DIR");

	//protProtAlignParams.addIfExists(ip,"CLUSTALW_EXE_DIR");
	if(ip.exists("EXE_DIR"))
	  protProtAlignParams.setValue("CLUSTALW_EXE_DIR", ip.getValue("EXE_DIR"));

	ostringstream tmp;   tmp << refProtIdx;
	protProtAlignParams.setValue("REFERENCE_PROTEIN_IDX",tmp.str());

	protProtAlignParams.addIfExists(ip,"CLUSTALW_MINSCORE");
	protProtAlignParams.addIfDoesntExist("CLUSTALW_MINSCORE", "250");

	protProtAlignParams.setValue("CLUSTALW_FILENAME_PREFIX",  getProjPath(ip, "homology/cwseqs_") );
	protProtAlignParams.setValue("CLUSTALW_INDEX",            getProjPath(ip, "homology/cwindex.txt") );
	protProtAlignParams.writeToFile("debug_protprotalign.params");

  DEBUG_TRACE;
  vector<string> alnFileNames;
  ExecProtProtAlign moduleProtProtAlign(protProtAlignParams,
                               &db,
                               &dbIndexes,
                               &alnFileNames);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleProtProtAlign.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleProtProtAlign.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecProtProtAlign");
  // Save output data
  returnStatus = moduleProtProtAlign.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecProtProtAlign");

	return true;
}

//-----------------------------------------------------------------------------

bool performHomologyAssembly(ParameterList & ip,
                             SpecSet & spectra,
                             DB_fasta & db,
                             SpecSet & contigShifts,
                             SpecSet & contigMatchedIndices,
                             vector<vector<int> > & contigMatchedProts)
{
  DEBUG_TRACE;
	ParameterList homologyAssemblyParams;
	homologyAssemblyParams.setValue("INPUT_HOMOLOGIES",         getProjPath(ip, "homology/cwindex.txt") );
//	homologyAssemblyParams.setValue("INPUT_SPECS_NAMES",      getProjPath(ip, "homology/ref_sps_names.txt") );
	homologyAssemblyParams.setValue("SPEC_TYPE_MSMS",           "0");
	homologyAssemblyParams.setValue("GRAPH_TYPE",               "2");
	homologyAssemblyParams.setValue("MIN_CONTIG_SET",           "1");
	homologyAssemblyParams.setValue("EDGE_SCORE_TYPE",          "1");
	homologyAssemblyParams.setValue("OUTPUT_SPECS",             getProjPath(ip, "homology/homglue_matches.pklbin") );
	homologyAssemblyParams.setValue("OUTPUT_MATCHES_REF_MIDX",  getProjPath(ip, "homology/homglue_ref_midx.pklbin") );
	homologyAssemblyParams.setValue("OUTPUT_MATCHES_REF_MP",    getProjPath(ip, "homology/homglue_ref_mp.bin") );
	homologyAssemblyParams.setValue("OUTPUT_MATCHES_CSPS_MIDX", getProjPath(ip, "homology/homglue_matches_midx.pklbin") );
	homologyAssemblyParams.setValue("OUTPUT_MATCHES_CSPS_MP",   getProjPath(ip, "homology/homglue_matches_mp.bin") );
	homologyAssemblyParams.setValue("OUTPUT_CSV",               getProjPath(ip, "homglue_matches.txt") );
	homologyAssemblyParams.addIfExists(ip,"RESOLUTION");
	homologyAssemblyParams.addIfExists(ip,"MAX_MOD_MASS");
	homologyAssemblyParams.addIfExists(ip,"TOLERANCE_PEAK");
	homologyAssemblyParams.addIfExists(ip,"TOLERANCE_PM");
	homologyAssemblyParams.addIfExists(ip,"MAX_AA_JUMP");
	homologyAssemblyParams.addIfExists(ip,"PENALTY_PTM");
	homologyAssemblyParams.addIfExists(ip,"PENALTY_SAME_VERTEX");

  homologyAssemblyParams.addIfExists(ip, "PROJECT_DIR");

	homologyAssemblyParams.writeToFile(getProjPath(ip, "debug_homology.params") );

  DEBUG_TRACE;
  ExecHomologyAssembly moduleHomologyAssembly(homologyAssemblyParams,
                               &spectra,
                               &db,
                               &contigShifts,
                               &contigMatchedIndices,
                               &contigMatchedProts);
  DEBUG_TRACE;

  string errorString;
  bool isValid = moduleHomologyAssembly.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = moduleHomologyAssembly.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecHomologyAssembly");

  returnStatus = moduleHomologyAssembly.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecHomologyAssembly");

	return true;
}

//-----------------------------------------------------------------------------

bool performReport(ParameterList & ip)
{
  DEBUG_MSG("Entering performReport()");

  bool spsMode = ip.getValueBool("PARTIAL_OVERLAPS",true);

  // Get current directory -- needed for absolute path composition
  string currentWorkindDir;
  currentWorkindDir = getcwd(NULL, 1024);


  // auxiliary string used for file name composition using project name
  string aux;

  // Exec Spsplot section

  ParameterList reportSpsplotParams;

  reportSpsplotParams.addIfExists(ip, "PROJECT_DIR");
  reportSpsplotParams.addIfExists(ip, "EXE_DIR");

  reportSpsplotParams.addIfExists(ip, "REPORT_JOB");
  reportSpsplotParams.addIfExists(ip, "REPORT_USER");
  // add dynamic report project directory (target upon relocation)
  reportSpsplotParams.addIfExists(ip, "REPORT_SERVER");
  reportSpsplotParams.addIfExists(ip, "REPORT_DIR_SERVER");
  reportSpsplotParams.addIfExists(ip, "REPORT_CELLS_PER_LINE");
  // add dynamic reports status
  reportSpsplotParams.addIfExists(ip, "REPORT_DYNAMIC");
  // specify if MS/MS images are shown
  reportSpsplotParams.addIfExists(ip, "REPORT_MSMS_IMAGES");

  // Get report dir -- for the new reports
  string report_dir;
  string report_dir_aux = "report";
  if(ip.exists("REPORT_DIR"))
    report_dir_aux = ip.getValue("REPORT_DIR");

  // If report directory is a relative path, make it absolute
  if(report_dir_aux[0] != '/') {
    report_dir = currentWorkindDir;
    if(report_dir[report_dir.size() - 1] != '/')
      report_dir += '/';
  }
  report_dir += report_dir_aux;

  // set the report directory parameter
  reportSpsplotParams.setValue("OUTDIR", report_dir);


  aux = currentWorkindDir;
  reportSpsplotParams.setValue("FONT_PATH",        aux);

  aux = currentWorkindDir ; aux += "/spectra/stars_only.pklbin";
  reportSpsplotParams.setValue("OUTPUT_STARS",        aux);
  aux = currentWorkindDir ; aux += "/spectra/stars_indices.bin";
  reportSpsplotParams.setValue("OUTPUT_STARS_INDEX",  aux);
  //aux = currentWorkindDir ; aux += "/spectra/stars.pklbin";
  //reportSpsplotParams.setValue("OUTPUT_STARS_ALL",    aux);

  aux = ip.getValue("FASTA_DATABASE");
  string aux2;
  if(aux[0] != '/') {
    aux2 = currentWorkindDir;
    if(aux2[aux2.size() - 1] != '/')
      aux2 += '/';
  }
  aux2 += aux;
  reportSpsplotParams.setValue("FILE_FASTA",       aux);

  if(spsMode) {

    DEBUG_MSG("sps selected");

    aux = currentWorkindDir ; aux += "/homology/ref_sps_names.txt";
    reportSpsplotParams.setValue("FILE_REFINDEX",    aux);
    aux = currentWorkindDir ; aux += "/assembly/component_info.bin";
    reportSpsplotParams.setValue("FILE_COMP",        aux);
    aux = currentWorkindDir ; aux += "/assembly/sps_seqs.pklbin";
    reportSpsplotParams.setValue("FILE_SEQS",        aux);

    // for the old reports
    aux = currentWorkindDir ; aux += "/homology/contigs_mp_all.bin";
    reportSpsplotParams.setValue("FILE_MP",          aux);
    aux = currentWorkindDir ; aux += "/homology/contigs_midx_all.pklbin";
    reportSpsplotParams.setValue("FILE_MIDX",        aux);

    // for the new reports
    aux = currentWorkindDir ; aux += "/homology/contigs_mp.bin";
    reportSpsplotParams.setValue("FILE_MP2",          aux);
    aux = currentWorkindDir ; aux += "/homology/contigs_midx.pklbin";
    reportSpsplotParams.setValue("FILE_MIDX2",        aux);

    aux = currentWorkindDir ; aux += "/homology/homglue_ref_mp.bin";
    reportSpsplotParams.setValue("FILE_REFMP",       aux);
    aux = currentWorkindDir ; aux += "/homology/homglue_ref_midx.pklbin";
    reportSpsplotParams.setValue("FILE_REFMIDX",     aux);

    aux = currentWorkindDir ; aux += "/spectra/contigs.pklbin";
    reportSpsplotParams.setValue("MATCHED_CONTIGS",     aux);
    aux = currentWorkindDir ; aux += "/spectra/contigs_indices.bin";
    reportSpsplotParams.setValue("MATCHED_CONTIGS_IDX",     aux);

  } else {

    DEBUG_MSG("snets selected");

    aux = currentWorkindDir ; aux += "/homology/ref_snets_names.txt";
    reportSpsplotParams.setValue("FILE_REFINDEX",    aux);
    aux = currentWorkindDir ; aux += "/homology/specnets_abruijn.bin";
    reportSpsplotParams.setValue("FILE_COMP",        aux);
//    aux = currentWorkindDir ; aux += "/specnets/snets_specs.pklbin";
    aux = currentWorkindDir ; aux += "/homology/homglue_matches.pklbin";
    reportSpsplotParams.setValue("FILE_SEQS",        aux);

    // for the old reports

//    aux = currentWorkindDir ; aux += "/specnets/snets_mp.bin";
    aux = currentWorkindDir ; aux += "/homology/homglue_matches_mp.bin";
    reportSpsplotParams.setValue("FILE_MP",          aux);
//    aux = currentWorkindDir ; aux += "/specnets/snets_midx.pklbin";
    aux = currentWorkindDir ; aux += "/homology/homglue_matches_midx.pklbin";
    reportSpsplotParams.setValue("FILE_MIDX",        aux);

    // for the new reports

//    aux = currentWorkindDir ; aux += "/specnets/snets_mp.bin";
    aux = currentWorkindDir ; aux += "/homology/homglue_matches_mp.bin";
    reportSpsplotParams.setValue("FILE_MP2",          aux);
//    aux = currentWorkindDir ; aux += "/specnets/snets_midx.pklbin";
    aux = currentWorkindDir ; aux += "/homology/homglue_matches_midx.pklbin";
    reportSpsplotParams.setValue("FILE_MIDX2",        aux);

    // Same as FILE_MP / FILE_MIDX - no homology mapping for spectral networks
    aux = currentWorkindDir ; aux += "/homology/homglue_matches_mp.bin";
    reportSpsplotParams.setValue("FILE_REFMP",       aux);
    aux = currentWorkindDir ; aux += "/homology/homglue_matches_midx.pklbin";
    reportSpsplotParams.setValue("FILE_REFMIDX",     aux);

    aux = currentWorkindDir ; aux += "/homology/homglue_matches.pklbin";
    reportSpsplotParams.setValue("MATCHED_CONTIGS",     aux);
//    aux = currentWorkindDir ; aux += "/spectra/contigs_indices.bin";
//    reportSpsplotParams.setValue("MATCHED_CONTIGS_IDX",     aux);
  }

  aux = currentWorkindDir ; aux += "/spectra/stars.pklbin";
  reportSpsplotParams.setValue("FILE_STARS",       aux);
  aux = currentWorkindDir ; aux += "/spectra/specs_ms.pklbin";
  reportSpsplotParams.setValue("FILE_MS",          aux);
  aux = currentWorkindDir ; aux += "/spectra/input_index.txt";
  reportSpsplotParams.setValue("FILE_INDEX",       aux);

//  aux = currentWorkindDir ; aux += "/spectra/cluster_files.txt";
  //aux = currentWorkindDir ; aux += "/spectra/out/clust/clusters_0_*.clust";
  aux = currentWorkindDir ; aux += "/spectra/out/clust/clusters_0_";
  reportSpsplotParams.setValue("FILE_CLUSTER",     aux);
  aux = currentWorkindDir ; aux += "/pklbin_files.txt";
  //reportSpsplotParams.setValue("FILE_CLUSTERMS",   "./pklbin_files.txt");
  reportSpsplotParams.setValue("FILE_CLUSTERMS",   aux);
  aux = currentWorkindDir ; aux += "/bin_files.txt";
  //reportSpsplotParams.setValue("FILE_CLUSTERSCAN", "./bin_files.txt");
  reportSpsplotParams.setValue("FILE_CLUSTERSCAN", aux);


  // write params to debug file
  reportSpsplotParams.writeToFile(getProjPath(ip, "debug_report.params") );

  // instatiate reports module
  DEBUG_TRACE;
  ExecReportSpsplot moduleReportSpsplot(reportSpsplotParams);
  DEBUG_TRACE;

  // test parameters
  string errorString;
  bool isValid = moduleReportSpsplot.validateParams(errorString);
  TEST_VALID;

  DEBUG_MSG("Invoking reports module");

  bool returnStatus = moduleReportSpsplot.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecReportSpsplot");

  // Copy CSS and Javascript files to report directory
  // skipped if reports are dynamic
  //if(!ip.exists("REPORT_DYNAMIC")) {
    string cssCopyCmd;
    cssCopyCmd = "cp -r ";
    cssCopyCmd += ip.getValue("EXE_DIR");
    cssCopyCmd += "/css/* ";
    cssCopyCmd += ip.getValue("REPORT_DIR");
    int dummy = system(cssCopyCmd.c_str());
  //}


  // ExecReportProteinCoverage section

  DEBUG_MSG("Entering Protein Coverage section");

  ParameterList reportProteinCoverageParams;

  // HTML_DEFS specifies CSS and Javascript root directory path
  reportProteinCoverageParams.addIfExists(ip, "EXE_DIR");
  reportProteinCoverageParams.addIfExists(ip, "HTML_DEFS");
  reportProteinCoverageParams.addIfExists(ip, "TOLERANCE_PEAK");
  reportProteinCoverageParams.addIfExists(ip, "REPORT_TITLE");
  reportProteinCoverageParams.setValue("AA_PER_LINE",              "20");
  //reportProteinCoverageParams.setValue("INPUT_FASTA",              getProjPath(ip, ip.getValue("FASTA_DATABASE")) );
  reportProteinCoverageParams.setValue("INPUT_FASTA",              ip.getValue("FASTA_DATABASE") );
  reportProteinCoverageParams.setValue("CSPS_CONTIG_SPECTRA",      getProjPath(ip, "./homology/homglue_matches.pklbin") );
  reportProteinCoverageParams.setValue("CSPS_CONTIG_MATCHES",      getProjPath(ip, "./homology/homglue_matches_mp.bin") );
  reportProteinCoverageParams.setValue("CSPS_CONTIG_MATCHES_IDX",  getProjPath(ip, "./homology/homglue_matches_midx.pklbin") );

  if(spsMode) {
    reportProteinCoverageParams.setValue("SPS_CONTIG_NAMES",         getProjPath(ip, "./homology/ref_sps_names.txt") );
    reportProteinCoverageParams.setValue("SPS_CONTIG_SPECTRA",       getProjPath(ip, "./spectra/contigs.pklbin") );
    reportProteinCoverageParams.setValue("SPS_CONTIG_MATCHES",       getProjPath(ip, "./homology/homglue_ref_mp.bin") );
    reportProteinCoverageParams.setValue("SPS_CONTIG_MATCHES_IDX",   getProjPath(ip, "./homology/homglue_ref_midx.pklbin") );
  } else {
    reportProteinCoverageParams.setValue("SPS_CONTIG_NAMES",         getProjPath(ip, "./homology/ref_snets_names.txt") );
    reportProteinCoverageParams.setValue("SPS_CONTIG_SPECTRA",       getProjPath(ip, "./homology/homglue_matches.pklbin") );
    reportProteinCoverageParams.setValue("SPS_CONTIG_MATCHES",       getProjPath(ip, "./homology/homglue_matches_mp.bin") );
    reportProteinCoverageParams.setValue("SPS_CONTIG_MATCHES_IDX",   getProjPath(ip, "./homology/homglue_matches_midx.pklbin") );
//    reportProteinCoverageParams.setValue("SPS_CONTIG_SPECTRA",       getProjPath(ip, "./specnets/snets_specs.pklbin") );
//    reportProteinCoverageParams.setValue("SPS_CONTIG_MATCHES",       getProjPath(ip, "./homology/homglue_ref_mp.bin") );
//    reportProteinCoverageParams.setValue("SPS_CONTIG_MATCHES_IDX",   getProjPath(ip, "./homology/homglue_ref_midx.pklbin") );
//    reportProteinCoverageParams.setValue("CSPS_CONTIG_MATCHES",      getProjPath(ip, "./specnets/snets_mp.bin") );
//    reportProteinCoverageParams.setValue("CSPS_CONTIG_MATCHES_IDX",  getProjPath(ip, "./specnets/snets_midx.pklbin") );
  }
//  reportProteinCoverageParams.setValue("OUTDIR",                   getProjPath(ip, "./report") );
  if(ip.exists("REPORT_DIR"))
  	reportProteinCoverageParams.setValue("OUTDIR",ip.getValue("REPORT_DIR"));
  else reportProteinCoverageParams.setValue("OUTDIR",     "report");


  // This section generates protein coverage pages for the static HTML reports (new).
  // If dynamic reports are specifid, coverage pages are not generated
  if(!ip.exists("REPORT_DYNAMIC")) {

    DEBUG_TRACE;
    ExecReportProteinCoverage moduleReportProteinCoverage(reportProteinCoverageParams);
    DEBUG_TRACE;

    errorString;
    isValid = moduleReportProteinCoverage.validateParams(errorString);
    TEST_VALID;

    returnStatus = moduleReportProteinCoverage.loadInputData();
    if (!returnStatus)
    {
      ERROR_MSG("loadInputData ExecReportProteinCoverage");
      return false;
    }

    DEBUG_MSG("Invoking Protein Coverage section");
    // Calling invoke()
    returnStatus = moduleReportProteinCoverage.invoke();
    // Test for return status
    TEST_RETURN_STATUS("ExecReportProteinCoverage");
    // Save output data
    returnStatus = moduleReportProteinCoverage.saveOutputData();
    TEST_SAVE_OUPUT_DATA("ExecReportProteinCoverage");
  }


  return true;
}

//-----------------------------------------------------------------------------
bool performExecMainSpecnets(ParameterList & ip,
    SpecSet * msSpectra,
    SpecSet * scoredSpectra,
    SpecSet * starSpectra,
    SpectrumPairSet * pairs,
    DB_fasta * db,
    PeptideSpectrumMatchSet * psms,
    SpecSet * psms_spectra,
    SpecSet * psms_midx,
    vector<vector<int> > * psms_mp,
    SpecSet * snets_contigs,
    SpecSet * snets_midx,
    vector<vector<int> > * snets_mp)
{
  ParameterList specnetsParams;
  specnetsParams.setValue("RESOLUTION",     ip.getValue( "RESOLUTION",     "0.1") );
  specnetsParams.setValue("MAX_MOD_MASS",   ip.getValue( "MAX_MOD_MASS",   "100") );
  specnetsParams.setValue("TOLERANCE_PEAK", ip.getValue( "TOLERANCE_PEAK", "0.45") );
  specnetsParams.setValue("TOLERANCE_PM",   ip.getValue( "TOLERANCE_PM",   "1.5") );
  specnetsParams.setValue("SPECNETS_PROJ_TYPE",   ip.getValue( "SPECNETS_PROJ_TYPE","all") ); // Possible values are "all" / "matched" to retain all/matched-only PRMs during propagation

  // Tagsearch params
  specnetsParams.addIfExists(ip,"INSPECT_PSMS");
  specnetsParams.setValue("OUTPUT_PEPTIDES","spectra/tagsearch_peptides.txt");

  // SVM params
  specnetsParams.setValue("STARS_SCORED_PRM_OFFSET", "-1.0072763");
  specnetsParams.setValue("STARS_SCORED_SRM_OFFSET", "-1.0072763");
  specnetsParams.setValue("SCAN_ZERO_INDEX",         "1");
  specnetsParams.addIfExists(ip,"MS2_SCORING_MODEL");               // To become "<trunk>/resources/dancik_model.txt
  specnetsParams.addIfExists(ip,"SPECS_MS_STATISTICS_CONFIG");      // To become "<trunk>/resources/spectra_stats_ms.txt
  specnetsParams.addIfExists(ip,"SPECS_SCORED_STATISTICS_CONFIG");  // To become "<trunk>/resources/spectra_stats_prm.txt
  specnetsParams.addIfExists(ip,"STARS_STATISTICS_CONFIG");         // To become "<trunk>/resources/spectra_stats_stars.txt
  specnetsParams.addIfExists(ip,"SVM_SCALE_PARAMS_CHARGE1");        // To become "<trunk>/resources/HEK293_charge1_model_svm_keys_range.txt
  specnetsParams.addIfExists(ip,"SVM_SCALE_PARAMS_CHARGE2");        // To become "<trunk>/resources/HEK293_charge2_model_svm_keys_range.txt
  specnetsParams.addIfExists(ip,"SVM_SCALE_PARAMS_CHARGE3");        // To become "<trunk>/resources/HEK293_charge3_model_svm_keys_range.txt
  specnetsParams.addIfExists(ip,"SVM_MODEL_CHARGE1");               // To become "<trunk>/resources/HEK293_charge1_model_SVM.model
  specnetsParams.addIfExists(ip,"SVM_MODEL_CHARGE2");               // To become "<trunk>/resources/HEK293_charge2_model_SVM.model
  specnetsParams.addIfExists(ip,"SVM_MODEL_CHARGE3");               // To become "<trunk>/resources/HEK293_charge3_model_SVM.model

  // FDR params
  specnetsParams.setValue("INPUT_RESULTS_TYPE", "snets");
  specnetsParams.setValue("OUTPUT_FDR_RESULTS", "spectra/psms_fdr.txt");

  specnetsParams.writeToFile("debug_specnets.params");

  MS2ScoringModel model;
  if (ip.exists("MS2_SCORING_MODEL"))
  {
    if (!model.LoadModel(ip.getValue("MS2_SCORING_MODEL").c_str()))
    {
      DEBUG_MSG("Could not load " << ip.getValue("MS2_SCORING_MODEL"));
      return false;
    }
  }
  DEBUG_MSG("SCORING MODEL LOADED");

  DEBUG_TRACE;
  ExecMainSpecnets mainSpecnets(specnetsParams,
                                msSpectra,
                                scoredSpectra,
                                starSpectra,
                                pairs,
                                & model,
                                db,
                                psms,
                                psms_spectra,
                                psms_midx,
                                psms_mp,
                                snets_contigs,
                                snets_midx,
                                snets_mp);
  DEBUG_TRACE;

  string errorString;
  bool isValid = mainSpecnets.validateParams(errorString);
  TEST_VALID;

  bool returnStatus = mainSpecnets.invoke();
  // Test for return status
  TEST_RETURN_STATUS("ExecMainSpecnets");

  // Saving output data
  returnStatus = mainSpecnets.saveOutputData();
  TEST_SAVE_OUPUT_DATA("ExecMainSpecnets");

  ofstream spsIndex(getProjPath(ip, "homology/ref_snets_names.txt").c_str(),ios_base::out);
  for(unsigned int i=0; i<snets_contigs->size(); i++)
    spsIndex<<"snet:"<<i+1<<endl;
  spsIndex.close();

  //---------------------------------------------------------------------------
  // REPORT STAGE
  //---------------------------------------------------------------------------

  ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(), ios_base::out);
  spsProj << "sps;.;"<<ip.getValue("TOLERANCE_PEAK")<<";"<<ip.getValue("TOLERANCE_PM")<<"\n";
  spsProj.close();

  if (!performReport(ip))
  {
    ERROR_MSG("Problem encountered during Report stage" );
    exit(-100-STAGE_REPORT);
  }

  return true;
}

bool writeStatusFile(string fileName, string status)
{

  FILE * f = fopen(fileName.c_str(),"w");

  if(f == NULL)
    {
      perror("The following error occurred");
    ERROR_MSG("Problem encountered creating status file");
      exit(-1);
    }

  int val = fwrite(status.c_str(),sizeof(char),strlen(status.c_str()),f);
  if(val != strlen(status.c_str()))
    {
      ERROR_MSG("Problem encountered writing status file");
ERROR_MSG(val);
ERROR_MSG(strlen(status.c_str()));
      exit(-1);
    }
  fclose(f);
  DEBUG_MSG("Updated program status to " + status);
}

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// MAIN
//-----------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  int logLevel = 5;
  Logger::setDefaultLogger(Logger::getLogger(logLevel));

  string runMode("");
  string initialStageString("");
  string statusFileName("status.txt");
  bool resumeFlag = false;
  bool gridExecutionFlag = false;
  bool runGenoMSFlag = false;
  bool runMergedFlag = false;

  bool showHelp = false;
  for (size_t i = 0; i < argc; i++) {
    string arg(argv[i]);
    if (arg.compare("--help") == 0) {
      showHelp = true;
    }
  }

  if (argc < 2 || showHelp)
  {
    if (showHelp) {
      cout << "Usage: main_specnets [PARAM FILE] [OPTION]..." << endl << endl;
      cout << "Optional arguments are listed below " << endl;
      cout << "  -i  <intialstage>   begin processing at specified stage:" << endl;
      cout << "                         mscluster, scoring, filterpairs" << endl;
      cout << "                         filteraligns, alignment, filterstarpairs" << endl;
      cout << "                         assembly, tagsearch, specnets, specprotalign" << endl;
      cout << "                         protprotalign, homologyassembly, genoms, report" << endl;
      cout << "  -g                  execution is on a grid" << endl;
      cout << "  -lf <filename>      name of log file for output" << endl;
      cout << "  -ll <loglevel>      log level for debug/warn/error output:" << endl;
      cout << "                         9 for errors only" << endl;
      cout << "                         5 for warnings and errors" << endl;
      cout << "                         0 for all debug output" << endl;
      cout << "  -r  <runmode>        run mode (sps or specnets):" << endl;
      cout << "       NOTE: THIS OPTION IS NOT IMPLEMENTED AT THIS TIME" << endl;
      cout << "  -s                   execute a single step then exit" << endl;
      cout << "  -z                   resume an earlier run" << endl;
      cout << "       NOTE: this option is only used when intermediate results" << endl;
      cout << "             were saved during a previous run where FilterPairs" << endl;
      cout << "             was executed on a remote grid" << endl;
      cout << "  -q                   Run in GenoMS mode" << endl;
      cout << "       NOTE: In this mode, only 4 stages are valid: mscluster, " << endl;
      cout << "             scoring, genoms, and report" << endl ;
      cout << "  -m                   Run in integrated mode (GenoMS + CSPS)" << endl << endl;
    }
    else {
      cerr << "main_specnets: insufficient arguments" << endl;
      cerr << "Try \'main_specnets --help\' for more information." << endl << endl;
    }



    cout << PROGRAM_NAME << endl;
    cout << "main_specnets 3.0." << XSTR(SPS_VERSION) << endl;
    cout << endl;

    cout << COPYRIGHT1 << endl;
    cout << COPYRIGHT2 << endl;
    cout << endl;

    return -1;
  }

  // Parse the command line parameters
  vector<CommandLineParser::Option> listOptions;
  listOptions.push_back(CommandLineParser::Option("i", "INITIAL_STAGE", 1));
  listOptions.push_back(CommandLineParser::Option("g", "GRID_EXECUTION", 0));
  listOptions.push_back(CommandLineParser::Option("lf", "LOG_FILE_NAME", 1));
  listOptions.push_back(CommandLineParser::Option("ll", "LOG_LEVEL", 1));
  listOptions.push_back(CommandLineParser::Option("r", "RUN_MODE", 1));
  listOptions.push_back(CommandLineParser::Option("s", "SINGLE_STEP", 0));
  listOptions.push_back(CommandLineParser::Option("z", "RESUME_FLAG", 0));
  listOptions.push_back(CommandLineParser::Option("q", "GENOMS_FLAG", 0));
  listOptions.push_back(CommandLineParser::Option("m", "MERGE_FLAG", 0));

  CommandLineParser clp(argc, argv, 1, listOptions);
  string parserError = "";
  if (!clp.validate(parserError))
  {
    cerr << "main_specnets: " << parserError << endl;
    cerr << "Try \'main_specnets --help\' for more information." << endl << endl;
    cout << PROGRAM_NAME << endl;
    cout << "main_specnets 3.0." << XSTR(SPS_VERSION) << endl;
    cout << endl;
    cout << COPYRIGHT1 << endl;
    cout << COPYRIGHT2 << endl;
    cout << endl;
    return -1;
  }



  ParameterList commandLineParams;
  clp.getOptionsAsParameterList(commandLineParams);
  if(commandLineParams.exists("GENOMS_FLAG"))
    commandLineParams.setValue("GENOMS_FLAG","1");

  if(commandLineParams.exists("MERGE_FLAG"))
    commandLineParams.setValue("MERGE_FLAG","1");


  ParameterList ip;
  ip.readFromFile(argv[1]);
  ip.writeToFile("debug_sps.params");

  // Combine the command line parameters to the file ones
  //   Command line parameters take precedence (hence the overwrite flag set)
  ip.addList(commandLineParams, true);
  ip.writeToFile("debug_wcommand.params");

  if(ip.getValueInt("GENOMS_FLAG",0) == 1)
    runGenoMSFlag = true;

  if(ip.getValueInt("MERGE_FLAG",0) == 1)
    runMergedFlag = true;

  logLevel = ip.getValueInt("LOG_LEVEL", 5);
  if (ip.exists("LOG_FILE_NAME"))
  {
    string logFileName = ip.getValue("LOG_FILE_NAME");
    Logger::setDefaultLogger(Logger::getLogger(logFileName, logLevel));
  } else {
    Logger::setDefaultLogger(Logger::getLogger(logLevel));
  }
  DEBUG_MSG("GenoMS_FLAG=" + ip.getValue("GENOMS_FLAG","X"));
  DEBUG_MSG("Merge_FLAG=" + ip.getValue("MERGE_FLAG","X"));

  DEBUG_TRACE;

  if (ip.exists("RUN_MODE"))
  {
    runMode = commandLineParams.getValue("RUN_MODE");
  }
  if (ip.exists("INITIAL_STAGE"))
  {
    initialStageString = commandLineParams.getValue("INITIAL_STAGE");
  }

  DEBUG_VAR(initialStageString);

  if (ip.exists("MIN_PVALUE")) {
	ERROR_MSG("MIN_PVALUE is a deprecated variable, must use MAX_PVALUE instead");
	return -1;
  }

  addDefaultParameterValues(ip);
  ip.writeToFile("debug_default.params");

  if (ip.exists("EXE_DIR"))
  {
    // In case there is a "/" at the end of EXE_DIR.. remove it
    string exeDir = ip.getValue("EXE_DIR");
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
    if (gridSgeExeDir.length() > 2 && gridSgeExeDir[gridSgeExeDir.length() - 1] == '/')
    {
      gridSgeExeDir = gridSgeExeDir.substr(0, gridSgeExeDir.length() - 1);
      ip.setValue("GRID_SGE_EXE_DIR", gridSgeExeDir);
    }
  }


  if (ip.exists("INITIAL_STAGE"))
  {
    initialStageString = ip.getValue("INITIAL_STAGE");
  }

  DEBUG_VAR(initialStageString);

  if (initialStageString.empty())
  {
    initialStageString = "begin";
  }

  DEBUG_VAR(initialStageString);

  map<string, Stage> map_stage;
  map_stage["begin"]            = STAGE_BEGIN;
  map_stage["mscluster"]        = STAGE_BEGIN;
  map_stage["scoring"]          = STAGE_SCORING;
  map_stage["filterpairs"]      = STAGE_FILTERPAIRS;
  map_stage["filteraligns"]     = STAGE_FILTERALIGNS;
  map_stage["alignment"]        = STAGE_ALIGNMENT;
  map_stage["filterstarpairs"]  = STAGE_FILTERSTARPAIRS;
  map_stage["assembly"]         = STAGE_ASSEMBLY;
  map_stage["tagsearch"]        = STAGE_TAGSEARCH;
  map_stage["specnets"]         = STAGE_SPECNETS;
  map_stage["specprotalign"]    = STAGE_SPECPROTALIGN;
  map_stage["protprotalign"]    = STAGE_PROTPROTALIGN;
  map_stage["homologyassembly"] = STAGE_HOMOLOGYASSEMBLY;
  map_stage["genoms"]           = STAGE_GENOMS;
  map_stage["merge"]            = STAGE_MERGE;
  map_stage["report"]           = STAGE_REPORT;

  if (map_stage.find(initialStageString) == map_stage.end())
  {
    ERROR_MSG("Unknown starting stage [" << initialStageString << "]");
    return -1;
  }

  //Start the status as "running" and write itout
  writeStatusFile(statusFileName,"Running");



  int initialStage = map_stage[initialStageString];
  DEBUG_VAR(initialStage);

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

  //
  // Initialize environment
  //

  DEBUG_MSG("Spectral Networks 2.0.0: session started on: " << getCurrentTimeString());
  DEBUG_MSG("Starting stage: [" << initialStageString << "] on: " << getCurrentTimeString());

  bool res;
  res = mkdir_if_not_exist("spectra");
  if (res) {
    DEBUG_MSG("Made directory \'spectra\'");
  }
  res = mkdir_if_not_exist("aligns");
  if (res) {
    DEBUG_MSG("Made directory \'aligns\'");
  }
  res = mkdir_if_not_exist("specnets");
  if (res) {
    DEBUG_MSG("Made directory \'specnets\'");
  }
  res = mkdir_if_not_exist("homology");
  if (res) {
    DEBUG_MSG("Made directory \'homology\'");
  }
  res = mkdir_if_not_exist("report");
  if (res) {
    DEBUG_MSG("Made directory \'report\'");
  }
  res = mkdir_if_not_exist("assembly");
  if (res) {
    DEBUG_MSG("Made directory \'assembly\'");
  }


  // get LD_LIBRARY_PATH from system
  char *curLibPath = getenv("LD_LIBRARY_PATH");


  // Build the needed library path
  string libPath;
  libPath = ip.getValue("EXE_DIR");
  // set LD_LIBRARY_PATH to EXE_DIR + /libs.
  libPath += "/libs";

  string fullLibPath;
  // check if LD_LIBRARY_PATH is already defined.
  if(curLibPath) {
    // if it is, check if it contains the path we want.
    fullLibPath = curLibPath;
    // if the library path IS NOT contained in the path variable, add it, and set the environment variable.
    if(fullLibPath.find(libPath) == string::npos) {
      fullLibPath += ':';
      fullLibPath += libPath;
      setenv("LD_LIBRARY_PATH", fullLibPath.c_str(), 1);
    }
  } else {
    // if LD_LIBRARY_PATH is not defined,, define it with what we want.
    setenv("LD_LIBRARY_PATH", libPath.c_str(), 1);
  }


  //---------------------------------
  // Load amino acid masses
  //---------------------------------
  AAJumps jumps(-1);	// Amino acid masses
  if (ip.exists("AMINO_ACID_MASSES"))
  {
    jumps.loadJumps(ip.getValue("AMINO_ACID_MASSES").c_str(), true);
  }

  DEBUG_TRACE;

  SpecSet ms2Spectra;  // MS/MS spectra to be used when creating PRM spectra
  string exeDir = ip.getValue("EXE_DIR");
  string convertCmd = exeDir + "/convert ";

  if (initialStage == STAGE_BEGIN)
  {
    DEBUG_TRACE;
    vector<string> spectrumNames; // Name of each spectrum for final report

    // Build the needed library path
    string libPath = exeDir + "/libs";
	  addEnvironmentVariable(convertCmd, "LD_LIBRARY_PATH", libPath);

    //---------------------------------
    // Parse/Load input MS/MS spectra
    //---------------------------------
    if(ip.getValueInt("CLUSTER_MIN_SIZE") > 0)
    {
      DEBUG_TRACE;
    	loadInitialData(ip,
    	                convertCmd,
    	                ip.getValue("INPUT_SPECS_MS"),
    	                getProjPath(ip, "spectra/specs_ms_") );
    }
    else
    {
      DEBUG_TRACE;
    	loadInitialData(ip,
    	                convertCmd,
    	                ip.getValue("INPUT_SPECS_MS"),
    	                getProjPath(ip, "spectra/specs_ms_" ),
    	                &ms2Spectra);
    }

    if(ip.getValueInt("CLUSTER_MIN_SIZE") > 0)
    {
      DEBUG_MSG("Starting stage McCluster on: " << getCurrentTimeString());
      if (!performMsCluster(ip, ms2Spectra, spectrumNames))
      {
        ERROR_MSG("Problem encountered during MsCluster stage" );
	writeStatusFile(statusFileName,"Error");
	exit(-1);
      }
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing MsCluster stage" );
    if (!ms2Spectra.LoadSpecSet_pklbin(getProjPath(ip,
                                                   "./spectra/specs_ms.pklbin").c_str())) {
      DEBUG_MSG(getProjPath(ip,"./spectra/specs_ms.pklbin").c_str() << " not found, attempting to load from INPUT_SPECS_MS instead");
      // Still try to load MS/MS spectra if bypassing MSCluster stage

      if (!loadInitialData(ip,
                           convertCmd,
                           ip.getValue("INPUT_SPECS_MS"),
                           getProjPath(ip, "spectra/specs_ms_"),
                           &ms2Spectra)) {

        ERROR_MSG("Could not find input MS/MS spectra!");
      }

    }

    DEBUG_VAR(ms2Spectra.size());
  }

  // Test for empty return data structures
  TEST_RETURN("moduleMsCluster", ms2Spectra);

  DEBUG_VAR(ms2Spectra.size());

  if (ms2Spectra.size() == 0)
  {
    ERROR_MSG("ms2Spectra size is 0!");
    writeStatusFile(statusFileName,"Error");
    exit(-2);
  }

  //---------------------------------------------------------------------------
  // PEPNOVO/SCORING STAGE
  //---------------------------------------------------------------------------
  SpecSet prmSpectra;
  if (initialStage <= STAGE_SCORING)
  {
    DEBUG_MSG("Starting stage Scoring on: " << getCurrentTimeString());
    if (!performScoring(ip, ms2Spectra, prmSpectra))
    {
      ERROR_MSG("Problem encountered during Scoring stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-3);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      prmSpectra.SaveSpecSet_pklbin(getProjPath(ip, "./spectra/specs_scored.pklbin").c_str() );
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing Scoring stage" );
    if (!prmSpectra.LoadSpecSet_pklbin(getProjPath(ip,
                                                   "./spectra/specs_scored.pklbin").c_str())) {
      DEBUG_MSG(getProjPath(ip,"./spectra/specs_scored.pklbin").c_str() << " not found, attempting to load from INPUT_SPECS_PRM instead");

      if ((!ip.exists("INPUT_SPECS_PRM"))
          || (!loadInitialData(ip,
                                convertCmd,
                                ip.getValue("INPUT_SPECS_PRM"),
                                getProjPath(ip, "spectra/specs_scored_"),
                                &prmSpectra))) {
        ERROR_MSG("Could not find input PRM spectra!");
      }
    }
    DEBUG_VAR(prmSpectra.size());
  }

  // Test for empty return data structures
  TEST_RETURN("moduleScoring", prmSpectra);

  DEBUG_TRACE;
  float peakTol = ip.getValueDouble("TOLERANCE_PEAK", 0.5);
  bool specTypeMSMS = ip.getValueBool("SPEC_TYPE_MSMS", false);
  float ionOffset = specTypeMSMS ? AAJumps::massHion : 0;
  for (unsigned int i = 0; i < prmSpectra.size(); i++)
  {
    prmSpectra[i].addZPMpeaks(peakTol, ionOffset, true);
  }


  //--------------------------------------------------------------------------
  // GenoMS STAGE
  //--------------------------------------------------------------------------

  if (runGenoMSFlag || runMergedFlag)
  {
    if(initialStage <= STAGE_GENOMS)
    {
      DEBUG_MSG("Starting GenoMS");

      if(!performGenoMS(ip)) {
	ERROR_MSG("Problem encountered during GenoMS");
	writeStatusFile(statusFileName,"Error");
	exit(-1);
      }
    } else {
      DEBUG_MSG("Bypassing GenoMS stage");
    }
    //SpecSet testSpectra;

    //testSpectra.LoadSpecSet_pklbin( getProjPath(ip, "./spectra/stars.pklbin").c_str() );
    //DEBUG_VAR(testSpectra.size());
    //testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/stars.mgf").c_str());
    //DEBUG_TRACE;

    //testSpectra.LoadSpecSet_pklbin( getProjPath(ip, "./spectra/specs_ms.pklbin").c_str() );
    //DEBUG_VAR(testSpectra.size());
    //testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/spec_ms2.mgf").c_str());
    //DEBUG_TRACE;

    //testSpectra.LoadSpecSet_pklbin( getProjPath(ip, "./assembly/sps_seqs.pklbin").c_str() );
    //DEBUG_VAR(testSpectra.size());
    //testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./assembly/sps_seqs.mgf").c_str());
    //DEBUG_TRACE;

    //testSpectra.LoadSpecSet_pklbin( getProjPath(ip, "./homology/contigs_midx_all.pklbin").c_str() );
    //DEBUG_VAR(testSpectra.size());
    //testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./homology/contigs_midx_all.mgf").c_str());
    //DEBUG_TRACE;

    //testSpectra.LoadSpecSet_pklbin( getProjPath(ip, "./spectra/contigs.pklbin").c_str() );
    //DEBUG_VAR(testSpectra.size());
    //testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./spectra/contigs.mgf").c_str());
    //DEBUG_TRACE;

    //testSpectra.LoadSpecSet_pklbin( getProjPath(ip, "./homology/homglue_matches.pklbin").c_str() );
    //DEBUG_VAR(testSpectra.size());
    //testSpectra.SaveSpecSet_mgf(getProjPath(ip, "./homology/homglue_matches.mgf").c_str());
    //DEBUG_TRACE;
  }

  if(runMergedFlag && initialStage <= STAGE_GENOMS) {

    DEBUG_MSG("Stage Genoms: " << STAGE_GENOMS);
    //We need to rename all of our results files
    if(rename("./spectra/stars.pklbin","./spectra/genoms.stars.pklbin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./assembly/sps_seqs.pklbin","./assembly/genoms.sps_seqs.pklbin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./assembly/component_info.bin","./assembly/genoms.component_info.bin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/contigs_mp.bin","./homology/genoms.contigs_mp.bin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/contigs_midx.pklbin","./homology/genoms.contigs_midx.pklbin") < 0) {
    	  ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/contigs_mp_all.bin","./homology/genoms.contigs_mp_all.bin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/contigs_midx_all.pklbin","./homology/genoms.contigs_midx_all.pklbin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./spectra/contigs.pklbin","./spectra/genoms.contigs.pklbin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./spectra/contigs_indices.bin","./spectra/genoms.contigs_indices.bin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/homglue_ref_mp.bin","./homology/genoms.homglue_ref_mp.bin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/homglue_ref_midx.pklbin","./homology/genoms.homglue_ref_midx.pklbin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/homglue_matches.pklbin","./homology/genoms.homglue_matches.pklbin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/homglue_matches_mp.bin","./homology/genoms.homglue_matches_mp.bin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/homglue_matches_midx.pklbin","./homology/genoms.homglue_matches_midx.pklbin") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }
    if(rename("./homology/ref_sps_names.txt","./homology/genoms.ref_sps_names.txt") < 0) {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
	  }

    if(rename("./protid.fasta","./genoms.protid.fasta") < 0)
    {
      ERROR_MSG("Problem encountered renaming GenoMS files");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_GENOMS);
    }

  } else if(runGenoMSFlag) {

    //If we are running GenoMS, then we need to update FASTA_DATABASE to point to our output

    ip.setValue("FASTA_DATABASE","./protid.fasta");


    if (initialStage <= STAGE_REPORT)
    {
      DEBUG_MSG("Starting Report stage");
      ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(), ios_base::out);
      spsProj << "sps;.;"<<ip.getValue("TOLERANCE_PEAK")<<";"<<ip.getValue("TOLERANCE_PM")<<"\n";
      spsProj.close();

      if (!performReport(ip)) {
        ERROR_MSG("Problem encountered during Report stage" );
	writeStatusFile(statusFileName,"Error");
	exit(-STAGE_REPORT);
      }
      //--------------------------------------------------------------------------
      //Set up for relaunching
      //--------------------------------------------------------------------------
      if(!generateRelaunchScript(ip))
	{
	  ERROR_MSG("Problem encountered during relaunch script creation");
	  writeStatusFile(statusFileName,"Error");
	  exit(-STAGE_REPORT);
	}

    } else  {
      DEBUG_MSG("Bypassing Report stage" );
    }
    writeStatusFile(statusFileName,"Finished");
    return 0;
  }

  SpectrumPairSet filteredPairs;

  vector<TwoValues<float> > ratios;
  vector<TwoValues<float> > means;
  vector<float> varTerms;
  list<vector<float> > alignStats;
  vector<vector<float> > specStats;
  vector<unsigned int> idxKept;
  vector<TwoValues<float> > pvalues;

  //---------------------------------------------------------------------------
  // FLITERPAIRS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERPAIRS)
  {
    DEBUG_MSG("Starting stage FilterPairs on: " << getCurrentTimeString());
    if (!performFilterPairs(ip,
                            prmSpectra,
                            ms2Spectra,
                            filteredPairs,
                            ratios,
                            means,
                            varTerms,
                            alignStats,
                            specStats,
                            idxKept,
                            pvalues,
                            gridExecutionFlag,
                            resumeFlag))
    {
      ERROR_MSG("Problem encountered during Filter Pairs stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-4);
    }

    if((!ip.exists("GRID_NUMNODES") || ip.getValueInt("GRID_NUMNODES") > 0) &&
       (!resumeFlag && !gridExecutionFlag))
    {
      // If we are doing a grid execution (and are not actually on the grid)
      //    and we aren't resuming... then exit (we'll resume execution later)
      DEBUG_MSG("Files for grid execution have been saved.");
      DEBUG_MSG("Restart with -z option when grid execution has been completed.");
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      filteredPairs.saveToBinaryFile( getProjPath(ip, "./aligns/pairs_raw.bin") );
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }
    else

	filteredPairs.saveToBinaryFile( getProjPath(ip, "./aligns/pairs_raw.bin") );
  }
  else
  {
    DEBUG_MSG("Bypassing Filter Pairs stage" );
    filteredPairs.loadFromBinaryFile( getProjPath(ip, "./aligns/pairs_raw.bin") );
    DEBUG_VAR(filteredPairs.size());
    Load_binArray("./aligns/ratios.bin", ratios);
    DEBUG_VAR(ratios.size());
    Load_binArray("./aligns/means.bin", means);
    DEBUG_VAR(means.size());
    Load_binArray("./aligns/vars.bin", varTerms);
    DEBUG_VAR(varTerms.size());
  }

  if(ip.getValue("PAIRS_MATCH_MODE","") == "cosine") {
	  DEBUG_MSG("Parameters include PAIRS_MATCH_MODE=cosine, exiting after ExecFilterPairs" );
	  return 0;
  }


  // Test for empty return data structures
//  TEST_RETURN("moduleFilterPairs", ratios);
//  TEST_RETURN("moduleFilterPairs", means);
//  TEST_RETURN("moduleFilterPairs", varTerms);
  TEST_RETURN("moduleFilterPairs", filteredPairs);
//  TEST_RETURN("moduleFilterPairs", pvalues);
//  TEST_RETURN("moduleFilterPairs", idxKept);
//  TEST_RETURN("moduleFilterPairs", alignStats);
//  TEST_RETURN("moduleFilterPairs", specStats);

  //---------------------------------------------------------------------------
  // FLITERALIGNS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERALIGNS)
  {
    DEBUG_MSG("Starting stage FilterAligns on: " << getCurrentTimeString());
    if (!performFilterAligns(ip,
                            filteredPairs,
                            ratios,
                            means,
                            varTerms,
                            idxKept,
                            pvalues))
    {
      ERROR_MSG("Problem encountered during Filter Aligns stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-4);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      filteredPairs.saveToBinaryFile( getProjPath(ip, "./aligns/pairs.bin") );
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }
    else
      filteredPairs.saveToBinaryFile( getProjPath(ip, "./aligns/pairs.bin") );
  }
  else
  {
    DEBUG_MSG("Bypassing Filter Aligns stage" );
    filteredPairs.loadFromBinaryFile( getProjPath(ip, "./aligns/pairs.bin") );
    DEBUG_VAR(filteredPairs.size());
  }

  // Test for empty return data structures
  TEST_RETURN("moduleFilterPairs", filteredPairs);

  SpecSet pairAlignments;
  SpecSet starSpectraOnly;
  SpecSet starSpectra;
  vector<unsigned int> alignedSpectra;

  //---------------------------------------------------------------------------
  // ALIGNMENT STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_ALIGNMENT)
  {
    DEBUG_MSG("Starting stage Alignment on: " << getCurrentTimeString());
    if (!performAlignment(ip,
                          jumps,
                          prmSpectra,
                          filteredPairs,
                          pairAlignments,
                          starSpectraOnly,
                          starSpectra,
                          alignedSpectra))
    {
      ERROR_MSG("Problem encountered during Alignment stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-5);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing Alignment stage" );
    starSpectra.LoadSpecSet_pklbin( getProjPath(ip, "./spectra/stars.pklbin").c_str() );
    DEBUG_VAR(starSpectra.size());
  }

  // Test for empty return data structures
  //TEST_RETURN("ExecAlignment", pairAlignments);
  //TEST_RETURN("ExecAlignment", starSpectraOnly);
  TEST_RETURN("ExecAlignment", starSpectra);
  //TEST_RETURN("ExecAlignment", alignedSpectra);

  vector<vector<float> > starRatios;
  SpecSet matchedPeaks;

  //---------------------------------------------------------------------------
  // FILTERSTARPAIRS STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_FILTERSTARPAIRS)
  {
    if (!performFilterStarPairs(ip,
                               filteredPairs,
                               starSpectra,
                               starRatios,
                               matchedPeaks))
    {
      ERROR_MSG("Problem encountered during Filter Star Pairs stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_FILTERSTARPAIRS);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }
  }
  else
  {
    DEBUG_MSG("Bypassing Filter Star Pairs stage" );
//    DEBUG_MSG("Not implemented yet" );
//    Load_binArray("./aligns/ratios.bin", starRatios);
//    DEBUG_VAR(starRatios.size());
//    matchedPeaks->LoadSpecSet_pklbin("");
//    DEBUG_VAR(matchedPeaks.size());

    if(!filteredPairs.loadFromBinaryFile("aligns/pairs_stars.bin")) {
      ERROR_MSG("Problem encountered during Filter Star Pairs stage (loading filteredPairs)" );
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_FILTERSTARPAIRS);
    }
  }


  DEBUG_MSG("Reloading Star specturm" );
  // Need to reload the star spectra because they were altered by FilterStarPairs
	if(starSpectra.LoadSpecSet_pklbin( getProjPath(ip, "./spectra/stars.pklbin").c_str() ) <=0 ) {
    ERROR_MSG("Problem encountered while reloading star specturm" );
    writeStatusFile(statusFileName,"Error");
    exit(-STAGE_ASSEMBLY);
	}
  DEBUG_VAR(starSpectra.size());


  if(ip.getValue("EXECUTION_MODE","") == string("signatures")) {
  	return 0;
  }


  // Test for empty return data structures
  //TEST_RETURN("ExecFilterStarPairs", starRatios);
  //TEST_RETURN("ExecFilterStarPairs", matchedPeaks);

  //---------------------------------------------------------------------------
  // Spectral Networks STAGE
  //---------------------------------------------------------------------------
  DB_fasta dbAll;
  string dbFileName;
  DEBUG_MSG(ip.getValue("FASTA_DATABASE"));
  if(ip.exists("FASTA_DATABASE")) {

    dbFileName = ip.getValue("FASTA_DATABASE");

    //If we are integrating, then add GenoMS protein sequences to the DB
    if(runMergedFlag && initialStage <= STAGE_FILTERSTARPAIRS)
      {
	DEBUG_MSG("Concatenating " << dbFileName << " and genoms.protid.fasta");
	if(!concatenateFiles("./genoms.protid.fasta",dbFileName,"./protid.fasta"))
	  {
	    ERROR_MSG("Problem encountered concatenating fasta files");
	    writeStatusFile(statusFileName,"Error");
	    exit(-100);
	  }
	dbFileName = "./protid.fasta";
	ip.setValue("FASTA_DATABASE","./protid.fasta");
      }

    if(dbAll.Load(dbFileName.c_str())<=0) {
      ERROR_MSG("Problem encountered during Spectral Networks stage (loading "<<dbFileName<<")" );
      writeStatusFile(statusFileName,"Error");
      exit(-100);
    }
  }

  if(not ip.getValueBool("PARTIAL_OVERLAPS",true)) {
    // Entering Spectral Networks mode
    PeptideSpectrumMatchSet * psms = new PeptideSpectrumMatchSet;
    SpecSet * psms_spectra = new SpecSet;
    SpecSet * psms_midx = new SpecSet;
    vector<vector<int> > * psms_mp = new vector<vector<int> >;
    SpecSet * snets_contigs = new SpecSet;
    SpecSet * snets_midx = new SpecSet;
    vector<vector<int> > * snets_mp = new vector<vector<int> >;

    if(not performExecMainSpecnets(ip,
                                   &ms2Spectra,
                                   &prmSpectra,
                                   &starSpectra,
                                   &filteredPairs,
                                   &dbAll,
                                   psms,
                                   psms_spectra,
                                   psms_midx,
                                   psms_mp,
                                   snets_contigs,
                                   snets_midx,
                                   snets_mp) )
    {
      ERROR_MSG("Problem encountered during ExecMainSpecnets stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-100);
    }

    delete psms;
    delete psms_spectra;
    delete psms_midx;
    delete psms_mp;
    delete snets_contigs;
    delete snets_midx;
    delete snets_mp;

    return 0;
  }

  Clusters contigShifts;
  abinfo_t contigAbinfo;
  //---------------------------------------------------------------------------
  // ASSEMBLY STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_ASSEMBLY and ip.getValueBool("PARTIAL_OVERLAPS",true))
  {
    if (!performAssembly(ip,
										    starSpectra,
										    filteredPairs,
                        contigShifts,
                        contigAbinfo))
    {
      ERROR_MSG("Problem encountered during Assembly stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_ASSEMBLY);
    }

    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }
    contigShifts.Save("assembly/path_spectra_as_cluster.txt");

  } else {
    DEBUG_MSG("Bypassing Assembly stage" );
    if(contigShifts.Load("assembly/path_spectra_as_cluster.txt") <= 0) {
      ERROR_MSG("Problem encountered while skipping Assembly stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_ASSEMBLY);
    }
    if(contigShifts.consensus.LoadSpecSet_pklbin("assembly/sps_seqs.pklbin") <= 0) {
      ERROR_MSG("Problem encountered while skipping Assembly stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_ASSEMBLY);
    }
  }

  // Test for empty return data structures
  TEST_RETURN("ExecAssembly", contigShifts);

  if(ip.exists("FASTA_DATABASE")) {

	  //---------------------------------------------------------------------------
	  // TAGSEARCH STAGE
	  //---------------------------------------------------------------------------
	  if (initialStage <= STAGE_TAGSEARCH)
	  {
		bool ok;
		if(ip.getValueBool("PARTIAL_OVERLAPS",true)) {
			ok = performTagsearch(ip,
									contigShifts.consensus,
									dbAll,
									(vector<unsigned int> *) 0);  //specsToSearch,

		} else {
			ok = performTagsearch(ip,
									contigShifts.consensus,
									dbAll,
									&alignedSpectra);
		}
		if(!ok)
		{
		  ERROR_MSG("Problem encountered during TagSearch stage" );
		  writeStatusFile(statusFileName,"Error");
		  exit(-STAGE_TAGSEARCH);
		}

		DEBUG_VAR(contigShifts.size());

		if (commandLineParams.exists("SINGLE_STEP"))
		{
		  DEBUG_MSG("Option -s given. Exiting after single step.");
		  writeStatusFile(statusFileName,"Finished");
		  exit(0);
		}
	  } else {
		// No need to load data here.. since we always have to reload after TagSearch
	  }
	  DEBUG_MSG("Reloading contig specturm" );
	  // Need to reload the contigs because they were altered by TagSearch
		if(contigShifts.Load(getProjPath(ip, "assembly/path_spectra_as_cluster.txt").c_str(),
							 getProjPath(ip, "assembly/sps_seqs.pklbin").c_str(),
							 getProjPath(ip, "assembly/tagsearchpsm.txt").c_str()) <= 0) {
		ERROR_MSG("Problem encountered while skipping Assembly stage" );
		writeStatusFile(statusFileName,"Error");
		exit(-1);
		}

	  DEBUG_VAR(contigShifts.size());
	  DEBUG_VAR(contigShifts.consensus.size());


	  // Test for empty return data structures
	  //TEST_RETURN("ExecTagSearch", alignedSpectra);

	  //---------------------------------------------------------------------------
	  // Spectrum/Protein Alignment STAGE
	  //---------------------------------------------------------------------------
	  SpecSet matchedContigs; // Matched versions of spectra significantly matched to a protein
	  vector<unsigned int> matchedContigIndices; // Indices of spectra significantly matched to a protein

	  if (initialStage <= STAGE_SPECPROTALIGN)
	  {
		DEBUG_VAR(contigShifts.size());
		if (!performSpecProtAlign(ip,
									  contigShifts,
									  dbAll,
								  &matchedContigs,
								  matchedContigIndices))
		{
		  ERROR_MSG("Problem encountered during SpecProtAlign stage" );
		  writeStatusFile(statusFileName,"Error");
		  exit(-STAGE_SPECPROTALIGN);
		}

		if (commandLineParams.exists("SINGLE_STEP"))
		{
		  DEBUG_MSG("Option -s given. Exiting after single step.");
		  writeStatusFile(statusFileName,"Finished");
		  exit(0);
		}
	  } else {
		DEBUG_MSG("Bypassing Spectrum/Protein Alignment stage" );
		if (matchedContigs.loadPklBin(getProjPath(ip, "spectra/contigs.pklbin").c_str(),
									  getProjPath(ip, "homology/contigs_psm.txt").c_str(),
									  getProjPath(ip, "homology/contigs_midx.pklbin").c_str()) <= 0) {
		  ERROR_MSG("Problem encountered while skipping SpecProtAlign stage" );
		  writeStatusFile(statusFileName,"Error");
		  exit(-STAGE_SPECPROTALIGN);
		}
	  }

	  //SpecSet matchedContigs; // Matched versions of spectra significantly matched to a protein
	  //SpecSet matchedMassIndices; // Per-spectrum indices of matched spectrum/protein masses

	  // Test for empty return data structures
	  //TEST_RETURN("ExecSpecProtAlign", matchedContigs);
	  //TEST_RETURN("ExecSpecProtAlign", matchedContigIndices);


	  //---------------------------------------------------------------------------
	  // Find reference protein for protein/protein sequence alignments
	  //---------------------------------------------------------------------------
		unsigned int refProtIdx=0, refProtMatchCount=0;
	  set<unsigned int> dbIndexes;

	  //---------------------------------------------------------------------------
	  // Assign the reference protein for protein/protein sequence alignments
	  //---------------------------------------------------------------------------
	  if (ip.exists("FORCE_REFERENCE")) {
		// Force use of the supplied reference
		refProtIdx = ip.getValueInt("FORCE_REFERENCE");

		for (unsigned int i = 0; i < matchedContigs.size(); i++) {
		  dbIndexes.insert(matchedContigs[i].psmList.front()->m_dbIndex);
		}
		DEBUG_VAR(refProtIdx);

	  } else {
		// Use the most common protein match as the reference
		vector<unsigned int> numMatches(dbAll.size());
		for (unsigned int i = 0; i < dbAll.size(); i++)
		  numMatches[i]=0;

		for (unsigned int i = 0; i < matchedContigs.size(); i++) {
		  int index = matchedContigs[i].psmList.front()->m_dbIndex;
		  dbIndexes.insert(index);
			if(index >= 0 and index < dbAll.size()) {
	  //  	numMatches[index]++;  // Select by highest number of matched contigs
				numMatches[index] += matchedContigs[i].psmList.front()->m_matchedPeaks.size();  // Select by highest number of matched masses
				if(numMatches[index] > refProtMatchCount) {
					refProtIdx = index;
					refProtMatchCount = numMatches[index];
				}
			}
		}
		DEBUG_VAR(refProtIdx);
		DEBUG_VAR(refProtMatchCount);

	  } // if (ip.exists("FORCE_REFERENCE")) {


	  //---------------------------------------------------------------------------
	  // Protein/Protein Alignment STAGE
	  //---------------------------------------------------------------------------
	  DEBUG_TRACE;
	  string indexFileName; // File names of clustalw .aln output files
	  if (initialStage <= STAGE_PROTPROTALIGN)
	  {
		if (!performProtProtAlign(ip, refProtIdx, dbIndexes, dbAll))
		{
		  ERROR_MSG("Problem encountered during ProtProtAlign stage" );
		  writeStatusFile(statusFileName,"Error");
		  exit(-STAGE_PROTPROTALIGN);
		}

		if (commandLineParams.exists("SINGLE_STEP"))
		{
		  DEBUG_MSG("Option -s given. Exiting after single step.");
		  writeStatusFile(statusFileName,"Finished");
		  exit(0);
		}
	  } else {
		DEBUG_MSG("Bypassing Protein/Protein Alignment stage" );
	  }

	  DEBUG_VAR(dbAll.size());
	  DEBUG_VAR(matchedContigs.size());

	  //---------------------------------------------------------------------------
	  // Homology-based Assembly STAGE
	  //---------------------------------------------------------------------------
	  SpecSet cspsContigs;
	  SpecSet cspsMatchedIndices;
	  vector<vector<int> > cspsMatchedProts;
	  if (initialStage <= STAGE_HOMOLOGYASSEMBLY)
	  {
		DEBUG_TRACE;
		if (!performHomologyAssembly(ip,
												   matchedContigs,
															   dbAll,
															   cspsContigs,
															   cspsMatchedIndices,
														   cspsMatchedProts))
		{
		  ERROR_MSG("Problem encountered during cSPS HomologyAssembly stage" );
		  writeStatusFile(statusFileName,"Error");
		  exit(-STAGE_HOMOLOGYASSEMBLY);
		}

		if (commandLineParams.exists("SINGLE_STEP"))
		{
		  DEBUG_MSG("Option -s given. Exiting after single step.");
		  writeStatusFile(statusFileName,"Finished");
		  exit(0);
		}
	  } else {
		DEBUG_MSG("Bypassing Homology Assembly stage" );
		if(cspsContigs.LoadSpecSet_pklbin(getProjPath(ip, "homology/homglue_matches.pklbin").c_str() ) <= 0 or
				cspsMatchedIndices.LoadSpecSet_pklbin(getProjPath(ip, "homology/homglue_matches_midx.pklbin").c_str() ) <= 0 or
				Load_binArray<int>(getProjPath(ip, "homology/homglue_matches_mp.bin").c_str(), cspsMatchedProts) <=0 )
		{
		  ERROR_MSG("Problem encountered while skipping SpecProtAlign stage" );
		  writeStatusFile(statusFileName,"Error");
		  exit(-STAGE_HOMOLOGYASSEMBLY);
		}
	  }

	  // Test for empty return data structures
	  TEST_RETURN("ExecHomologyAssembly", cspsContigs);
	  TEST_RETURN("ExecHomologyAssembly", cspsMatchedIndices);
	  TEST_RETURN("ExecHomologyAssembly", cspsMatchedProts);

	  //---------------------------------------------------------------------------
	  // MERGE RESULTS IF RUNNING IN INTEGRATIVE MODE
	  //---------------------------------------------------------------------------
	  if(runMergedFlag && initialStage <= STAGE_HOMOLOGYASSEMBLY)
	  {
		//We need to rename all of our results files
		if(rename("./spectra/stars.pklbin","./spectra/csps.stars.pklbin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }
		if(rename("./assembly/sps_seqs.pklbin","./assembly/csps.sps_seqs.pklbin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }
		if(rename("./assembly/component_info.bin","./assembly/csps.component_info.bin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }
		if(rename("./homology/contigs_mp.bin","./homology/csps.contigs_mp.bin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }
		if(rename("./homology/contigs_midx.pklbin","./homology/csps.contigs_midx.pklbin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }
		if(rename("./homology/contigs_mp_all.bin","./homology/csps.contigs_mp_all.bin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./homology/contigs_midx_all.pklbin","./homology/csps.contigs_midx_all.pklbin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./spectra/contigs.pklbin","./spectra/csps.contigs.pklbin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./spectra/contigs_indices.bin","./spectra/csps.contigs_indices.bin") < 0)	  {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./homology/homglue_ref_mp.bin","./homology/csps.homglue_ref_mp.bin") < 0)	  {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./homology/homglue_ref_midx.pklbin","./homology/csps.homglue_ref_midx.pklbin") < 0) {
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./homology/homglue_matches.pklbin","./homology/csps.homglue_matches.pklbin") < 0)	{
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./homology/homglue_matches_mp.bin","./homology/csps.homglue_matches_mp.bin") < 0)	{
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./homology/homglue_matches_midx.pklbin","./homology/csps.homglue_matches_midx.pklbin") < 0)	{
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

		if(rename("./homology/ref_sps_names.txt","./homology/csps.ref_sps_names.txt") < 0)	{
			ERROR_MSG("Problem encountered renaming CSPS files");
			writeStatusFile(statusFileName,"Error");
			exit(-STAGE_MERGE);
		  }

	  } // if(runMergedFlag && initialStage <= STAGE_HOMOLOGYASSEMBLY)

	  if (runMergedFlag && initialStage <= STAGE_MERGE) {
		if(!performMergeOfCSPSAndGenoMS(ip)) {
		  ERROR_MSG("Problem encounctered during merge of CSPS and GenoMS");
		  writeStatusFile(statusFileName,"Error");
		  exit(-STAGE_MERGE);
		}
	  }
  }
  else // if(ip.exists("FASTA_DATABASE"))
  {
//
//  --- need to fix performReports to allow reporting de novo results only
//
//	WARN_MSG("No FASTA_DATABASE specified, reporting only de novo results.\n");
	WARN_MSG("No FASTA_DATABASE specified, returning without database matches.\n");
	return 0;
  }

  //---------------------------------------------------------------------------
  // REPORT STAGE
  //---------------------------------------------------------------------------
  if (initialStage <= STAGE_REPORT)
  {
  	ofstream spsProj(getProjPath(ip, "sps_projects.txt").c_str(), ios_base::out);
  	spsProj << "sps;.;"<<ip.getValue("TOLERANCE_PEAK")<<";"<<ip.getValue("TOLERANCE_PM")<<"\n";
  	spsProj.close();


    if (!performReport(ip))
    {
      ERROR_MSG("Problem encountered during Report stage" );
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_REPORT);
    }
    if (commandLineParams.exists("SINGLE_STEP"))
    {
      DEBUG_MSG("Option -s given. Exiting after single step.");
      writeStatusFile(statusFileName,"Finished");
      exit(0);
    }

  }

  //--------------------------------------------------------------------------
  //Set up for relaunching
  //--------------------------------------------------------------------------
  if(!generateRelaunchScript(ip))
    {
      ERROR_MSG("Problem encountered during relaunch script creation");
      writeStatusFile(statusFileName,"Error");
      exit(-STAGE_REPORT);
    }
  //---------------------------------------------------------------------------
  // END
  //---------------------------------------------------------------------------
  writeStatusFile(statusFileName,"Finished");
  return 0;
} // END

