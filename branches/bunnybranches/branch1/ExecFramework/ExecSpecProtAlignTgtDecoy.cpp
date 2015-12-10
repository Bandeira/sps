#include "ExecSpecProtAlignTgtDecoy.h"

#include "alignment_modmut.h"
#include "AlignmentPenaltyBased.h"
#include "FdrPeptide.h"
#include "FileUtils.h"
#include "Logger.h"
#include "ParallelPbsExecution.h"
#include "ParallelSgeExecution.h"

// SpecNets Includes
#include "tags.h"
#include <limits.h>


#define DEBUG_SPECPROTALIGN 0

using namespace std;
using namespace specnets;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(void) :
        m_inputSpectra(0), m_db(0x0), m_dbDecoy(0x0), 
        m_psmSetTag(0x0), m_psmSetTagDecoy(0x0), ownInput(true),
        m_matchedSpectraAll(0x0), 
        m_psmSet(0x0), m_psmSetDecoy(0x0), m_psmSetFdr(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputSpectra(0x0), 
        m_db(0x0), m_dbDecoy(0x0), m_psmSetTag(0x0), m_psmSetTagDecoy(0x0), ownInput(true),
        m_matchedSpectraAll(0x0), 
        m_psmSet(0x0), m_psmSetDecoy(0x0), m_psmSetFdr(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::ExecSpecProtAlignTgtDecoy(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       DB_fasta * db,
                                       DB_fasta * dbDecoy,
                                       PeptideSpectrumMatchSet * psmSetTag,
                                       PeptideSpectrumMatchSet * psmSetTagDecoy,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods,                                             
                                       SpecSet * matchedSpectraAll,
                                       SpecSet * matchedSpectra,
                                       PeptideSpectrumMatchSet * psmSet,
                                       PeptideSpectrumMatchSet * psmSetDecoy,
                                       PeptideSpectrumMatchSet * psmSetFdr) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), 
        m_db(db), m_dbDecoy(dbDecoy), m_psmSetTag(psmSetTag), m_psmSetTagDecoy(psmSetTagDecoy),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum), 
        m_penaltyMatrixMods(penaltyMatrixMods), ownInput(false),
        m_matchedSpectraAll(matchedSpectraAll), m_matchedSpectra(matchedSpectra),
        m_psmSet(psmSet), m_psmSetDecoy(psmSetDecoy), m_psmSetFdr(psmSetFdr), ownOutput(false)
  {
    m_name = "ExecSpecProtAlignTgtDecoy";
    m_type = "ExecSpecProtAlignTgtDecoy";
  }
                			

  // -------------------------------------------------------------------------

  ExecSpecProtAlignTgtDecoy::~ExecSpecProtAlignTgtDecoy(void)
  {
    if (ownInput) {
      if (m_inputSpectra)
        delete m_inputSpectra;
      if (m_db)
        delete m_db;
      if (m_dbDecoy)
        delete m_dbDecoy;
      if (m_psmSetTag)
        delete m_psmSetTag;
      if (m_psmSetTagDecoy)
        delete m_psmSetTagDecoy;
      if (m_penaltyMatrixBlosum)
        delete m_penaltyMatrixBlosum;
      if (m_penaltyMatrixMods)
        delete m_penaltyMatrixMods;
    }
    if (ownOutput) {
      if (m_matchedSpectraAll)
        delete m_matchedSpectraAll;
      if (m_matchedSpectra)
        delete m_matchedSpectra;
      if (m_psmSet)
        delete m_psmSet;
      if (m_psmSetDecoy)
        delete m_psmSetDecoy;
      if (m_psmSetFdr)
        delete m_psmSetFdr;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecSpecProtAlignTgtDecoy::clone(const ParameterList & inputParams) const
  {
    return new ExecSpecProtAlignTgtDecoy(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::invoke(void)
  {
    if (!m_inputSpectra || m_inputSpectra->size() == 0) {
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

    if (!m_db or m_db->size() == 0) {
      ERROR_MSG("ERROR: empty database");
      return false;
    }

    if (!m_matchedSpectraAll or !m_matchedSpectra) {
      ERROR_MSG("ERROR: empty returned data pointer");
      return false;
    }

    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");
    DEBUG_VAR(peakTol);
    bool gridExecutionFlag = (bool)m_params.getValueInt("GRID_EXECUTION_FLAG", 0);
    DEBUG_VAR(gridExecutionFlag);
    bool resume = (bool)m_params.getValueInt("GRID_RESUME_FLAG", 0);
    DEBUG_VAR(resume );

    m_params.setValue("GRID_DATA_DIR", m_params.getValue("GRID_DATA_DIR_TARGET"));
    m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_TGT"));
    m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", m_params.getValue("OUTPUT_MATCHED_SPECS_TGT"));
    m_params.setValue("OUTPUT_PSM_ALL", m_params.getValue("OUTPUT_PSM_TGT"));

    ExecSpecProtAlign moduleSpecProtAlign(m_params,
                                          m_inputSpectra,
                                          m_db,
                                          m_penaltyMatrixBlosum,
                                          m_penaltyMatrixMods,
                                          m_matchedSpectraAll);

    bool returnStatus;
    if (!m_params.exists("GRID_NUMNODES") or m_params.getValueInt("GRID_NUMNODES") <= 0) {
      returnStatus = moduleSpecProtAlign.invoke();
    } else {
      DEBUG_TRACE;
      int numNodes = m_params.getValueInt("GRID_NUMNODES");

      string gridType = m_params.getValue("GRID_TYPE");
      if (gridType == "pbs") {
        ParallelPbsExecution exec(&moduleSpecProtAlign,
                                  true,
                                  false,
                                  resume);
        returnStatus = exec.invoke(numNodes);
      } else if (gridType == "sge") {
        ParallelSgeExecution exec(&moduleSpecProtAlign,
                                  true,
                                  false,
                                  resume);
        returnStatus = exec.invoke(numNodes);
      }
    }

    returnStatus = moduleSpecProtAlign.saveOutputData();
    DEBUG_VAR(returnStatus);

    m_params.setValue("GRID_DATA_DIR", m_params.getValue("GRID_DATA_DIR_DECOY"));
    m_params.setValue("OUTPUT_MATCHED_PEAKS_IDX_ALL", m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_DEC"));
    m_params.setValue("OUTPUT_MATCHED_SPECS_ALL", m_params.getValue("OUTPUT_MATCHED_SPECS_DEC"));
    m_params.setValue("OUTPUT_PSM_ALL", m_params.getValue("OUTPUT_PSM_DEC"));

    ExecSpecProtAlign moduleSpecProtAlign2(m_params,
                                          m_inputSpectra,
                                          m_dbDecoy,
                                          m_penaltyMatrixBlosum,
                                          m_penaltyMatrixMods,
                                          &m_matchedSpectraAllDecoy);

    if (!m_params.exists("GRID_NUMNODES") or m_params.getValueInt("GRID_NUMNODES") <= 0) {
      returnStatus = moduleSpecProtAlign2.invoke();
    } else {
      DEBUG_TRACE;
      int numNodes = m_params.getValueInt("GRID_NUMNODES");

      string gridType = m_params.getValue("GRID_TYPE");
      if (gridType == "pbs") {
        ParallelPbsExecution exec(&moduleSpecProtAlign2,
                                  true,
                                  false,
                                  resume);
        returnStatus = exec.invoke(numNodes);
      } else if (gridType == "sge") {
        ParallelSgeExecution exec(&moduleSpecProtAlign2,
                                  true,
                                  false,
                                  resume);
        returnStatus = exec.invoke(numNodes);
      }
    }

    returnStatus = moduleSpecProtAlign2.saveOutputData();
    DEBUG_VAR(returnStatus);

    float fdrThreshold = m_params.getValueFloat("FDR_THRESHOLD", -1.0);
    DEBUG_VAR(fdrThreshold);

    if (fdrThreshold >= 0.0) {

      m_psmSet->getPSMSet(m_matchedSpectraAll);
      m_psmSetDecoy->getPSMSet(&m_matchedSpectraAllDecoy);

      PeptideSpectrumMatchSet psmSetCombined;
      if (!FdrPeptide::mergeTargetDecoy(*m_psmSet, *m_psmSetDecoy, psmSetCombined)) {
        ERROR_MSG("FDR concatenation failed.");
      }
      if (DEBUG_SPECPROTALIGN) psmSetCombined.saveToFile("debug_psm_combined.txt");
      DEBUG_VAR(psmSetCombined.size());

      if (!FdrPeptide::concatenatedTargetDecoy(psmSetCombined, *m_psmSetFdr)) {
        ERROR_MSG("FDR sorting failed.");
      }
      if (DEBUG_SPECPROTALIGN) m_psmSetFdr->saveToFile("debug_psm_sorted.txt");
      DEBUG_VAR(m_psmSetFdr->size());

      double fdrCutoff  = 1.0;
      if (m_params.getValueFloat("FDR_THRESHOLD", -1.0) >= 0.0) {
        fdrCutoff = m_params.getValueFloat("FDR_THRESHOLD", 1.0);  // Purposely 1.0 not -1.0 (should exist anyway)
      }
      DEBUG_VAR(fdrCutoff);

      if (!FdrPeptide::filterByPValue(*m_psmSetFdr, fdrCutoff)) {
        ERROR_MSG("Filter by FDR value failed.");
      }
      if (DEBUG_SPECPROTALIGN) m_psmSetFdr->saveToFile("debug_psm_cutoff.txt");
      DEBUG_VAR(m_psmSetFdr->size());
    }

    // Resize to make enough room
    m_matchedSpectraIndices.resize(m_matchedSpectraAll->size());
    m_matchedSpectra->resize(m_matchedSpectraAll->size());

    int keepIdx = 0;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(i);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].size());
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
      if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
        m_matchedSpectraIndices[keepIdx] = i;
        (*m_matchedSpectra)[keepIdx++] = (*m_matchedSpectraAll)[i];
      }
    }

    DEBUG_VAR(keepIdx);
    m_matchedSpectra->resize(keepIdx);
    m_matchedSpectraIndices.resize(keepIdx);

    DEBUG_VAR(m_matchedSpectra->size());
    DEBUG_VAR(m_matchedSpectraIndices.size());

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::loadInputData(void)
  {
    if (ownInput) {
      if (!m_inputSpectra)
        m_inputSpectra = new SpecSet;
      if (!m_db)
        m_db = new DB_fasta;
      if (!m_dbDecoy)
        m_dbDecoy = new DB_fasta;
      if (!m_psmSetTag)
        m_psmSetTag = new PeptideSpectrumMatchSet;
      if (!m_psmSetTagDecoy)
        m_psmSetTagDecoy = new PeptideSpectrumMatchSet;
    }
    m_inputSpectra->resize(0);

    if (ownOutput) {
      if (!m_matchedSpectraAll)
        m_matchedSpectraAll = new SpecSet;
      if (!m_psmSet)
        m_psmSet = new PeptideSpectrumMatchSet;
      if (!m_psmSetDecoy)
        m_psmSetDecoy = new PeptideSpectrumMatchSet;
      if (!m_psmSetFdr)
        m_psmSetFdr = new PeptideSpectrumMatchSet;
    }
    m_matchedSpectraAll->resize(0);

    if (m_params.exists("INPUT_SPECS_PKLBIN") and m_params.exists("INPUT_PSM")) {
      if (m_inputSpectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str(),
                                     m_params.getValue("INPUT_PSM").c_str(),
                                     m_params.getValue("INPUT_MATCHED_PEAKS_IDX").c_str()) < 0) {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    if (!m_params.exists("INPUT_FASTA")) {
      ERROR_MSG("Parameters are incomplete. INPUT_FASTA is missing.");
      return false;
    }
    else if (m_db->Load(m_params.getValue("INPUT_FASTA").c_str()) <= 0) {
      ERROR_MSG("Error reading database sequences from "
          << m_params.getValue("INPUT_FASTA"));
      return false;
    }

    if (!m_params.exists("INPUT_FASTA_DECOY")) {
      if (m_dbDecoy->Load(m_params.getValue("INPUT_FASTA_DECOY").c_str()) <= 0) {
        ERROR_MSG("Error reading decoy database sequences from "
            << m_params.getValue("INPUT_FASTA_DECOY"));
        return false;
      }
    }

    if (m_params.exists("INPUT_PSM")) {
      m_psmSetTag->loadFromFile(m_params.getValue("INPUT_PSM").c_str());
    }

    if (m_params.exists("INPUT_PSM_DECOY")) {
      m_psmSetTagDecoy->loadFromFile(m_params.getValue("INPUT_PSM_DECOY").c_str());
    }

    //---------------------------------------------------------------------------
    // Load penalty matrices for new alignment
    //---------------------------------------------------------------------------

    bool penaltyAlign = (bool) m_params.exists("PENALTY_ALIGNMENT") ?
    		((int) m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;

    if (penaltyAlign) {

      // Load amino acid masses
      AAJumps jumps(1);
      if (!m_params.exists("AMINO_ACID_MASSES")) {
        jumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true);
      }

      float resolution = m_params.getValueDouble("RESOLUTION", 0.1);

      if (!m_penaltyMatrixBlosum)
        m_penaltyMatrixBlosum = new PenaltyMatrix(jumps, resolution);
      if (!m_penaltyMatrixMods)
        m_penaltyMatrixMods = new PenaltyMatrix(jumps, resolution);

      DEBUG_TRACE;
      if (!m_params.exists("BLOSUM_PENALTY_FILE")) {
        ERROR_MSG("Parameters are incomplete. BLOSUM_FILE is missing and penalty alignment was specified.");
        return false;
      } else {
        string blosumFileName = m_params.getValue("BLOSUM_PENALTY_FILE"); 
        if (!m_penaltyMatrixBlosum->load(blosumFileName) ) {
          ERROR_MSG("Error loading blosum penalties from " << blosumFileName);
          return false;
        }
      }

      DEBUG_TRACE;
      if (!m_params.exists("KNOWN_MODS_FILE")) {
        ERROR_MSG("Parameters are incomplete. KNOWN_MODS_FILE is missing and penalty alignment was specified.");
        return false;
      } else {
        string knowmModsFileName = m_params.getValue("KNOWN_MODS_FILE"); 
        m_penaltyMatrixMods->loadKnownModifications(knowmModsFileName);
      }

      DEBUG_TRACE;
      if (!m_params.exists("MODS_PENALTY_FILE")) {
        ERROR_MSG("Parameters are incomplete. MODS_PENALTY_FILE is missing and penalty alignment was specified.");
        return false;
      } else {
        string modFileName = m_params.getValue("MODS_PENALTY_FILE"); 
        if (!m_penaltyMatrixMods->load(modFileName) ) {
          ERROR_MSG("Error loading mod penalties from " << modFileName);
          return false;
        }
      }

    } // if (penaltyAlign)

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::saveOutputData(void)
  {
    bool penaltyAlign = (bool) m_params.exists("PENALTY_ALIGNMENT") ?
  	  	((int) m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;
    DEBUG_VAR(penaltyAlign);

    if (m_params.exists("OUTPUT_MATCHED_PEAKS_IDX")) {
      int keepIdx = 0;
      for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
        if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
          (*m_matchedSpectra)[keepIdx++] = (*m_matchedSpectraAll)[i];
        }
      }

      DEBUG_MSG("Saving matched peaks...");
      SpecSet tempMatchedPeaks;
      tempMatchedPeaks.resize(m_matchedSpectra->size());
      for (int i = 0; i < m_matchedSpectra->size(); i++) {
        if ((*m_matchedSpectra)[i].psmList.size() != 0) {
          list<psmPtr>::iterator litr = (*m_matchedSpectra)[i].psmList.begin();
          int peakListSize = (*litr)->m_matchedPeaks.size();
          //DEBUG_VAR(peakListSize);
          tempMatchedPeaks[i].resize(peakListSize);
          for (int j = 0; j < peakListSize; j++) {
            tempMatchedPeaks[i][j].set((*litr)->m_matchedPeaks[j][0],
                                       (*litr)->m_matchedPeaks[j][1]);
            //DEBUG_MSG((*litr)->m_matchedPeaks[j][0] << "  " << (*litr)->m_matchedPeaks[j][1]);
          }
        }
      }
      tempMatchedPeaks.SaveSpecSet_pklbin(m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX").c_str());
    }

    DEBUG_VAR(m_matchedSpectra->size());
    if (m_matchedSpectra and m_params.exists("OUTPUT_MATCHED_SPECS")) {
      if (m_params.exists("OUTPUT_PSM") && m_params.exists("OUTPUT_MATCHED_PEAKS_IDX")) {
        DEBUG_MSG("Saving matched specs (3 files)...");
        m_matchedSpectra->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS").c_str(),
                                     m_params.getValue("OUTPUT_PSM").c_str());

      } else if (m_params.exists("OUTPUT_PSM")) {
        DEBUG_MSG("Saving matched specs (2 files)...");
        m_matchedSpectra->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS").c_str(),
                                     m_params.getValue("OUTPUT_PSM").c_str());

      } else {
        DEBUG_MSG("Saving matched specs (1 file)...");
        m_matchedSpectra->SaveSpecSet_pklbin(m_params.getValue("OUTPUT_MATCHED_SPECS").c_str());
      }
    }

    //=================================================
    // LARS - TEMPORARY MEASURE SO REPORTS KEEP WORKING
    //=================================================
    DEBUG_VAR(m_matchedSpectraAll->size());
    if (m_params.exists("OUTPUT_MATCHED_PROTS_ALL")) {
      DEBUG_MSG("Saving matched prots all...");
      m_matchedSpectraAll->saveMatchedProts(m_params.getValue("OUTPUT_MATCHED_PROTS_ALL").c_str());
    }
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_MATCHED_PROTS")) {
      DEBUG_MSG("Saving matched prots...");
      m_matchedSpectra->saveMatchedProts(m_params.getValue("OUTPUT_MATCHED_PROTS").c_str());
    }
    DEBUG_TRACE;
    if (m_params.exists("OUTPUT_MATCHED_SPECS_IDX")) {
      DEBUG_MSG("Saving matched specs idx...");
      Save_binArray(m_params.getValue("OUTPUT_MATCHED_SPECS_IDX").c_str(),
                    m_matchedSpectraIndices);
    }
    DEBUG_TRACE;
    //=================================================

    if (m_params.exists("OUTPUT_PSM_MOD_FREQS")) {
      PeptideSpectrumMatchSet psmSetTemp;
      psmSetTemp.getPSMSet(m_matchedSpectraAll);
      psmSetTemp.saveModMatrix(m_params.getValue("OUTPUT_PSM_MOD_FREQS").c_str(), penaltyAlign);
    }

    if (m_params.exists("OUTPUT_PSM_FDR")) {
      m_psmSetFdr->saveToFile(m_params.getValue("OUTPUT_PSM_FDR").c_str());
    }

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecSpecProtAlignTgtDecoy::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::merge(void)
  {
    return false;
  }


  // -------------------------------------------------------------------------

  bool ExecSpecProtAlignTgtDecoy::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");

    m_isValid = true;
    return true;
  }

// -------------------------------------------------------------------------


}

