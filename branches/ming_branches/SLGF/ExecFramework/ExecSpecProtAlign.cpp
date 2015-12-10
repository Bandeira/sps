#include "ExecSpecProtAlign.h"

// Module Includes
#include "Logger.h"
//#include "FileUtils.h"

// SpecNets Includes
#include "tags.h"
#include <limits.h>

using namespace std;
using namespace specnets;

namespace specnets
{
  class LessThanPredicate
  {
  public:
    LessThanPredicate(float value)
      : m_value(value)
      {
      }
    bool operator() (const psmPtr & value) { return value->m_score < m_value; }
  private:
    float m_value;
  };

  
  bool PsmDbIndexSort(psmPtr left, psmPtr right) {
    return left->m_dbIndex < right->m_dbIndex;
  }

  bool PsmDbIndexUnique(psmPtr left, psmPtr right) { 
    return left->m_dbIndex == right->m_dbIndex;
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::ExecSpecProtAlign(void) :
    m_contigs(0), m_db(0x0), ownInput(true),
        m_matchedSpectra(0x0), m_matchedSpectraIndices(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams) :
    ExecBase(inputParams), m_contigs(0), m_db(0x0), ownInput(true),
        m_matchedSpectra(0x0), m_matchedSpectraIndices(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams,
                                       Clusters * contigs,
                                       DB_fasta * db,
                                       SpecSet * matchedSpectra,
                                       vector<unsigned int> * matchedSpectraIndices) :
    ExecBase(inputParams), m_contigs(contigs), m_db(db), ownInput(false), 
        m_matchedSpectra(matchedSpectra),
        m_matchedSpectraIndices(matchedSpectraIndices), ownOutput(false)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";

    if (m_matchedSpectraIndices)
      m_matchedSpectraIndices->resize(m_contigs->size());
  }
                			

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::~ExecSpecProtAlign(void)
  {
    if (ownInput) {
      if (m_contigs)
        delete m_contigs;
    }
    if (ownInput) {
      if (m_db)
        delete m_db;
    }
    if (ownOutput) {
      if (m_matchedSpectra)
        delete m_matchedSpectra;
      if (m_matchedSpectraIndices)
        delete m_matchedSpectraIndices;
    }
  }

  // -------------------------------------------------------------------------

  ExecBase * ExecSpecProtAlign::clone(const ParameterList & inputParams) const
  {
    return new ExecSpecProtAlign(inputParams);
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::invoke(void)
  {
    if (!m_contigs or m_contigs->size() == 0) {
      DEBUG_VAR(m_contigs);
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

    if (!m_db or m_db->size() == 0) {
      ERROR_MSG("ERROR: empty database");
      return false;
    }

    if (!m_matchedSpectra or !m_matchedSpectraIndices) {
      ERROR_MSG("ERROR: empty returned data pointer");
      return false;
    }

    float pmTol = (float) m_params.getValueDouble("TOLERANCE_PM");
    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");

    unsigned int startIdx = m_params.exists("IDX_START") ? 
            max(0, m_params.getValueInt("IDX_START")) : 
            0;
    unsigned int endIdx = m_params.exists("IDX_END") ? 
            min(m_contigs->size() - 1, (unsigned int) m_params.getValueInt("IDX_END")) : 
            m_contigs->size() - 1;
    unsigned int startIdx_db = m_params.exists("IDX_START_DB") ? 
            max(0, m_params.getValueInt("IDX_START_DB")) : 
            0;
    unsigned int endIdx_db = m_params.exists("IDX_END_DB") ? 
            min(m_db->size() - 1, (unsigned int) m_params.getValueInt("IDX_END_DB")) : 
            m_db->size() - 1;
    if (startIdx > m_contigs->size())
      return false;
    if (startIdx_db > m_db->size())
      return false;

    bool enforceEndpeaks = m_params.exists("ENFORCE_ENDPEAKS") ? 
            m_params.getValueInt("ENFORCE_ENDPEAKS") > 0 : false;
    float maxModMass = m_params.exists("MAX_MOD_MASS") ? 
            (float) max(0.0, m_params.getValueDouble("MAX_MOD_MASS")) : 372.2;
    float minModMass = m_params.exists("MIN_MOD_MASS") ? 
            (float) max(-372.2, m_params.getValueDouble("MIN_MOD_MASS")) : -372.2;
    int showMatchedPRMs = m_params.exists("SHOW_MATCHED_PRMS") ? 
            m_params.getValueInt("SHOW_MATCHED_PRMS") : 0;
    int maxNumMods = m_params.exists("MAX_NUM_MODS") ? 
            max(0, m_params.getValueInt("MAX_NUM_MODS")) : 1;
    int minNumMatchPeaks = m_params.exists("MIN_NUM_MATCH_PEAKS") ? 
            max(0, m_params.getValueInt("MIN_NUM_MATCH_PEAKS")) : 4;
    int matchesPerMod = m_params.exists("MATCHES_PER_MOD") ? 
            max(0, m_params.getValueInt("MATCHES_PER_MOD")) : 2;
    int specType = m_params.exists("SPEC_TYPE_MSMS") ? 
            m_params.getValueInt("SPEC_TYPE_MSMS") > 0 : 0;
    bool tagParsimony =
            (bool) m_params.exists("MAX_PARSIMONY") ? ((int) m_params.getValueInt("MAX_PARSIMONY") ? 1 : 0) : 0;
    DEBUG_VAR(tagParsimony);
    float ionOffset = specType ? AAJumps::massHion : 0;
    int modIdx, specIdx, protIdx;

    Spectrum tmpSpec;
    tmpSpec.reserve(1024);
    Spectrum cSpec;
    cSpec.reserve(1024);
    Spectrum cSpecRev;
    cSpecRev.reserve(1024);
    AMM_match *bestMatch;
    
    DEBUG_VAR(m_contigs->size());

    for (specIdx = startIdx; specIdx <= endIdx; specIdx++) {
      DEBUG_VAR(specIdx);
      
      //DEBUG_VAR(m_contigs->consensus[specIdx].size());
      
      // Apparently it is possible to have a contig that is empty
      if (m_contigs->consensus[specIdx].size() == 0) {
        m_matchedSpectra->push_back(m_contigs->consensus[specIdx]);
        continue;
      }
      DEBUG_VAR(cSpec.size());
      DEBUG_VAR(peakTol);
      m_contigs->getSpecIfB(specIdx, cSpec, peakTol);
      cSpec.psmList.clear();
      DEBUG_VAR(cSpec.size());

      m_contigs->getSpecIfY(specIdx, cSpecRev, peakTol, ionOffset);
      cSpecRev.psmList.clear();
      DEBUG_VAR(cSpecRev.size());
      
      // Set the reverse scan number (not set by getSpecIfY)
      cSpecRev.scan = cSpec.scan;

      map<int, set<float> > matchedProtMap;
      if (m_contigs->consensus[specIdx].psmList.size() != 0) {
        // Coalesce the starting positions on all proteins
        list<psmPtr>::iterator itr = m_contigs->consensus[specIdx].psmList.begin();
        list<psmPtr>::iterator itr_end = m_contigs->consensus[specIdx].psmList.end();
        for (; itr != itr_end; itr++) {
          int idx = (*itr)->m_dbIndex;
          // Only use the forward orientation tags(if tag parsimony was not done)
          if (tagParsimony || (*itr)->m_matchOrientation == 0) {
            if ( matchedProtMap.find(idx) == matchedProtMap.end()) {
              set<float> newList;
              matchedProtMap[idx] = newList;
            }
            matchedProtMap[idx].insert((*itr)->m_startMass);
          }
        } // for (; itr != itrEnd; itr++) {
      } else {
         // If there are no PSM's insert dummy start mass for every db idx
         // This will cause alignment to every protein in database
        for (int idx = 0; idx < m_db->size(); idx++) {
          matchedProtMap[idx].insert(0.0);
        }
      }
            
      DEBUG_MSG("Matching as b...");
      //cSpec.output(cerr);

      //itrMap = matchedProtMap.begin();
      map<int, set<float> >::iterator itrMap = matchedProtMap.begin();
      map<int, set<float> >::iterator itrMapEnd = matchedProtMap.end();
      for (; itrMap !=  itrMapEnd; itrMap++) {
        int protIdx = itrMap->first;

        // If parsimony was done in TagSearch then don't use start positions
        if (tagParsimony) {
          itrMap->second.clear();
        }
        
        //DEBUG_VAR(protIdx);
        scoreOverlapAMM(cSpec,
                        m_db->getMassesSpec(protIdx),
                        protIdx,
                        0,
                        itrMap->second,
                        maxNumMods,
                        minNumMatchPeaks,
                        pmTol,
                        peakTol,
                        57,
                        maxModMass,
                        minModMass,
                        enforceEndpeaks);
      } // for (; itrMap !=  itrMapEnd; itrMap++) {

      DEBUG_MSG("Matching as y...");
      //cSpecRev.output(cerr);

      matchedProtMap.clear();
      // Coalesce the starting positions on all proteins
      if (m_contigs->consensus[specIdx].psmList.size() != 0) {
        list<psmPtr>::iterator itr = m_contigs->consensus[specIdx].psmList.begin();
        list<psmPtr>::iterator itr_end = m_contigs->consensus[specIdx].psmList.end();
        for (; itr != itr_end; itr++) {
          int idx = (*itr)->m_dbIndex;
          // Only use the reverse orientation tags (if tag parsimony was not done)
          if (tagParsimony || (*itr)->m_matchOrientation == 1) {
            if ( matchedProtMap.find(idx) == matchedProtMap.end()) {
              set<float> newList;
              matchedProtMap[idx] = newList;
            }
            matchedProtMap[idx].insert((*itr)->m_startMass);
          }
        } // for (; itr != itrEnd; itr++) {
      } else {
         // If there are no PSM's insert dummy start mass for every db idx
         // This will cause alignment to every protein in database
        for (int idx = 0; idx < m_db->size(); idx++) {
          matchedProtMap[idx].insert(0.0);
        }
      }
      
      itrMap = matchedProtMap.begin();
      for (; itrMap !=  itrMapEnd; itrMap++) {
        int protIdx = itrMap->first;

        // If parsimony was done in TagSearch then don't use start positions
        if (tagParsimony) {
          itrMap->second.clear();
        }
        
        //DEBUG_VAR(protIdx);
        scoreOverlapAMM(cSpecRev,
                        m_db->getMassesSpec(protIdx),
                        protIdx,
                        1,
                        itrMap->second,
                        maxNumMods,
                        minNumMatchPeaks,
                        pmTol,
                        peakTol,
                        57,
                        maxModMass,
                        minModMass,
                        enforceEndpeaks);
                        
      } // for (; itr != itrEnd; itr++) {

      DEBUG_MSG("Matching Complete");

      // Find the top score 
      float bestScore = -(float)INT_MAX;
      DEBUG_VAR(cSpec.psmList.size());
      list<psmPtr>::iterator litr = cSpec.psmList.begin();
      list<psmPtr>::iterator litrEnd = cSpec.psmList.end();
      for (; litr != litrEnd; litr++) {
        if ((*litr)->m_score > bestScore)
          bestScore = (*litr)->m_score;
      }
      DEBUG_VAR(bestScore);
      
      DEBUG_VAR(cSpecRev.psmList.size());
      litr = cSpecRev.psmList.begin();
      litrEnd = cSpecRev.psmList.end();
      for (; litr != litrEnd; litr++) {
        if ((*litr)->m_score > bestScore)
          bestScore = (*litr)->m_score;
      }
      DEBUG_VAR(bestScore);
      
      cSpec.psmList.remove_if(LessThanPredicate(bestScore));
      cSpecRev.psmList.remove_if(LessThanPredicate(bestScore));
      if (cSpecRev.psmList.size() == 0) {
        m_matchedSpectra->push_back(cSpec);
      } else {
        m_matchedSpectra->push_back(cSpecRev);
        //DEBUG_VAR(m_matchedSpectra->back().size());
      }

    } // for (specIdx = startIdx; specIdx <= endIdx; specIdx++)

    DEBUG_VAR(m_matchedSpectra->size());

    // We could have multiple hits on one protein
    // So we perform a sort() and unique() to eliminate them
    for (int i = 0; i < m_matchedSpectra->size(); i++) {
      (*m_matchedSpectra)[i].psmList.sort(PsmDbIndexSort);
      (*m_matchedSpectra)[i].psmList.unique(PsmDbIndexUnique);
    }
    
    // Set backwards pointers/scans for PSM--->Spectrum link
    list<psmPtr>::iterator iter;
    for (int i = 0; i < m_matchedSpectra->size(); i++) {
      (*m_matchedSpectra)[i].scan = i + 1;
      for (iter = (*m_matchedSpectra)[i].psmList.begin(); 
           iter != (*m_matchedSpectra)[i].psmList.end(); 
           iter++) {
        (*iter)->m_spectrum = &(*m_matchedSpectra)[i];
        (*iter)->m_scanNum = (*m_matchedSpectra)[i].scan;
      }
    }

    DEBUG_TRACE;
    
    m_matchedSpectra->maximumParsimony();   

    DEBUG_TRACE;
    
    
    //=================================================
    // LARS - TEMPORARY MEASURE SO REPORTS KEEP WORKING
    //=================================================
    if (m_params.exists("OUTPUT_MATCHED_PROTS_ALL")) {
      m_matchedSpectra->saveMatchedProts(m_params.getValue("OUTPUT_MATCHED_PROTS_ALL").c_str());
    }                   
    //=================================================

    DEBUG_TRACE;
    
    if (m_params.exists("OUTPUT_MATCHED_PEAKS_IDX_ALL")) {
      SpecSet tempMatchedPeaks;
      tempMatchedPeaks.resize(m_matchedSpectra->size());
      for (int i = 0; i < m_matchedSpectra->size(); i++) {
        if ((*m_matchedSpectra)[i].psmList.size() != 0) {
          list<psmPtr>::iterator litr = (*m_matchedSpectra)[i].psmList.begin();
          int peakListSize = (*litr)->m_matchedPeaks.size();
          tempMatchedPeaks[i].resize(peakListSize);
          for (int j = 0; j < peakListSize; j++) {
            tempMatchedPeaks[i][j].set((*litr)->m_matchedPeaks[j][0],
                                       (*litr)->m_matchedPeaks[j][1]);
            //DEBUG_MSG((*litr)->m_matchedPeaks[j][0] << "  " << (*litr)->m_matchedPeaks[j][1]);
          }
        }
      }
      tempMatchedPeaks.SaveSpecSet_pklbin(m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str());
    }

    DEBUG_TRACE;

    // Resize to make enough room
    m_matchedSpectraIndices->resize(m_matchedSpectra->size());

    int keepIdx = 0;
    for (int i = 0; i < m_matchedSpectra->size(); i++) {
      if ((*m_matchedSpectra)[i].psmList.size() != 0) {
        (*m_matchedSpectraIndices)[keepIdx] = i;
        (*m_matchedSpectra)[keepIdx++] = (*m_matchedSpectra)[i];
      }
    }
    DEBUG_VAR(keepIdx);
    m_matchedSpectra->resize(keepIdx);
    m_matchedSpectraIndices->resize(keepIdx);

    DEBUG_VAR(m_matchedSpectra->size());
    DEBUG_VAR(m_matchedSpectraIndices->size());
  
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::loadInputData(void)
  {
    if (ownInput) {
      if (!m_contigs)
        m_contigs = new Clusters;
      if (!m_db)
        m_db = new DB_fasta;
    }
    m_contigs->resize(0);

    if (ownOutput) {
      if (!m_matchedSpectra)
        m_matchedSpectra = new SpecSet;
      if (!m_matchedSpectraIndices)
        m_matchedSpectraIndices = new vector<unsigned int> ;
    }
    m_matchedSpectra->resize(0);
    m_matchedSpectraIndices->resize(0);

    if (m_params.exists("AMINO_ACID_MASSES")) {
      AAJumps tmpJumps(-1);
      tmpJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), true); // Set global defaults for amino acid masses
    }

    if (!m_params.exists("INPUT_CLUSTERS")
        and !m_params.exists("INPUT_SPECS_PKLBIN")) {
      ERROR_MSG("Parameters are incomplete. Both INPUT_CLUSTERS and INPUT_SPECS_PKLBIN are missing.");
      return false;
    }
    else {

      if (m_params.exists("INPUT_CLUSTERS")
          and (m_contigs->Load(m_params.getValue("INPUT_CLUSTERS").c_str())
              <= 0 or m_contigs->size() == 0)) {
        ERROR_MSG("Error reading input spectra from "
            << m_params.getValue("INPUT_CLUSTERS"));
        return false;
      }

      if (m_params.exists("INPUT_SPECS_PKLBIN")
          and (m_contigs->Load_pklbin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str())
              <= 0 or m_contigs->size() == 0)) {
        ERROR_MSG("Error reading input spectra from "
            << m_params.getValue("INPUT_SPECS_PKLBIN"));
        return false;
      }
    }

    if (m_params.exists("OUTPUT_MATCHED_SPECS_IDX")) {
      m_matchedSpectraIndices->resize(m_contigs->size());
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

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::saveOutputData(void)
  {
    if (m_matchedSpectra and m_params.exists("OUTPUT_MATCHED_SPECS")) {
        m_matchedSpectra->SaveSpecSet_pklbin(m_params.getValue("OUTPUT_MATCHED_SPECS").c_str());
    }
    
    m_matchedSpectra->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS").c_str(),
                                 m_params.getValue("OUTPUT_PSM").c_str(),
                                 m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX").c_str());

    //=================================================
    // LARS - TEMPORARY MEASURE SO REPORTS KEEP WORKING
    //=================================================
    if (m_params.exists("OUTPUT_MATCHED_PROTS")) {
      m_matchedSpectra->saveMatchedProts(m_params.getValue("OUTPUT_MATCHED_PROTS").c_str());
    }                   
    //=================================================

    if (m_matchedSpectraIndices and m_params.exists("OUTPUT_MATCHED_SPECS_IDX"))
      Save_binArray(m_params.getValue("OUTPUT_MATCHED_SPECS_IDX").c_str(),
                    *m_matchedSpectraIndices);

    if (m_matchedSpectraIndices and m_params.exists("OUTPUT_REF_SPS_NAMES")) {
	    ofstream spsIndex(m_params.getValue("OUTPUT_REF_SPS_NAMES").c_str(), ios_base::out);
	    for(unsigned int i = 0; i < m_matchedSpectraIndices->size(); i++) {
		    spsIndex << "sps:" << (*m_matchedSpectraIndices)[i] + 1 << endl;
		  }
	    spsIndex.close();
    }
    
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecSpecProtAlign::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TOLERANCE_PM");

    m_isValid = true;
    return true;
  }

// -------------------------------------------------------------------------


}

