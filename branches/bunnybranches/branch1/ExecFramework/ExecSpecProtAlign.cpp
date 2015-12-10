#include "ExecSpecProtAlign.h"

#include "alignment_modmut.h"
#include "AlignmentPenaltyBased.h"
#include "FdrPeptide.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"

// SpecNets Includes
#include "tags.h"
#include <limits.h>


#define DEBUG_SPECPROTALIGN 0

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
    m_inputSpectra(0), m_db(0x0), ownInput(true),
        m_matchedSpectraAll(0x0), m_matchedSpectra(0x0),
        m_matchedSpectraIndices(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------
  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       DB_fasta * db,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), m_db(db),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum), m_penaltyMatrixMods(penaltyMatrixMods), ownInput(true),
        m_matchedSpectraAll(0x0), m_matchedSpectra(0x0),
        m_matchedSpectraIndices(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams) :
    ExecBase(inputParams), m_inputSpectra(0x0), m_db(0x0), ownInput(true),
        m_matchedSpectraAll(0x0), m_matchedSpectra(0x0),
        m_matchedSpectraIndices(0x0), ownOutput(true)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";
  }

  // -------------------------------------------------------------------------

  ExecSpecProtAlign::ExecSpecProtAlign(const ParameterList & inputParams,
                                       SpecSet * inputSpectra,
                                       DB_fasta * db,
                                       PenaltyMatrix * penaltyMatrixBlosum,
                                       PenaltyMatrix * penaltyMatrixMods,
                                       SpecSet * matchedSpectraAll,
                                       SpecSet * matchedSpectra,
                                       vector<unsigned int> * matchedSpectraIndices) :
    ExecBase(inputParams), m_inputSpectra(inputSpectra), m_db(db),
        m_penaltyMatrixBlosum(penaltyMatrixBlosum),
        m_penaltyMatrixMods(penaltyMatrixMods), ownInput(false),
        m_matchedSpectraAll(matchedSpectraAll),
        m_matchedSpectra(matchedSpectra),
        m_matchedSpectraIndices(matchedSpectraIndices), ownOutput(false)
  {
    m_name = "ExecSpecProtAlign";
    m_type = "ExecSpecProtAlign";

    if (m_matchedSpectraIndices)
      m_matchedSpectraIndices->resize(m_inputSpectra->size());
  }


  // -------------------------------------------------------------------------

  ExecSpecProtAlign::~ExecSpecProtAlign(void)
  {
    if (ownInput) {
      if (m_inputSpectra)
        delete m_inputSpectra;
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
    if (!m_inputSpectra || m_inputSpectra->size() == 0) {
      ERROR_MSG("ERROR: empty set of input spectra");
      return false;
    }

   // More readable code if we use a reference instead of pointer
    SpecSet & inputSpectra = *m_inputSpectra;

    if (!m_db or m_db->size() == 0) {
      ERROR_MSG("ERROR: empty database");
      return false;
    }

    if (!m_matchedSpectraAll or !m_matchedSpectraIndices) {
      ERROR_MSG("ERROR: empty returned data pointer");
      return false;
    }

    float pmTol = (float) m_params.getValueDouble("TOLERANCE_PM");
    DEBUG_VAR(pmTol);
    float peakTol = (float) m_params.getValueDouble("TOLERANCE_PEAK");
    DEBUG_VAR(peakTol);

    unsigned int startIdx = m_params.exists("IDX_START") ?
            max(0, m_params.getValueInt("IDX_START")) :
            0;
    DEBUG_VAR(startIdx);
    unsigned int endIdx = m_params.exists("IDX_END") ?
            min(inputSpectra.size() - 1, (unsigned int) m_params.getValueInt("IDX_END")) :
            inputSpectra.size() - 1;
    DEBUG_VAR(endIdx);
    unsigned int startIdx_db = m_params.exists("IDX_START_DB") ?
            max(0, m_params.getValueInt("IDX_START_DB")) :
            0;
    unsigned int endIdx_db = m_params.exists("IDX_END_DB") ?
            min(m_db->size() - 1, (unsigned int) m_params.getValueInt("IDX_END_DB")) :
            m_db->size() - 1;
    if (startIdx > inputSpectra.size())
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
    int maxGapSize = m_params.exists("MAX_ALIGN_GAP_SIZE") ?
            max(0, m_params.getValueInt("MAX_ALIGN_GAP_SIZE")) : 8;
    int maxNumMods = m_params.exists("MAX_NUM_MODS") ?
            max(0, m_params.getValueInt("MAX_NUM_MODS")) : 1;
    int minNumMatchPeaks = m_params.exists("MIN_MATCHED_PEAKS_DB") ?
            max(0, m_params.getValueInt("MIN_MATCHED_PEAKS_DB")) : 7;
    int matchesPerMod = m_params.exists("MATCHES_PER_MOD") ?
            max(0, m_params.getValueInt("MATCHES_PER_MOD")) : 2;
    int specType = m_params.exists("SPEC_TYPE_MSMS") ?
            m_params.getValueInt("SPEC_TYPE_MSMS") > 0 : 0;
    bool tagParsimony =
            (bool) m_params.exists("MAX_PARSIMONY") ? ((int) m_params.getValueInt("MAX_PARSIMONY") ? 1 : 0) : 0;
    DEBUG_VAR(tagParsimony);

    bool penaltyAlign = (bool) m_params.exists("PENALTY_ALIGNMENT") ?
    		((int) m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;
    DEBUG_VAR(penaltyAlign);
    float penaltyAlpha = m_params.exists("PENALTY_ALIGNMENT_ALPHA") ?
            m_params.getValueFloat("PENALTY_ALIGNMENT_ALPHA") : 1.0;
    DEBUG_VAR(penaltyAlpha);
    float penaltyBeta = m_params.exists("PENALTY_ALIGNMENT_BETA") ?
            m_params.getValueFloat("PENALTY_ALIGNMENT_BETA") : 1.0;
    DEBUG_VAR(penaltyBeta);

    float ionOffset = specType ? AAJumps::massHion : 0;
    DEBUG_VAR(ionOffset);
    int modIdx, specIdx, protIdx;

    Spectrum tmpSpec;
    tmpSpec.reserve(1024);
    Spectrum cSpec;
    cSpec.reserve(1024);
    Spectrum cSpecRev;
    cSpecRev.reserve(1024);
    AMM_match *bestMatch;

    DEBUG_VAR(inputSpectra.size());

    // If no spectrum contain PSMs, then match to all proteins
    bool matchAll = true;
    for (specIdx = startIdx; specIdx <= endIdx; specIdx++) {
      if (inputSpectra[specIdx].psmList.size() != 0) {
        matchAll = false;
        break;
      }
    }

    for (specIdx = startIdx; specIdx <= endIdx; specIdx++) {

      DEBUG_VAR(specIdx);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[specIdx].parentMass);
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[specIdx].size());

      // Apparently it is possible to have a contig that is empty
      if (inputSpectra[specIdx].size() == 0) {
        m_matchedSpectraAll->push_back(inputSpectra[specIdx]);
        continue;
      }

      cSpec = inputSpectra[specIdx];
      cSpec.psmList.clear();

      // Get the reverse spectrum
      cSpec.reverse(0, &cSpecRev);
      cSpecRev.psmList.clear();

#if 0
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.size());
      if (DEBUG_SPECPROTALIGN) DEBUG_MSG("Original spectrum");
      if (DEBUG_SPECPROTALIGN) cSpec.output(cerr);
      // Clean up any pre-existing endpoints so we don't add too many
      list<int> peaksToRemove;
      for (int iPeak = 0; iPeak < cSpec.size(); iPeak++) {
        if (cSpec[iPeak][0] < 57.0 || cSpec[iPeak][0] > cSpec.parentMass - 57.0) {
          peaksToRemove.push_back(iPeak);
        }
      }
      if (peaksToRemove.size() != 0) {
        cSpec.removePeaks(peaksToRemove);
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.size());

      if (DEBUG_SPECPROTALIGN) DEBUG_MSG("Cleaned spectrum");
      if (DEBUG_SPECPROTALIGN) cSpec.output(cerr);
      // Add end peaks if necessary
      if (enforceEndpeaks) {
        cSpec.addZPMpeaks(peakTol, ionOffset, false);
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_MSG("Cleaned spectrum with peaks");
      if (DEBUG_SPECPROTALIGN) cSpec.output(cerr);

      if (DEBUG_SPECPROTALIGN) DEBUG_MSG("Reversed spectrum");
      if (DEBUG_SPECPROTALIGN) cSpecRev.output(cerr);
      // Add end peaks if necessary
      if (enforceEndpeaks) {
        cSpecRev.addZPMpeaks(peakTol, ionOffset, false);
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_MSG("Reversed spectrum with peaks");
      if (DEBUG_SPECPROTALIGN) cSpecRev.output(cerr);
      if (DEBUG_SPECPROTALIGN) DEBUG_TRACE;
#endif

      map<int, set<float> > matchedProtMap;
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[specIdx].psmList.size())
      if (inputSpectra[specIdx].psmList.size() != 0) {
        // Coalesce the starting positions on all proteins
        list<psmPtr>::iterator itr = inputSpectra[specIdx].psmList.begin();
        list<psmPtr>::iterator itr_end = inputSpectra[specIdx].psmList.end();
        for (; itr != itr_end; itr++) {
          int idx = (*itr)->m_dbIndex;
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_protein)
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_annotation)
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_dbIndex)
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR(tagParsimony)
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_matchOrientation)
          // Only use the forward orientation tags(if tag parsimony was not done)
          if (tagParsimony || (*itr)->m_matchOrientation == 0) {
            if ( matchedProtMap.find(idx) == matchedProtMap.end()) {
              set<float> newList;
              matchedProtMap[idx] = newList;
            }
            matchedProtMap[idx].insert((*itr)->m_startMass);
            if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_startMass)
          }
        } // for (; itr != itrEnd; itr++) {
      } else if (matchAll) {
    	 DEBUG_MSG("Inserting dummy start mass for every DB index")
        // If there are no PSM's insert dummy start mass for every db idx
        // This will cause alignment to every protein in database
        for (int idx = 0; idx < m_db->size(); idx++) {
          matchedProtMap[idx].insert(0.0);
        }
      }

      DEBUG_MSG("Matching as b...");

      //itrMap = matchedProtMap.begin();
      map<int, set<float> >::iterator itrMap = matchedProtMap.begin();
      map<int, set<float> >::iterator itrMapEnd = matchedProtMap.end();
      for (; itrMap !=  itrMapEnd; itrMap++) {
        int protIdx = itrMap->first;
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(protIdx);
        // If parsimony was done in TagSearch then don't use start positions
        if (tagParsimony) {
          itrMap->second.clear();
        }

        if (!penaltyAlign) {
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
        } else {
          alignmentPenaltyBased(cSpec,
                                m_db->getMassesSpec(protIdx),
                                m_db->getSequence(protIdx),
                                protIdx,
                                0,
                                itrMap->second,
                                m_penaltyMatrixMods,
                                m_penaltyMatrixBlosum,
                                penaltyAlpha,
                                penaltyBeta,
                                maxNumMods,
                                minNumMatchPeaks,
                                maxGapSize,
                                pmTol,
                                peakTol,
                                57,
                                maxModMass,
                                minModMass,
                                enforceEndpeaks);
        }

      } // for (; itrMap !=  itrMapEnd; itrMap++) {

      DEBUG_MSG("Matching as y...");
      //cSpecRev.output(cerr);

      matchedProtMap.clear();
      // Coalesce the starting positions on all proteins
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(inputSpectra[specIdx].psmList.size())
      if (inputSpectra[specIdx].psmList.size() != 0) {
        list<psmPtr>::iterator itr = inputSpectra[specIdx].psmList.begin();
        list<psmPtr>::iterator itr_end = inputSpectra[specIdx].psmList.end();
        for (; itr != itr_end; itr++) {
          int idx = (*itr)->m_dbIndex;
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_protein)
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_annotation)
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_dbIndex)
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR(tagParsimony)
          if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_matchOrientation)
          // Only use the reverse orientation tags (if tag parsimony was not done)
          if (tagParsimony || (*itr)->m_matchOrientation == 1) {
            if ( matchedProtMap.find(idx) == matchedProtMap.end()) {
              set<float> newList;
              matchedProtMap[idx] = newList;
            }
            matchedProtMap[idx].insert((*itr)->m_startMass);
            if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*itr)->m_startMass)
          }
        } // for (; itr != itrEnd; itr++) {
      } else if (matchAll) {
         // If there are no PSM's insert dummy start mass for every db idx
         // This will cause alignment to every protein in database
        for (int idx = 0; idx < m_db->size(); idx++) {
          matchedProtMap[idx].insert(0.0);
        }
      }

      itrMap = matchedProtMap.begin();
      for (; itrMap !=  itrMapEnd; itrMap++) {
        int protIdx = itrMap->first;
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR(protIdx);

        // If parsimony was done in TagSearch then don't use start positions
        if (tagParsimony) {
          itrMap->second.clear();
        }

        if (!penaltyAlign) {
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
        } else {
          alignmentPenaltyBased(cSpecRev,
                                m_db->getMassesSpec(protIdx),
                                m_db->getSequence(protIdx),
                                protIdx,
                                1,
                                itrMap->second,
                                m_penaltyMatrixMods,
                                m_penaltyMatrixBlosum,
                                penaltyAlpha,
                                penaltyBeta,
                                maxNumMods,
                                minNumMatchPeaks,
                                maxGapSize,
                                pmTol,
                                peakTol,
                                57,
                                maxModMass,
                                minModMass,
                                enforceEndpeaks);
        }

      } // for (; itr != itrEnd; itr++) {

      DEBUG_MSG("Matching Complete");

      // Find the top score
      float bestScore = -(float)INT_MAX;
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpec.psmList.size());
      list<psmPtr>::iterator litr = cSpec.psmList.begin();
      list<psmPtr>::iterator litrEnd = cSpec.psmList.end();
      for (; litr != litrEnd; litr++) {
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*litr)->m_scanNum);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*litr)->m_origAnnotation);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*litr)->m_annotation);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*litr)->m_score);
        if ((*litr)->m_score > bestScore)
          bestScore = (*litr)->m_score;
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(bestScore);

      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(cSpecRev.psmList.size());
      litr = cSpecRev.psmList.begin();
      litrEnd = cSpecRev.psmList.end();
      for (; litr != litrEnd; litr++) {
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*litr)->m_scanNum);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*litr)->m_origAnnotation);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*litr)->m_annotation);
        if (DEBUG_SPECPROTALIGN) DEBUG_VAR((*litr)->m_score);
        if ((*litr)->m_score > bestScore)
          bestScore = (*litr)->m_score;
      }
      if (DEBUG_SPECPROTALIGN) DEBUG_VAR(bestScore);

      cSpec.psmList.remove_if(LessThanPredicate(bestScore));
      cSpecRev.psmList.remove_if(LessThanPredicate(bestScore));
      if (cSpecRev.psmList.size() == 0) {
        m_matchedSpectraAll->push_back(cSpec);
      } else {
        m_matchedSpectraAll->push_back(cSpecRev);
        //DEBUG_VAR(m_matchedSpectraAll->back().size());
      }

    } // for (specIdx = startIdx; specIdx <= endIdx; specIdx++)

    DEBUG_VAR(m_matchedSpectraAll->size());

    // We could have multiple hits on one protein
    // So we perform a sort() and unique() to eliminate them
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      (*m_matchedSpectraAll)[i].psmList.sort(PsmDbIndexSort);
      (*m_matchedSpectraAll)[i].psmList.unique(PsmDbIndexUnique);
    }

    // Set backwards pointers/scans for PSM--->Spectrum link
    list<psmPtr>::iterator iter;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      (*m_matchedSpectraAll)[i].scan = i + 1;
      for (iter = (*m_matchedSpectraAll)[i].psmList.begin();
           iter != (*m_matchedSpectraAll)[i].psmList.end();
           iter++) {
        (*iter)->m_spectrum = &(*m_matchedSpectraAll)[i];
        (*iter)->m_scanNum = (*m_matchedSpectraAll)[i].scan;
        (*iter)->m_protein = m_db->getID((*iter)->m_dbIndex);

        if (tagParsimony) {
          string annotation;
          (*iter)->getAnnotationFromMatchedPeaks(m_db->getMassesSpec((*iter)->m_dbIndex),
                                                 m_db->getSequence((*iter)->m_dbIndex),
                                                 annotation);
          (*iter)->m_annotation = (*iter)->m_origAnnotation;
          (*iter)->m_annotation = annotation;
          DEBUG_VAR((*iter)->m_annotation);
        }

        // For now copy annotation to original annotation(we'll fix this later)
        (*iter)->m_origAnnotation = (*iter)->m_annotation;
      }
    }

    if (m_params.getValueFloat("FDR_THRESHOLD", -1.0) >= 0.0) {

      vector<pair<unsigned int,bool> > ntermMods;
      ntermMods.resize(2);
      ntermMods[0].first = 42;   ntermMods[0].second = false;  // N-term acetylation
      ntermMods[1].first = 111;  ntermMods[1].second = true;   // N-term PyroGlu

      // Turn known mods from penalty matrix into known mods for spectral probabilities
      vector<unsigned int> mods;
      if (m_penaltyMatrixMods != 0) {
        const map<string, set<float> > & knownMods = m_penaltyMatrixMods->getKnownMods();
        map<string, set<float> >::const_iterator itr = knownMods.begin();
        map<string, set<float> >::const_iterator itrEnd = knownMods.end();
        for ( ; itr != itrEnd; itr++) {
          string stringAA = itr->first;
          //DEBUG_VAR(stringAA);
          float aaMass = m_penaltyMatrixMods->getMass(itr->first);
          //DEBUG_VAR(aaMass);
          if (aaMass == 0) continue; // Sanity check
          const set<float> & setMods = itr->second;
          set<float>::const_iterator itrSet = setMods.begin();
          set<float>::const_iterator itrSetEnd = setMods.end();
          for ( ; itrSet != itrSetEnd; itrSet++) {
            float mod = *itrSet;
            //DEBUG_VAR(mod);
            float totalMass = aaMass + mod;
            //DEBUG_VAR(totalMass);
            mods.push_back(totalMass);
          }
        }
        DEBUG_VAR(mods.size());
      }

      for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
        (*m_matchedSpectraAll)[i].computeSpectralProbabilities(ntermMods, mods, peakTol);
      }
    }

    DEBUG_VAR(m_matchedSpectraAll->size());

    if (m_params.exists("GRID_CLONE_INTERNAL_FLAG")) {
      DEBUG_MSG("Exiting early for grid execution");
      return true;
    }

    CompressResults();

    return true;
  }

  // -------------------------------------------------------------------------

  void ExecSpecProtAlign::CompressResults(void)
  {
    DEBUG_VAR(m_matchedSpectraAll->size());

    if (m_params.getValueFloat("FDR_THRESHOLD", -1.0) >= 0.0) {

      PeptideSpectrumMatchSet psmSetTemp;
      PeptideSpectrumMatchSet psmSetTempOut;
      psmSetTemp.getPSMSet(m_matchedSpectraAll);
      psmSetTemp.saveToFile("psm1.txt");
      DEBUG_VAR(psmSetTemp.size());
      if (!FdrPeptide::concatenatedTargetDecoy(psmSetTemp, psmSetTempOut, 1.0)) {
        ERROR_MSG("FDR sorting failed.");
      }
      psmSetTempOut.saveToFile("psm2.txt");

      double fdrCutoff  = 1.0;
      if (m_params.getValueFloat("FDR_THRESHOLD", -1.0) >= 0.0) {
        fdrCutoff = m_params.getValueFloat("FDR_THRESHOLD", 1.0);  // Purposely 1.0 not -1.0 (should exist anyway)
      }
      DEBUG_VAR(fdrCutoff);

      if (!FdrPeptide::filterByPValue(psmSetTempOut, fdrCutoff)) {
        ERROR_MSG("Filter by FDR value failed.");
      }
      psmSetTempOut.saveToFile("psm3.txt");
      DEBUG_VAR(psmSetTempOut.size());
      m_matchedSpectraAll->clearPsms();
      DEBUG_TRACE;
      psmSetTempOut.addSpectra(m_matchedSpectraAll);
    }

    DEBUG_TRACE;
    m_matchedSpectraAll->maximumParsimony();

    // Resize to make enough room
    m_matchedSpectraIndices->resize(m_matchedSpectraAll->size());
    m_matchedSpectra->resize(m_matchedSpectraAll->size());

    int keepIdx = 0;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      //DEBUG_VAR(i);
      //DEBUG_VAR((*m_matchedSpectraAll)[i].size());
      //DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
      if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
        (*m_matchedSpectraIndices)[keepIdx] = i;
        (*m_matchedSpectra)[keepIdx++] = (*m_matchedSpectraAll)[i];
      }
    }

    DEBUG_VAR(keepIdx);
    m_matchedSpectra->resize(keepIdx);
    m_matchedSpectraIndices->resize(keepIdx);

    DEBUG_VAR(m_matchedSpectra->size());
    DEBUG_VAR(m_matchedSpectraIndices->size());

    return;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::loadInputData(void)
  {
    if (ownInput) {
      if (!m_inputSpectra)
        m_inputSpectra = new SpecSet;
      if (!m_db)
        m_db = new DB_fasta;
    }
    m_inputSpectra->resize(0);

    if (ownOutput) {
      if (!m_matchedSpectraAll)
        m_matchedSpectraAll = new SpecSet;
      if (!m_matchedSpectra)
        m_matchedSpectra = new SpecSet;
      if (!m_matchedSpectraIndices)
        m_matchedSpectraIndices = new vector<unsigned int> ;
    }
    m_matchedSpectraAll->resize(0);
    m_matchedSpectra->resize(0);
    m_matchedSpectraIndices->resize(0);

    if (m_params.exists("INPUT_SPECS_PKLBIN") and m_params.exists("INPUT_PSM")) {
      if (m_inputSpectra->loadPklBin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str(),
                                     m_params.getValue("INPUT_PSM").c_str(),
                                     m_params.getValue("INPUT_MATCHED_PEAKS_IDX").c_str()) < 0) {
        ERROR_MSG("Error reading input spectra files.");
        return false;
      }
    }

    if (m_params.exists("OUTPUT_MATCHED_SPECS_IDX")) {
      m_matchedSpectraIndices->resize(m_inputSpectra->size());
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

  bool ExecSpecProtAlign::saveOutputData(void)
  {
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
    //=================================================

    bool penaltyAlign = (bool) m_params.exists("PENALTY_ALIGNMENT") ?
  	  	((int) m_params.getValueInt("PENALTY_ALIGNMENT") ? 1 : 0) : 0;
    DEBUG_VAR(penaltyAlign);

    if (m_params.exists("OUTPUT_MATCHED_PEAKS_IDX_ALL")) {

      AAJumps aaJumps(1);
      if (m_params.exists("AMINO_ACID_MASSES")) {
        DEBUG_MSG("Loading amino acid masses...");
        if (!aaJumps.loadJumps(m_params.getValue("AMINO_ACID_MASSES").c_str(), false)) {
          ERROR_MSG("Error reading input amino acid mass file " << m_params.getValue("AMINO_ACID_MASSES"));
          return false;
        }
      }

      DEBUG_MSG("Saving matched peaks all...");
      SpecSet tempMatchedPeaks;
      tempMatchedPeaks.resize(m_matchedSpectraAll->size());
      for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
        DEBUG_VAR(i);
        DEBUG_VAR((*m_matchedSpectraAll)[i].psmList.size());
        if ((*m_matchedSpectraAll)[i].psmList.size() != 0) {
          list<psmPtr>::iterator litr = (*m_matchedSpectraAll)[i].psmList.begin();

          DEBUG_VAR((*litr)->m_annotation);
          //if (penaltyAlign) DEBUG_VAR((*litr)->m_origAnnotation);
          (*litr)->getMatchedPeaksFromAnnotation(m_db->getMassesSpec((*litr)->m_dbIndex),
                                                 aaJumps);

          int peakListSize = (*litr)->m_matchedPeaks.size();
          DEBUG_VAR(peakListSize);
          tempMatchedPeaks[i].resize(peakListSize);
          for (int j = 0; j < peakListSize; j++) {
            tempMatchedPeaks[i][j].set((*litr)->m_matchedPeaks[j][0],
                                       (*litr)->m_matchedPeaks[j][1]);
            DEBUG_MSG((*litr)->m_matchedPeaks[j][0] << "  " << (*litr)->m_matchedPeaks[j][1]);
          }
        }
      }
      tempMatchedPeaks.SaveSpecSet_pklbin(m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str());
    }

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

    DEBUG_VAR(m_matchedSpectraAll->size());
    if (m_matchedSpectraAll and m_params.exists("OUTPUT_MATCHED_SPECS_ALL")) {
      if (m_params.exists("OUTPUT_PSM_ALL") && m_params.exists("OUTPUT_MATCHED_PEAKS_IDX_ALL")) {
        DEBUG_MSG("Saving matched specs all (3 files)...");
        m_matchedSpectraAll->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str(),
                                        m_params.getValue("OUTPUT_PSM_ALL").c_str(),
                                        m_params.getValue("OUTPUT_MATCHED_PEAKS_IDX_ALL").c_str());

      } else if (m_params.exists("OUTPUT_PSM_ALL")) {
        DEBUG_MSG("Saving matched specs all (2 files)...");
        m_matchedSpectraAll->savePklBin(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str(),
                                        m_params.getValue("OUTPUT_PSM_ALL").c_str());

      } else {
        DEBUG_MSG("Saving matched specs all (1 file)...");
        m_matchedSpectraAll->SaveSpecSet_pklbin(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str());
      }
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

    DEBUG_TRACE;
    if (m_matchedSpectraIndices and m_params.exists("OUTPUT_MATCHED_SPECS_IDX"))
      Save_binArray(m_params.getValue("OUTPUT_MATCHED_SPECS_IDX").c_str(),
                    *m_matchedSpectraIndices);

    if (m_params.exists("OUTPUT_PSM_MOD_FREQS")) {
      PeptideSpectrumMatchSet psmSetTemp;
      psmSetTemp.getPSMSet(m_matchedSpectra);
      psmSetTemp.saveModMatrix(m_params.getValue("OUTPUT_PSM_MOD_FREQS").c_str(), penaltyAlign);
    }

    DEBUG_VAR(m_matchedSpectraIndices->size());
    if (m_matchedSpectraIndices and m_params.exists("OUTPUT_REF_SPS_NAMES")) {
	    ofstream spsIndex(m_params.getValue("OUTPUT_REF_SPS_NAMES").c_str(), ios_base::out | ios_base::binary);
	    for(unsigned int i = 0; i < m_matchedSpectraIndices->size(); i++) {
		    spsIndex << "sps:" << (*m_matchedSpectraIndices)[i] + 1 << endl;
		  }
	    spsIndex.close();
    }

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::saveInputData(std::vector<std::string> & filenames)
  {
    string dataDir = m_params.getValue("GRID_DATA_DIR");
    if (dataDir.empty())
    {
      dataDir = ".";
    }
    string baseDirectory = dataDir + "/";
    string baseFilename = baseDirectory + getName();

    string specFilename = baseDirectory + "sps_seqs.pklbin";
    string psmFilename = baseDirectory + "tagsearchpsm.txt";
//    if (!fileExists(specFilename))
    {
      DEBUG_MSG("Saving " << specFilename);
      m_inputSpectra->savePklBin(specFilename.c_str(),
                                 psmFilename.c_str());
    }
//    else
//    {
//      DEBUG_MSG("Not Saving " << specFilename << " (already exists)");
//    }

    DEBUG_VAR(m_db);
    string dbFilename = baseDirectory + "db.fasta";
//    if (!fileExists(dbFilename))
    {
      DEBUG_MSG("Saving " << dbFilename );
      m_db->Save(dbFilename.c_str());
    }
//    else
//    {
//      DEBUG_MSG("Not Saving " << dbFilename  << " (already exists)");
//    }

    string aaFilename = baseDirectory + "aminoacids.txt";
//    if (!fileExists(aaFilename))
    {
      DEBUG_MSG("Saving " << aaFilename);
      m_penaltyMatrixMods->saveAminoAcids(aaFilename);
    }
//    else
//    {
//      DEBUG_MSG("Not Saving " << aaFilename << " (already exists)");
//    }

    string knownmodsFilename = baseDirectory + "knownmods.txt";
//    if (!fileExists(knownmodsFilename ))
    {
      DEBUG_MSG("Saving " << knownmodsFilename);
      m_penaltyMatrixMods->saveKnownMods(knownmodsFilename );
    }
//    else
//    {
//      DEBUG_MSG("Not Saving " << knownmodsFilename << " (already exists)");
//    }

    string modPenaltiesFilename = baseDirectory + "modpenalties.txt";
//    if (!fileExists(modPenaltiesFilename))
    {
      DEBUG_MSG("Saving " << modPenaltiesFilename);
      m_penaltyMatrixMods->saveMatrix(modPenaltiesFilename);
    }
//    else
//    {
//      DEBUG_MSG("Not Saving " << modPenaltiesFilename << " (already exists)");
//    }

    string blossumFilename = baseDirectory + "blosumpenalties.txt";
//    if (!fileExists(blossumFilename))
    {
      DEBUG_MSG("Saving " << blossumFilename);
      m_penaltyMatrixBlosum->saveMatrix(blossumFilename);
    }
//    else
//    {
//      DEBUG_MSG("Not Saving " << blossumFilename << " (already exists)");
//    }

    m_params.setValue("INPUT_SPECS_PKLBIN", specFilename);
    m_params.setValue("INPUT_PSM", psmFilename);
    m_params.setValue("INPUT_FASTA", dbFilename);

    m_params.setValue("AMINO_ACID_MASSES", aaFilename);
    m_params.setValue("KNOWN_MODS_FILE", knownmodsFilename);
    m_params.setValue("MODS_PENALTY_FILE", modPenaltiesFilename);
    m_params.setValue("BLOSUM_PENALTY_FILE", blossumFilename);

    string paramFilename = baseFilename + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(specFilename);
    filenames.push_back(psmFilename);
    filenames.push_back(dbFilename);
    filenames.push_back(modPenaltiesFilename);
    filenames.push_back(aaFilename);
    filenames.push_back(knownmodsFilename );
    filenames.push_back(blossumFilename);

    return true;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::loadOutputData(void)
  {
    if (m_matchedSpectraAll == 0x0) {
      ownOutput = true;
      m_matchedSpectraAll = new SpecSet;
      m_matchedSpectra = new SpecSet;
      m_matchedSpectraIndices = new vector<unsigned int> ;
    }

    if (m_params.exists("OUTPUT_MATCHED_SPECS_ALL")) {
      m_matchedSpectraAll->loadPklBin(m_params.getValue("OUTPUT_MATCHED_SPECS_ALL").c_str(),
                                      m_params.getValue("OUTPUT_PSM_ALL").c_str() );

    }

    DEBUG_TRACE;
    return true;
  }

  // -------------------------------------------------------------------------

  vector<ExecBase *> const & ExecSpecProtAlign::split(int numSplit)
  {
    DEBUG_VAR(numSplit);

    if (numSplit < 2) {
      DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
      return m_subModules;
    }

    int spectraSize = m_inputSpectra->size();
    DEBUG_VAR(spectraSize);
    if (spectraSize == 0) {
      DEBUG_MSG("Must have at least one spectra");
      return m_subModules;
    }

    int startBaseIdx;
    if (m_params.exists("IDX_START")) {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    } else {
      startBaseIdx = 0;
    }
    DEBUG_VAR(startBaseIdx);
    int endBaseIdx;
    if (m_params.exists("IDX_END")) {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
      endBaseIdx = min(endBaseIdx , (int)spectraSize - 1);
    } else {
      endBaseIdx = spectraSize - 1;
    }
    DEBUG_VAR(endBaseIdx);

    int numSpectraPerSplit = (endBaseIdx - startBaseIdx) / numSplit;
    int modSpectraPerSplit = (endBaseIdx - startBaseIdx) % numSplit;
    if (modSpectraPerSplit != 0) {
      numSpectraPerSplit++;
    }

    int indexStart = startBaseIdx;
    int indexEnd = startBaseIdx + numSpectraPerSplit - 1;

    for (int i = 0; i < numSplit; i++) {

      if (startBaseIdx >= spectraSize) {
        break;
      }

      //DEBUG_VAR(i);

      // Copy the parameters
      ParameterList childParams(m_params);
      // Set the start and end indices
      char buf[128];
      sprintf(buf, "%d", indexStart);
      childParams.setValue("IDX_START", buf);
      sprintf(buf, "%d", indexEnd);
      childParams.setValue("IDX_END", buf);

      DEBUG_MSG("Start [" << indexStart << "] End [" << indexEnd << "] Split [" << i << "]");

      // Make a clone of this module
      ExecBase * theClone = new ExecSpecProtAlign(childParams,
                                                  m_inputSpectra,
                                                  m_db,
                                                  m_penaltyMatrixBlosum,
                                                  m_penaltyMatrixMods);

      // Give it a new name based on the split
      theClone->setName(makeName(m_name, i));

      // Have to set up the output files also so the params will be correct on reload
      string dataDir = m_params.getValue("GRID_DATA_DIR");
      if (dataDir.empty()) {
        dataDir = ".";
      }
      string baseName = dataDir + "/" + theClone->getName();
      //DEBUG_VAR(baseName);

      theClone->m_params.setValue("OUTPUT_MATCHED_SPECS_ALL",     baseName + "_contigs_all.pklbin");
      theClone->m_params.setValue("OUTPUT_PSM_ALL",               baseName + "_psm_all.txt");
      theClone->m_params.setValue("GRID_CLONE_INTERNAL_FLAG", "yes");

      theClone->m_params.removeParam("OUTPUT_MATCHED_PEAKS_IDX");
      theClone->m_params.removeParam("OUTPUT_MATCHED_PEAKS_IDX_ALL");
      theClone->m_params.removeParam("OUTPUT_MATCHED_SPECS");
      theClone->m_params.removeParam("OUTPUT_REF_SPS_NAMES");

      std::string suffix("");
      char bufSplit[128];
      sprintf(bufSplit, "%d", i + 1);
      theClone->m_params.setValue("NUM_SPLIT", bufSplit);

      m_subModules.push_back(theClone);

      indexStart += numSpectraPerSplit;
      indexEnd += numSpectraPerSplit;

    } // for (int i = 0; i < numSplit; i++)

    DEBUG_VAR(m_subModules.size());

    return m_subModules;
  }

  // -------------------------------------------------------------------------

  bool ExecSpecProtAlign::merge(void)
  {
    if (m_matchedSpectraAll == 0x0)
    {
      ownOutput = true;
      m_matchedSpectraAll = new SpecSet;
      m_matchedSpectra = new SpecSet;
      m_matchedSpectraIndices = new vector<unsigned int> ;
    }

    DEBUG_VAR(m_subModules.size());
    int iSpec = 0;
    for (int i = 0; i < m_subModules.size(); i++)
    {
      ExecSpecProtAlign * espa = (ExecSpecProtAlign*)m_subModules[i];
      //DEBUG_VAR(espa);
      //DEBUG_VAR(espa->m_matchedSpectraAll->size());
      for (int j = 0; j < espa->m_matchedSpectraAll->size(); j++) {
        //DEBUG_VAR(iSpec);
        m_matchedSpectraAll->push_back(espa->m_matchedSpectraAll->operator[](j));
        //DEBUG_VAR((*m_matchedSpectraAll)[iSpec].size());
        //DEBUG_VAR((*m_matchedSpectraAll)[iSpec].psmList.size());
        iSpec++;
      }
    }

    // Set backwards pointers/scans for PSM--->Spectrum link
    list<psmPtr>::iterator iter;
    for (int i = 0; i < m_matchedSpectraAll->size(); i++) {
      (*m_matchedSpectraAll)[i].scan = i + 1;
      for (iter = (*m_matchedSpectraAll)[i].psmList.begin();
           iter != (*m_matchedSpectraAll)[i].psmList.end();
           iter++) {
        (*iter)->m_spectrum = &(*m_matchedSpectraAll)[i];
        (*iter)->m_scanNum = (*m_matchedSpectraAll)[i].scan;
        //DEBUG_VAR((*iter)->m_origAnnotation);
        //DEBUG_VAR((*iter)->m_annotation);
      }
    }


    CompressResults();

    return true;
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

