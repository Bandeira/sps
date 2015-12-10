// Header Includes
#include "ExecPrmScoring.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "ParameterList.h"
#include "ExecMergeConvert.h"
#include "ExecFilterPairs.h"
#include "DeconvSpectrum.h"
#include "PairedSpecSet.h"
#include "UnionFind.h"
#include "prm_alignment.h"
#include "mzrange.h"

// External Includes
#include "Specific.h"

// System Includes
#include <stdlib.h>
#include <fstream>

using namespace std;
using namespace specnets;
using namespace spsReports;

namespace specnets
{

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(void) :
    ownInput(true), m_inputSpectra(0x0), ownOutput(true), m_outputSpectra(0x0),
        m_outputClusters(0x0), enforceDaTolerance(false), m_loader(0x0),
        m_loadedIndices(0x0), m_specIDs(0x0), m_fasta(0x0)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
    m_scanToIdx.resize(0);
    m_scanToInputScan.resize(0);
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(const ParameterList & inputParams) :
    ExecBase(inputParams), ownInput(true), m_inputSpectra(0x0),
        ownOutput(true), m_outputSpectra(0x0), m_outputClusters(0x0),
        enforceDaTolerance(false), m_loader(0x0), m_loadedIndices(0x0),
        m_specIDs(0x0), m_fasta(0x0)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
    m_scanToIdx.resize(0);
    m_scanToInputScan.resize(0);
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(const ParameterList & inputParams,
                                 SpecSet * inputSpectra,
                                 SpecSet * outputSpectra,
                                 ClusterData * outputClusters) :
    ExecBase(inputParams), ownInput(false), m_inputSpectra(inputSpectra),
        ownOutput(false), m_outputSpectra(outputSpectra),
        m_outputClusters(outputClusters), enforceDaTolerance(false),
        m_loader(0x0), m_loadedIndices(0x0), m_specIDs(0x0), m_fasta(0x0)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
    m_scanToIdx.resize(0);
    m_scanToInputScan.resize(0);
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::ExecPrmScoring(const ParameterList & inputParams,
                                 SpecSet * outputSpectra,
                                 ClusterData * outputClusters) :
    ExecBase(inputParams), ownInput(true), m_inputSpectra(0x0),
        ownOutput(false), m_outputSpectra(outputSpectra),
        m_outputClusters(outputClusters), enforceDaTolerance(false),
        m_loader(0x0), m_loadedIndices(0x0), m_specIDs(0x0), m_fasta(0x0)
  {
    m_name = "ExecPrmScoring";
    m_type = "ExecPrmScoring";
    m_scanToIdx.resize(0);
    m_scanToInputScan.resize(0);
  }

  // -------------------------------------------------------------------------
  ExecPrmScoring::~ExecPrmScoring(void)
  {
    if (ownInput && ownOutput && m_inputSpectra == m_outputSpectra
        && m_outputSpectra != 0)
    {
      delete m_outputSpectra;
      m_inputSpectra = 0;
      m_outputSpectra = 0;
    }
    else if (m_inputSpectra == m_outputSpectra && m_outputSpectra != 0)
    {
      m_inputSpectra = 0;
      m_outputSpectra = 0;
    }

    if (ownInput && m_inputSpectra != 0)
    {
      delete m_inputSpectra;
    }
    if (ownOutput && m_outputSpectra != 0)
    {
      delete m_outputSpectra;
    }
    if (ownOutput && m_outputClusters != 0)
    {
      delete m_outputClusters;
    }

    if (m_loader != 0)
    {
      delete m_loader;
    }

    if (m_loadedIndices != 0)
    {
      delete m_loadedIndices;
    }
    if (m_fasta != 0)
    {
      delete m_fasta;
    }
    if (m_specIDs != 0)
    {
      delete m_specIDs;
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecPrmScoring::clone(const ParameterList & inputParams) const
  {
    return new ExecPrmScoring(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::invoke(void)
  {
    if (m_inputSpectra == 0)
    {
      ERROR_MSG("No input spectra!");
      return false;
    }

    if (ownOutput && m_outputSpectra == 0)
    { // initialize output spectra if needed
      m_outputSpectra = new SpecSet;
    }

    if (ownOutput && m_outputClusters == 0)
    {
      m_outputClusters = new ClusterData;
    }

    if (m_loader == 0)
    { // initialize converter/loader
      m_loadedIndices = new vector<pair<int, int> > ;
      m_loader = new ExecMergeConvert(m_params,
                                      m_loadedIndices,
                                      m_inputSpectra,
                                      m_outputSpectra);
    }

    // Do any pre-processing here
    if (!m_loader->invoke())
    {
      return false;
    }

    delete m_loader;
    m_loader = 0;

    m_scanToIdx.resize(m_outputSpectra->size());
    m_scanToInputScan.resize(m_outputSpectra->size());
    for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
    {
      m_scanToInputScan[i] = (*m_outputSpectra)[i].scan;
      (*m_outputSpectra)[i].scan = i;
      m_scanToIdx[i].first = 0;
      m_scanToIdx[i].second = i;
    }

    // if true, just enforce Da tolerance rather than capturing PPM tolerance from input MS/MS spectra
    enforceDaTolerance = m_params.getValueInt("ENFORCE_DA_TOL", 1) > 0;

    // if true, do clustering of PRM spectra
    bool mergeSamePrec = m_params.getValueInt("MERGE_SAME_PREC", 0) > 0;

    // number of consecutive spectra from the same precursor
    int numConsec = m_params.getValueInt("NUM_CONSECUTIVE", 1);

    vector<vector<bool> >* reversedPeaks = (vector<vector<bool> >*)0;

    if (mergeSamePrec && !enforceDaTolerance)
    {
      // details which peaks were derived from suffix fragments and thus need their tolerances adjusted
      reversedPeaks = new vector<vector<bool> > (m_outputSpectra->size());
    }

    vector<vector<unsigned int> > pairedScans;
    m_outputSpectra->extractPairedScans(pairedScans, numConsec);

    DEBUG_VAR(pairedScans.size());

    SpecSet prmSpecs;
    SpecSet ms2Specs;

    if (m_params.exists("PEPNOVO_MODEL")
        && m_params.exists("PEPNOVO_INPUT_MGF")
        && m_params.exists("PEPNOVO_OUTPUT_PRMS"))
    {
      try
      {
        bool res = invokePepNovo(m_params,
                                 *m_outputSpectra,
                                 prmSpecs,
                                 reversedPeaks);
        if (!res)
        {
          throw 20;
        }
      }
      catch (exception& e)
      {
        ERROR_MSG("Exception occurred: " << e.what());
        if (reversedPeaks)
        {
          delete reversedPeaks;
        }
        return false;
      }
      DEBUG_TRACE;
    }
    else
    {

      // extract spectra for each fragmentation mode
      SpecSet CIDms2specs;
      SpecSet HCDms2specs;
      SpecSet ETDms2specs;
      m_outputSpectra->swapExtractSpectra(CIDms2specs, Spectrum::FragType_CID);
      m_outputSpectra->swapExtractSpectra(HCDms2specs, Spectrum::FragType_HCD);
      m_outputSpectra->swapExtractSpectra(ETDms2specs, Spectrum::FragType_ETD);

      // Don't need these in memory any more
      m_outputSpectra->clear();

      DEBUG_VAR(CIDms2specs.size());
      DEBUG_VAR(HCDms2specs.size());
      DEBUG_VAR(ETDms2specs.size());
      DEBUG_VAR(reversedPeaks);

      unsigned int idxUse = 0;
      bool execSuccess = true;

      if (CIDms2specs.size() > 0)
      {

        SpecSet CIDprms(CIDms2specs.size());

        // propagate reversed peaks
        vector<vector<bool> >* reversedPeaksCID = (vector<vector<bool> >*)0;
        if (reversedPeaks)
        {
          reversedPeaksCID = new vector<vector<bool> > (CIDms2specs.size());
        }

        // generate PRM spectra w/ CID model
        bool res = invokePepNovoCID(m_params,
                                    CIDms2specs,
                                    CIDprms,
                                    reversedPeaksCID);
        if (!res)
        {
          if (reversedPeaksCID)
          {
            delete reversedPeaksCID;
          }
          return false;
        }

        ms2Specs.swapAppendSpecSet(CIDms2specs, false);
        CIDms2specs.clear();

        // append CID PRMs
        prmSpecs.appendSpecSet(CIDprms, false);

        if (reversedPeaks)
        { // get reversed peaks
          for (unsigned int i = 0; i < reversedPeaksCID->size(); i++)
          {
            (*reversedPeaks)[idxUse] = (*reversedPeaksCID)[i];
            ++idxUse;
          }
          delete reversedPeaksCID;
        }
      }
      if (ETDms2specs.size() > 0)
      {
        SpecSet ETDprms(ETDms2specs.size());
        vector<vector<bool> >* reversedPeaksETD = (vector<vector<bool> >*)0;
        if (reversedPeaks)
        {
          reversedPeaksETD = new vector<vector<bool> > (ETDms2specs.size());
        }

        // generate PRM spectra with ETD model
        bool res = invokePepNovoETD(m_params,
                                    ETDms2specs,
                                    ETDprms,
                                    reversedPeaksETD);
        if (!res)
        {
          if (reversedPeaksETD)
          {
            delete reversedPeaksETD;
          }
          return false;
        }

        ms2Specs.swapAppendSpecSet(ETDms2specs, false);
        ETDms2specs.clear();

        prmSpecs.appendSpecSet(ETDprms, false);

        if (reversedPeaks)
        {
          for (unsigned int i = 0; i < reversedPeaksETD->size(); i++)
          {
            (*reversedPeaks)[idxUse] = (*reversedPeaksETD)[i];
            ++idxUse;
          }
          delete reversedPeaksETD;
        }

      }
      if (HCDms2specs.size() > 0)
      {
        SpecSet HCDprms(HCDms2specs.size());
        vector<vector<bool> >* reversedPeaksHCD = (vector<vector<bool> >*)0;
        if (reversedPeaks)
        {
          reversedPeaksHCD = new vector<vector<bool> > (HCDms2specs.size());
        }

        // generate PRM spectra with HCD model
        bool res = invokePepNovoHCD(m_params,
                                    HCDms2specs,
                                    HCDprms,
                                    reversedPeaksHCD);
        DEBUG_TRACE;
        if (!res)
        {
          if (reversedPeaksHCD)
          {
            delete reversedPeaksHCD;
          }
          return false;
        }

        ms2Specs.swapAppendSpecSet(HCDms2specs, false);
        HCDms2specs.clear();

        prmSpecs.appendSpecSet(HCDprms, false);

        if (reversedPeaks)
        {
          for (unsigned int i = 0; i < reversedPeaksHCD->size(); i++)
          {
            (*reversedPeaks)[idxUse] = (*reversedPeaksHCD)[i];
            ++idxUse;
          }
          delete reversedPeaksHCD;
        }
      }
    }

    m_outputSpectra->operator =(prmSpecs);

    DEBUG_TRACE;

    // Keep spectra in the same order as input MS/MS
    SpecSet orderedSpecs(m_outputSpectra->size());
    SpecSet orderedMs2Specs(ms2Specs.size());
    DEBUG_VAR(m_outputSpectra->size());
    DEBUG_VAR(ms2Specs.size())

    if (m_outputSpectra->size() != ms2Specs.size())
    {
      ERROR_MSG("Number of PRM spectra does not equal # of MS/MS spectra!!");
      abort();
    }

    for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
    {
      orderedSpecs[(*m_outputSpectra)[i].scan].swap((*m_outputSpectra)[i]);
      orderedMs2Specs[ms2Specs[i].scan].swap(ms2Specs[i]);
    }
    m_outputSpectra->swap(orderedSpecs);
    ms2Specs.swap(orderedMs2Specs);

    DEBUG_TRACE;

    if (m_params.getValueInt("BOOST_SILAC_PRMS", 0) > 0)
    {
      vector<vector<unsigned int> >* newClusters = new vector<vector<
          unsigned int> > (pairedScans.size());
      boostSilacPairs(m_outputSpectra, &ms2Specs, pairedScans, *newClusters);
      pairedScans = *newClusters;
      delete newClusters;
    }

    list<pair<unsigned, unsigned> > refList;
    pair<unsigned, unsigned> refIdx;
    if (mergeSamePrec)
    {
      // cluster spectra from the same precursor
      DEBUG_TRACE;
      vector<vector<unsigned int> >* newClusters = new vector<vector<
          unsigned int> > (pairedScans.size());

      DEBUG_VAR(m_outputSpectra->size());

      mergeSamePrecursorStaged(m_outputSpectra,
                               pairedScans,
                               *newClusters,
                               reversedPeaks);

      DEBUG_VAR(m_outputSpectra->size());
      m_outputClusters->data.clear();

      for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
      {
        m_outputClusters->data[i] = refList;
        for (unsigned int j = 0; j < (*newClusters)[i].size(); j++)
        {
          refIdx.first = 0;
          refIdx.second = (*newClusters)[i][j];
          m_outputClusters->data[i].push_back(refIdx);
        }
      }
      DEBUG_TRACE;
      delete newClusters;

    }
    else
    {
      DEBUG_TRACE;
      // otherwise, just generate single-ton clusters
      for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
      {
        m_outputClusters->data[i] = refList;
        refIdx.first = 0;
        refIdx.second = i;
        m_outputClusters->data[i].push_back(refIdx);
      }
    }

    if (reversedPeaks)
    {
      delete reversedPeaks;
    }

    DEBUG_TRACE;
    /*
     if (m_params.exists("PRM_SPEC_IDS"))
     {
     m_specIDs = new PeptideSpectrumMatchSet;
     string specIDFormat = m_params.getValue("PRM_SPEC_ID_FORMAT");
     bool loadSuccess = true;
     DEBUG_MSG("Loading " << specIDFormat << " results from <" << m_params.getValue("PRM_SPEC_IDS") << "> ...");
     if (specIDFormat == "InspecT")
     {
     if (!m_specIDs->loadInspectResultsFile(m_params.getValue("PRM_SPEC_IDS").c_str(),
     true))
     {
     loadSuccess = false;
     }
     }
     else if (specIDFormat == "MSGFDB")
     {
     if (!m_specIDs->loadMSGFDBResultsFile(m_params.getValue("PRM_SPEC_IDS").c_str(),
     false))
     {
     loadSuccess = false;
     }
     }
     else
     {
     ERROR_MSG("Found unsupported spectrum ID format \'" << specIDFormat << "\'");
     return false;
     }

     if (!loadSuccess)
     {
     ERROR_MSG("Failed to load " << m_params.getValue("PRM_SPEC_IDS") << " in " << specIDFormat << " format.");
     return false;
     }

     AAJumps jumps(1);
     DEBUG_TRACE;
     map<int, list<int> > clusterInfoMap;
     for (map<unsigned, list<pair<unsigned, unsigned> > >::iterator clustIt =
     m_outputClusters->data.begin(); clustIt
     != m_outputClusters->data.end(); clustIt++)
     {
     list<int> scans;

     for (list<pair<unsigned, unsigned> >::iterator scanIt =
     clustIt->second.begin(); scanIt != clustIt->second.end(); scanIt++)
     {
     scans.push_back(((int)scanIt->second) + 1);
     }
     clusterInfoMap[((int)clustIt->first) + 1] = scans;
     }
     DEBUG_TRACE;
     m_specIDs->cluster(clusterInfoMap);

     DEBUG_TRACE;

     int lastSpecIdx = 0;

     for (int i = 0; i < m_specIDs->size(); i++)
     {
     if ((*m_specIDs)[i]->m_annotation.length() == 0)
     {
     continue;
     }
     else
     {
     //(*m_specIDs)[i]->m_scanNum += 1;
     int specIdx = -1;//(*m_specIDs)[i]->m_scanNum;
     int j = lastSpecIdx;
     bool init = true;
     while (init || j != lastSpecIdx)
     {
     init = false;

     if ((j + 1) == (*m_specIDs)[i]->m_scanNum)
     {
     specIdx = j;
     lastSpecIdx = j;
     break;
     }
     j = (j == m_outputSpectra->size() - 1) ? 0 : j + 1;
     }
     if (specIdx < 0)
     {
     ERROR_MSG("Failed to locate spectrum " << (*m_specIDs)[i]->m_scanNum);
     abort();
     }
     string peptide = makeBracketsMods((*m_specIDs)[i]->m_annotation);
     const char* id = peptide.c_str();
     char just_letters[strlen(id) + 1];

     //get rid of special characters and modifications so we can match to protein
     getPepSeq(id, &just_letters[0]);

     string AAPeptideSeq = &just_letters[0];

     DEBUG_MSG("Assigning PRM spectrum for \'" << peptide << "\' to spectrum " << specIdx);

     vector<float> prmMasses;
     jumps.getPRMMasses(peptide, prmMasses, 0, true);
     vector<float> unModPrmMasses;
     jumps.getPRMMasses(AAPeptideSeq, unModPrmMasses, 0, true);

     (*m_outputSpectra)[specIdx].resize(prmMasses.size());
     for (unsigned int x = 0; x < prmMasses.size(); x++)
     {
     (*m_outputSpectra)[specIdx][x][0] = prmMasses[x];
     (*m_outputSpectra)[specIdx][x][1] = 10.0;
     }

     (*m_outputSpectra)[specIdx].setPeakTolerance(0.001);

     for (unsigned int x = 0; x < unModPrmMasses.size(); x++)
     {
     if ((*m_outputSpectra)[specIdx].findPeaks(unModPrmMasses[x]) < 0)
     {
     (*m_outputSpectra)[specIdx].insertPeak(unModPrmMasses[x],
     10.0,
     0.001);
     }
     }

     (*m_outputSpectra)[specIdx].addZPMpeaks(-1.0, 0, true);
     }
     }
     }
     */

    int rankFilt = m_params.getValueInt("PRM_RANK_FILTER", -1);

    if (rankFilt > 0)
    {
      DEBUG_MSG("Rank-filtering PRM spectra with K = " << rankFilt);
    }

    for (unsigned int i = 0; i < m_outputSpectra->size(); i++)
    {
      // assign spectra their native scan numbers
      (*m_outputSpectra)[i].scan
          = m_scanToInputScan[(*m_outputSpectra)[i].scan];

      (*m_outputSpectra)[i].msFragType = Spectrum::FragType_PRM;

      if (rankFilt > 0)
      {
        (*m_outputSpectra)[i].rankFilterPeaks(rankFilt);
      }

      if (m_params.getValueInt("CLUSTER_MIN_SIZE", 0) > 0)
      {
        (*m_outputSpectra)[i].scan = i + 1;
      }
    }

    DEBUG_TRACE;

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecPrmScoring::loadInputData(void)
  {

    if (ownOutput && m_outputSpectra == 0)
    {
      m_outputSpectra = new SpecSet;
    }

    // If we own the input, there is no point copying the input spectra during pre-processing
    if (ownInput && m_inputSpectra == 0)
    {
      m_inputSpectra = m_outputSpectra;
    }

    if (m_loader == 0)
    {
      m_loadedIndices = new vector<pair<int, int> > ;

      m_loader = new ExecMergeConvert(m_params,
                                      m_loadedIndices,
                                      m_inputSpectra,
                                      m_outputSpectra);
    }

    if (!m_loader->loadInputData())
    {
      return false;
    }

    if (m_inputSpectra->size() == 0)
    {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_SPECTRA?");
      return false;
    }

    /*
     m_scanToIdx.resize(m_inputSpectra->size());
     m_scanToInputScan.resize(m_inputSpectra->size());
     unsigned int total = 0;
     for (unsigned int i = 0; i < m_loader->m_recordedInput.size(); i++)
     {
     for (unsigned int j = 0; j < m_loader->m_recordedInput[i].second; j++)
     {
     m_scanToInputScan[total] = (*m_inputSpectra)[total].scan;
     (*m_inputSpectra)[total].scan = total;
     m_scanToIdx[total].first = i;
     m_scanToIdx[total].second = j;
     }
     }
     */

    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::saveInputData(std::vector<std::string> & filenames)
  {
    //SpecSet m_inputSpectra; // the input spectra
    /*
     if (m_params.exists("INPUT_SPECTRA")) {
     if (!ExecMergeConvert::saveSpecset(m_params.getValue("INPUT_SPECTRA"),
     m_inputSpectra)) {
     return false;
     }
     }
     */
    std::string paramFilename = getName() + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector

    return true;
  }

  //-----------------------------------------------------------------------------
  bool ExecPrmScoring::saveOutputData(void)
  {

    if (m_outputSpectra && m_params.exists("OUTPUT_SPECTRA"))
    {
      if (!ExecMergeConvert::saveSpecsetMultiple(m_params.getValue("OUTPUT_SPECTRA"),
                                                 m_outputSpectra))
      {
        return false;
      }
    }
    if (m_outputClusters && m_params.exists("OUTPUT_CLUSTERS_DIR"))
    {
      if (m_outputClusters->saveData(m_params.getValue("OUTPUT_CLUSTERS_DIR"))
          != 1)
      {
        ERROR_MSG("Failed to save spectra/clusterData.bin!");
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::loadOutputData(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecPrmScoring::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecPrmScoring::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("CORRECT_PM");
    VALIDATE_PARAM_EXIST("GUESS_CHARGE");
    if (m_params.exists("PEPNOVO_EXE_DIR"))
    {
      VALIDATE_PARAM_EXIST("PEPNOVO_EXE_DIR");
    }
    else
    {
      VALIDATE_PARAM_EXIST("EXE_DIR");
    }
    if (m_params.exists("PEPNOVO_MODEL"))
    {
      VALIDATE_PARAM_EXIST("PEPNOVO_INPUT_MGF");
      VALIDATE_PARAM_EXIST("PEPNOVO_OUTPUT_PRMS");
    }
    else
    {
      VALIDATE_PARAM_EXIST("PEPNOVO_OUTDIR");
    }

    m_isValid = true;
    return true;
  }

  bool ExecPrmScoring::invokePepNovoCID(ParameterList& params,
                                        SpecSet& inputSpectra,
                                        SpecSet& outputSpectra,
                                        vector<vector<bool> >* reversedPeaks)
  {

    DEBUG_MSG("Invoking PepNovo for input CID spectra...");

    string ftType = "FT";
    bool highAcc = ftType == params.getValue("INSTRUMENT_TYPE", "IT");

    ParameterList pepnovoParams(params);

    if (highAcc)
    {
      pepnovoParams.setValue("PEPNOVO_MODEL", "CID_FT_Tryp_z");
    }
    else
    {
      pepnovoParams.setValue("PEPNOVO_MODEL", "CID_IT_TRYP");
    }

    if ((!enforceDaTolerance) && highAcc)
    {
      pepnovoParams.setValue("TOLERANCE_PEAK", "0.04");
    }
    else if ((!enforceDaTolerance) && (!highAcc))
    {
      pepnovoParams.setValue("TOLERANCE_PEAK", "0.5");

    }

    string inFile(pepnovoParams.getValue("PEPNOVO_OUTDIR"));
    inFile += "/pepnovo_in_CID.mgf";
    pepnovoParams.setValue("PEPNOVO_INPUT_MGF", inFile);
    string outFile(pepnovoParams.getValue("PEPNOVO_OUTDIR"));
    outFile += "/pepnovo_out_CID.prms";
    pepnovoParams.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    bool res = invokePepNovo(pepnovoParams,
                             inputSpectra,
                             outputSpectra,
                             reversedPeaks);

    if (!res)
    {
      return false;
    }

    return true;
  }

  bool ExecPrmScoring::invokePepNovoHCD(ParameterList& params,
                                        SpecSet& inputSpectra,
                                        SpecSet& outputSpectra,
                                        vector<vector<bool> >* reversedPeaks)
  {

    DEBUG_MSG("Invoking PepNovo for input HCD spectra...");

    ParameterList pepnovoParams(params);

    string ftType = "FT";

    SpecSet HCDms2charge2(inputSpectra.size());
    unsigned int charge2Idx = 0;
    SpecSet HCDms2charge3(inputSpectra.size());
    unsigned int charge3Idx = 0;

    for (unsigned int i = 0; i < inputSpectra.size(); i++)
    {
      // extract charge 2 and charge >= 3 spectra so they can be scored with different models
      Spectrum* hcdSpec = &inputSpectra[i];
      if (hcdSpec->parentCharge <= 2)
      {
        HCDms2charge2[charge2Idx] = *hcdSpec;
        ++charge2Idx;
      }
      else
      {
        HCDms2charge3[charge3Idx] = *hcdSpec;
        ++charge3Idx;
      }
      // Only keep one copy of each MS/MS spectrum in memory
      hcdSpec->resize(0);
    }
    inputSpectra.resize(0);
    HCDms2charge2.resize(charge2Idx);
    HCDms2charge3.resize(charge3Idx);

    SpecSet HCDprms2(HCDms2charge2.size());
    SpecSet HCDprms3(HCDms2charge3.size());

    vector<vector<bool> >* reversedPeaksHCD2 = (vector<vector<bool> >*)0;
    vector<vector<bool> >* reversedPeaksHCD3 = (vector<vector<bool> >*)0;

    if (reversedPeaks)
    {
      reversedPeaksHCD2 = new vector<vector<bool> > (HCDms2charge2.size());
      reversedPeaksHCD3 = new vector<vector<bool> > (HCDms2charge3.size());
    }

    DEBUG_VAR(HCDms2charge2.size());
    DEBUG_VAR(HCDms2charge3.size());

    bool returnVal = true;

    ParameterList charge2Params(pepnovoParams);
    ParameterList charge3Params(pepnovoParams);

    charge2Params.setValue("PEPNOVO_MODEL", "HCD_Tryp_z");
    charge3Params.setValue("PEPNOVO_MODEL", "HCD_LysC_z");

    string inFile(charge2Params.getValue("PEPNOVO_OUTDIR"));
    inFile += "/pepnovo_in_HCD_2.mgf";
    charge2Params.setValue("PEPNOVO_INPUT_MGF", inFile);
    string outFile(charge2Params.getValue("PEPNOVO_OUTDIR"));
    outFile += "/pepnovo_out_HCD_2.prms";
    charge2Params.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    bool res = invokePepNovo(charge2Params,
                             HCDms2charge2,
                             HCDprms2,
                             reversedPeaksHCD2);

    if (!res)
    {
      return false;
    }

    inFile = charge3Params.getValue("PEPNOVO_OUTDIR");
    inFile += "/pepnovo_in_HCD_3.mgf";
    charge3Params.setValue("PEPNOVO_INPUT_MGF", inFile);
    outFile = charge3Params.getValue("PEPNOVO_OUTDIR");
    outFile += "/pepnovo_out_HCD_3.prms";
    charge3Params.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    res = invokePepNovo(charge3Params,
                        HCDms2charge3,
                        HCDprms3,
                        reversedPeaksHCD3);

    if (!res)
    {
      return false;
    }

    outputSpectra = HCDprms2;
    outputSpectra.appendSpecSet(HCDprms3, false);

    if (reversedPeaks)
    {
      reversedPeaks->resize(HCDprms2.size() + HCDprms3.size());
      unsigned int idxUse = 0;
      for (unsigned int i = 0; i < HCDprms2.size(); i++)
      {
        (*reversedPeaks)[idxUse] = (*reversedPeaksHCD2)[i];
        ++idxUse;
      }
      for (unsigned int i = 0; i < HCDprms3.size(); i++)
      {
        (*reversedPeaks)[idxUse] = (*reversedPeaksHCD3)[i];
        ++idxUse;
      }
    }

    if (reversedPeaksHCD2)
    {
      delete reversedPeaksHCD2;
    }
    if (reversedPeaksHCD3)
    {
      delete reversedPeaksHCD3;
    }

    // Restore input MS/MS spectra
    inputSpectra.resize(HCDms2charge2.size() + HCDms2charge3.size());
    for (unsigned int i = 0; i < HCDms2charge2.size(); i++)
    {
      inputSpectra[i] = HCDms2charge2[i];
      HCDms2charge2[i].resize(0);
    }
    HCDms2charge2.resize(0);

    for (unsigned int i = 0; i < HCDms2charge3.size(); i++)
    {
      inputSpectra[HCDms2charge2.size() + i] = HCDms2charge3[i];
      HCDms2charge3[i].resize(0);
    }
    return true;
  }

  bool ExecPrmScoring::invokePepNovoETD(ParameterList& params,
                                        SpecSet& inputSpectra,
                                        SpecSet& outputSpectra,
                                        vector<vector<bool> >* reversedPeaks)
  {

    DEBUG_MSG("Invoking PepNovo for input ETD spectra...");

    ParameterList pepnovoParams(params);

    string ftType = "FT";
    bool highAcc = ftType == pepnovoParams.getValue("INSTRUMENT_TYPE", "IT");

    SpecSet ETDms2charge2(inputSpectra.size());
    unsigned int charge2Idx = 0;
    SpecSet ETDms2charge3(inputSpectra.size());
    unsigned int charge3Idx = 0;

    for (unsigned int i = 0; i < inputSpectra.size(); i++)
    {
      // extract charge 2 and charge >= 3 spectra so they can be scored with different models
      Spectrum* etdSpec = &inputSpectra[i];
      if (etdSpec->parentCharge <= 2)
      {
        ETDms2charge2[charge2Idx] = *etdSpec;
        ++charge2Idx;
      }
      else
      {
        ETDms2charge3[charge3Idx] = *etdSpec;
        ++charge3Idx;
      }
      // Only keep one copy of each MS/MS spectrum in memory
      etdSpec->resize(0);
    }
    inputSpectra.resize(0);

    ETDms2charge2.resize(charge2Idx);
    ETDms2charge3.resize(charge3Idx);

    SpecSet ETDprms2(ETDms2charge2.size());
    SpecSet ETDprms3(ETDms2charge3.size());

    vector<vector<bool> >* reversedPeaksETD2 = (vector<vector<bool> >*)0;
    vector<vector<bool> >* reversedPeaksETD3 = (vector<vector<bool> >*)0;

    if (reversedPeaks)
    {
      reversedPeaksETD2 = new vector<vector<bool> > (ETDms2charge2.size());
      reversedPeaksETD3 = new vector<vector<bool> > (ETDms2charge3.size());
    }

    DEBUG_VAR(ETDms2charge2.size());
    DEBUG_VAR(ETDms2charge3.size());

    bool returnVal = true;

    ParameterList charge2Params(pepnovoParams);
    ParameterList charge3Params(pepnovoParams);

    charge2Params.setValue("PEPNOVO_MODEL", "ETD_FT_Tryp_z");
    charge3Params.setValue("PEPNOVO_MODEL", "ETD_FT_LysC_z");

    string inFile(charge2Params.getValue("PEPNOVO_OUTDIR"));
    inFile += "/pepnovo_in_ETD_2.mgf";
    charge2Params.setValue("PEPNOVO_INPUT_MGF", inFile);
    string outFile(charge2Params.getValue("PEPNOVO_OUTDIR"));
    outFile += "/pepnovo_out_ETD_2.prms";
    charge2Params.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    bool res = invokePepNovo(charge2Params,
                             ETDms2charge2,
                             ETDprms2,
                             reversedPeaksETD2);

    if (!res)
    {
      return false;
    }

    inFile = charge3Params.getValue("PEPNOVO_OUTDIR");
    inFile += "/pepnovo_in_ETD_3.mgf";
    charge3Params.setValue("PEPNOVO_INPUT_MGF", inFile);
    outFile = charge3Params.getValue("PEPNOVO_OUTDIR");
    outFile += "/pepnovo_out_ETD_3.prms";
    charge3Params.setValue("PEPNOVO_OUTPUT_PRMS", outFile);

    res = invokePepNovo(charge3Params,
                        ETDms2charge3,
                        ETDprms3,
                        reversedPeaksETD3);

    if (!res)
    {
      return false;
    }

    outputSpectra = ETDprms2;
    outputSpectra.appendSpecSet(ETDprms3, false);

    if (reversedPeaks)
    {
      reversedPeaks->resize(ETDprms2.size() + ETDprms3.size());
      unsigned int idxUse = 0;
      for (unsigned int i = 0; i < ETDprms2.size(); i++)
      {
        (*reversedPeaks)[idxUse] = (*reversedPeaksETD2)[i];
        ++idxUse;
      }
      for (unsigned int i = 0; i < ETDprms3.size(); i++)
      {
        (*reversedPeaks)[idxUse] = (*reversedPeaksETD3)[i];
        ++idxUse;
      }
    }

    if (reversedPeaksETD2)
    {
      delete reversedPeaksETD2;
    }
    if (reversedPeaksETD3)
    {
      delete reversedPeaksETD3;
    }

    // Restore input MS/MS spectra
    inputSpectra.resize(ETDms2charge2.size() + ETDms2charge3.size());
    for (unsigned int i = 0; i < ETDms2charge2.size(); i++)
    {
      inputSpectra[i] = ETDms2charge2[i];
      ETDms2charge2[i].resize(0);
    }
    ETDms2charge2.resize(0);

    for (unsigned int i = 0; i < ETDms2charge3.size(); i++)
    {
      inputSpectra[ETDms2charge2.size() + i] = ETDms2charge3[i];
      ETDms2charge3[i].resize(0);
    }
    return true;
  }

  bool ExecPrmScoring::invokePepNovo(ParameterList& params,
                                     SpecSet& inputSpectra,
                                     SpecSet& outputSpectra,
                                     vector<vector<bool> >* reversedPeaks)
  {

    DEBUG_MSG("Invoking PepNovo ...");
    //DEBUG_VAR(inputSpectra.size());

    //DEBUG_VAR(params.getValue("TOLERANCE_PEAK"));

    string pepnovoMode;
    if (params.exists("MIN_SPECTRUM_QUALITY"))
    {
      pepnovoMode = "-min_filter_prob ";
      pepnovoMode += params.getValue("MIN_SPECTRUM_QUALITY");
    }
    else
    {
      pepnovoMode = "-no_quality_filter";
    }

    // TODO: Lower case the value?
    bool guessedPM = false;
    if (params.getValue("CORRECT_PM") == "yes")
    {
      pepnovoMode += " -correct_pm";
      guessedPM = true;
    }
    else
    {
      pepnovoMode += " -use_spectrum_mz";
    }

    // TODO: Lower case the value?
    bool guessedCharge = true;
    if (params.getValue("GUESS_CHARGE") == "no")
    {
      pepnovoMode += " -use_spectrum_charge";
      guessedCharge = false;
    }

    if (params.exists("PEPNOVO_PTMS"))
    {
      pepnovoMode += " -PTMs M+16:C+57:";
      pepnovoMode += params.getValue("PEPNOVO_PTMS");
    }
    else
    {
      pepnovoMode += " -PTMs M+16:C+57";
    }

    //DEBUG_VAR(pepnovoMode);
    string modelName = params.getValue("PEPNOVO_MODEL");

    string exeDir = params.getValue("EXE_DIR");
    string pepNovoModelDir(exeDir);

    // pepnovo models directory
    pepNovoModelDir = exeDir;
    if (params.exists("PEPNOVO_MODEL_DIR"))
    {
      pepNovoModelDir = params.getValue("PEPNOVO_MODEL_DIR");
    }
    else
    {
      pepNovoModelDir += "/resources/Models_pepnovo";
    }

    // pepnovo executable directory
    if (params.exists("PEPNOVO_EXE_DIR"))
    {
      exeDir = params.getValue("PEPNOVO_EXE_DIR");
    }

    string pepnovoCmd(exeDir);
    pepnovoCmd += "/PepNovo_bin -prm_norm -model_dir ";
    pepnovoCmd += pepNovoModelDir;
    pepnovoCmd += " -model ";
    pepnovoCmd += modelName;
    pepnovoCmd += " -file ";
    pepnovoCmd += params.getValue("PEPNOVO_INPUT_MGF");
    pepnovoCmd += " -fragment_tolerance ";
    pepnovoCmd += params.getValue("TOLERANCE_PEAK");
    if (params.exists("TOLERANCE_PM"))
    {
      pepnovoCmd += " -pm_tolerance ";
      pepnovoCmd += params.getValue("TOLERANCE_PM");
    }
    pepnovoCmd += " -digest NON_SPECIFIC ";
    pepnovoCmd += pepnovoMode;
    pepnovoCmd += " > ";
    pepnovoCmd += params.getValue("PEPNOVO_OUTPUT_PRMS");

    // Remove spectra with 5 or less peaks.
    checkSpecset(inputSpectra);

    map<unsigned int, unsigned int> scanRef;
    for (unsigned int i = 0; i < inputSpectra.size(); i++)
    {
      scanRef[inputSpectra[i].scan] = i;
    }

    // Prepare a copy of MS/MS spectra for PepNovo input

    vector<unsigned short> oldCharges(inputSpectra.size());

    for (int i = 0; i < inputSpectra.size(); i++)
    {
      Spectrum* spec = &inputSpectra[i];
      oldCharges[i] = spec->parentCharge;
      // ETD spectra have a charge 4 model while HCD spectra only have a charge 3 model
      if (modelName.find("ETD") != string::npos && spec->parentCharge > 4)
      {
        spec->setCharge(4);
      }
      else if (modelName.find("ETD") == string::npos && spec->parentCharge > 3)
      {
        spec->setCharge(3);
      }
    }

    DEBUG_VAR(enforceDaTolerance);

    ifstream infile;
    infile.open(params.getValue("PEPNOVO_OUTPUT_PRMS").c_str(), ios::binary);

    if (params.getValueInt("SKIP_PEPNOVO_INVOKE", 0) == 0 || !infile.is_open())
    {
      if (infile.is_open())
      {
        infile.close();
      }

      string tempPklbin = params.getValue("PEPNOVO_INPUT_MGF");
      tempPklbin += ".pklbin";

      // Save a copy of the input spectra so we can clear memory before calling PepNovo
      DEBUG_MSG("Saving temporary file \'" << tempPklbin << "\'");
      if (!inputSpectra.SaveSpecSet_pklbin(tempPklbin.c_str()))
      {
        ERROR_MSG("Failed to save \'" << tempPklbin << "\'");
        return false;
      }

      // invoke PepNovo via command line
      DEBUG_MSG("Saving temporary file \'" << params.getValue("PEPNOVO_INPUT_MGF") << "\'");
      if (!inputSpectra.SaveSpecSet_mgf(params.getValue("PEPNOVO_INPUT_MGF").c_str()))
      {
        return false;
      }

      inputSpectra.clear();

      DEBUG_VAR(pepnovoCmd);
      int status = spsSystem(pepnovoCmd.c_str());
      if (status != 0)
      {
        string errorString = "Executing ";
        errorString += exeDir;
        errorString += "/Pepnovo_bin!";
        ERROR_MSG(errorString);
        return false;
      }

      // Reload MS/MS spectra
      DEBUG_MSG("Loading temporary file \'" << tempPklbin << "\'");
      if (!inputSpectra.loadPklBin(tempPklbin.c_str()))
      {
        ERROR_MSG("Failed to load \'" << tempPklbin << "\'");
        return false;
      }

      // Remove temporary file
      if (remove(tempPklbin.c_str()) != 0)
      {
        int aa = errno;
        string aux = strerror(aa);
        ERROR_MSG("Failed to remove \'" << tempPklbin << "\'");
        ERROR_MSG("error = " << aux);
        //        return false;
      }
      else
      {
        DEBUG_MSG("Successfully removed \'" << tempPklbin << "\'");
      }
    }
    if (infile.is_open())
    {
      infile.close();
    }

    SpecSet * pepnovoSpectra = new SpecSet(inputSpectra.size());

    vector<vector<string> >* prmOrigins = 0;
    if (reversedPeaks)
    {
      prmOrigins = new vector<vector<string> > (pepnovoSpectra->size());
      reversedPeaks->resize(inputSpectra.size());
    }

    DEBUG_MSG("Loading temporary file \'" << params.getValue("PEPNOVO_OUTPUT_PRMS") << "\'");
    if (!pepnovoSpectra->LoadSpecSet_prmsv3(params.getValue("PEPNOVO_OUTPUT_PRMS").c_str(),
                                            prmOrigins))
    {
      return false;
    }

    if (params.getValueInt("SKIP_PEPNOVO_INVOKE", 0) == 0)
    {
      // remove temporary input/output files if they were generated
      if (remove(params.getValue("PEPNOVO_INPUT_MGF").c_str()) != 0)
      {
        int aa = errno;
        string aux = strerror(aa);
        ERROR_MSG("Failed to remove \'" << params.getValue("PEPNOVO_INPUT_MGF") << "\'");
        ERROR_MSG("error = " << aux);
        //        return false;
      }
      else
      {
        DEBUG_MSG("Successfully removed \'" << params.getValue("PEPNOVO_INPUT_MGF") << "\'");
      }

      if (remove(params.getValue("PEPNOVO_OUTPUT_PRMS").c_str()) != 0)
      {
        int aa = errno;
        string aux = strerror(aa);
        ERROR_MSG("Failed to remove \'" << params.getValue("PEPNOVO_OUTPUT_PRMS") << "\'");
        ERROR_MSG("error = " << aux);
        //        return false;
      }
      else
      {
        DEBUG_MSG("Successfully removed \'" << params.getValue("PEPNOVO_OUTPUT_PRMS") << "\'");
      }
    }

    //DEBUG_VAR(pepnovoSpectra->size());

    outputSpectra.resize(inputSpectra.size());

    for (unsigned int i = 0; i < inputSpectra.size(); i++)
    {
      inputSpectra[i].setCharge(oldCharges[i]);
      outputSpectra[i].resize(0);
      outputSpectra[i].copyNP(inputSpectra[i]);
    }
    for (unsigned int i = 0; i < pepnovoSpectra->size(); i++)
    {
      unsigned int specIdx = scanRef[(*pepnovoSpectra)[i].scan];

      // Copy PepNovo output spectra

      outputSpectra[specIdx].mergeClosestPeaks((*pepnovoSpectra)[i], 0);

      if (guessedPM)
      {
        // Let PepNovo set parent mass if we asked it to guess it
        outputSpectra[specIdx].setParentMass((*pepnovoSpectra)[i].parentMass);
      }

      if (guessedCharge || outputSpectra[specIdx].parentCharge == 0)
      {
        // Let PepNovo set the parent charge if we asked it to guess it
        outputSpectra[specIdx].setCharge((*pepnovoSpectra)[i].parentCharge);
      }

      Spectrum* ms2Spec = &inputSpectra[specIdx];
      Spectrum* prmSpec = &outputSpectra[specIdx];

      if (reversedPeaks && !enforceDaTolerance)
      {
        // Correct peak tolerances by reading last column of PepNovo output
        vector<bool>* reversedPeaksLoc = &(*reversedPeaks)[specIdx];
        vector<string>* prmOriginsLoc = &(*prmOrigins)[i];
        setOriginMasses(prmSpec, reversedPeaksLoc, ms2Spec, prmOriginsLoc);
        unsigned int revIdxUse = 0;
        list<int> removeIdxs;
        for (unsigned int j = 0; j < prmSpec->size(); j++)
        {
          if ((*prmSpec)[j][1] < -1.0)
          {
            removeIdxs.push_back(j);
          }
          else
          {
            (*reversedPeaksLoc)[revIdxUse] = (*reversedPeaksLoc)[j];
            ++revIdxUse;
          }
        }
        reversedPeaksLoc->resize(revIdxUse);
        prmSpec->removePeaks(removeIdxs);
      }
      else if ((!reversedPeaks) && !enforceDaTolerance)
      {
        vector<string>* prmOriginsLoc = &(*prmOrigins)[i];
        setOriginTolerances(prmSpec, ms2Spec, prmOriginsLoc);
        prmSpec->filterLowIntensity(-1.0);
      }
      else
      {
        prmSpec->setPeakTolerance(params.getValueFloat("TOLERANCE_PEAK"));
        prmSpec->filterLowIntensity(-1.0);
      }
    }

    if (prmOrigins)
    {
      delete prmOrigins;
    }

    delete pepnovoSpectra;

    //DEBUG_VAR(outputSpectra.size());

    return true;
  }

  void ExecPrmScoring::setOriginMasses(Spectrum* prmSpec,
                                       vector<bool>* reversedPeaks,
                                       Spectrum* ms2Spec,
                                       vector<string>* prmOriginPeaks)
  {

    reversedPeaks->resize(prmSpec->size());

    if (prmSpec->size() == 0)
    {
      return;
    }
    /*
     DEBUG_VAR(prmSpec->size());
     DEBUG_VAR(reversedPeaks->size());
     DEBUG_VAR(prmOriginPeaks->size());
     DEBUG_VAR(ms2Spec->size());

     prmSpec->output(cout);
     */

    Spectrum expectedPrmOffsets;
    expectedPrmOffsets.resize(8);
    expectedPrmOffsets[0][0] = 1.00728;
    expectedPrmOffsets[1][0] = -17.00328;
    expectedPrmOffsets[2][0] = -26.98710;
    expectedPrmOffsets[3][0] = -16.01927;
    expectedPrmOffsets[4][0] = 18.03430;
    expectedPrmOffsets[5][0] = 17.02232;
    expectedPrmOffsets[6][0] = -25.98199;
    expectedPrmOffsets[7][0] = 9.52079;
    expectedPrmOffsets.sortPeaks();

    Spectrum expectedSrmOffsets;
    expectedSrmOffsets.resize(5);
    expectedSrmOffsets[0][0] = 19.01840;
    expectedSrmOffsets[1][0] = 1.00784;
    expectedSrmOffsets[2][0] = 10.01284;
    expectedSrmOffsets[3][0] = 2.99000;
    expectedSrmOffsets[4][0] = 1.99864;
    expectedSrmOffsets.sortPeaks();

    prmSpec->setTolerance(0, 0.001);
    prmSpec->setTolerance(prmSpec->size() - 1, ms2Spec->parentMassTol);

    map<float, bool> reversedPeakRef;
    reversedPeakRef[prmSpec->front()->operator [](0)] = false;
    reversedPeakRef[prmSpec->back()->operator [](0)] = false;

    for (unsigned int i = 1; i < prmSpec->size() - 1; i++)
    {

      vector<string> originInfo(4);
      splitText((*prmOriginPeaks)[i].c_str(), originInfo, ",");
      float ms2mass = getFloat(originInfo[0].c_str());
      int orient = getInt(originInfo[1].c_str());
      int charge = getInt(originInfo[2].c_str());
      float massOffset = getFloat(originInfo[3].c_str());

      float srmOffset = 0;

      bool isSuf = false;

      if (orient == 0)
      { // from prefix ion
        int closestIdx = expectedPrmOffsets.findClosest(massOffset);
        if (MZRange::EqualWithinRange(massOffset,
                                      expectedPrmOffsets[closestIdx][0],
                                      0.001))
        {
          massOffset = expectedPrmOffsets[closestIdx][0];
        }
      }
      else
      { // suffix ion
        isSuf = true;
        int closestIdx = expectedSrmOffsets.findClosest(massOffset);
        if (MZRange::EqualWithinRange(massOffset,
                                      expectedSrmOffsets[closestIdx][0],
                                      0.001))
        {
          massOffset = expectedSrmOffsets[closestIdx][0];
        }
        srmOffset = (ms2Spec->msFragType == Spectrum::FragType_ETD) ? 0.0
            - AAJumps::massNH : AAJumps::massH2O;

      }

      int m2Idx = ms2Spec->findClosest(ms2mass);
      float massUse = (*ms2Spec)[m2Idx][0];

      if (charge > 1)
      {
        massUse = DeconvSpectrum::GetMonoisotopicMass(massUse, charge);
      }
      (*prmSpec)[i][0] = massUse - massOffset + srmOffset;
      prmSpec->setTolerance(i, ms2Spec->getTolerance(m2Idx) * ((float)charge));

      reversedPeakRef[(*prmSpec)[i][0]] = isSuf;

    }
    prmSpec->sortPeaks();

    for (unsigned int i = 0; i < prmSpec->size(); i++)
    {
      (*reversedPeaks)[i] = reversedPeakRef[(*prmSpec)[i][0]];
    }
  }

  void ExecPrmScoring::setOriginTolerances(Spectrum* prmSpec,
                                           Spectrum* ms2Spec,
                                           vector<string>* prmOriginPeaks)
  {
    if (prmSpec->size() == 0)
    {
      return;
    }

    prmSpec->setTolerance(0, 0.001);
    prmSpec->setTolerance(prmSpec->size() - 1, ms2Spec->parentMassTol);

    for (unsigned int i = 1; i < prmSpec->size() - 1; i++)
    {

      vector<string> originInfo(4);
      splitText((*prmOriginPeaks)[i].c_str(), originInfo, ",");
      float ms2mass = getFloat(originInfo[0].c_str());
      int orient = getInt(originInfo[1].c_str());
      int charge = getInt(originInfo[2].c_str());

      int m2Idx = ms2Spec->findClosest(ms2mass);

      float locms2Tol = (ms2Spec->getTolerance(m2Idx) + ((float)charge))
          + 0.001;

      if (orient == 0)
      { // from prefix ion
        prmSpec->setTolerance(i, locms2Tol);
      }
      else
      { // suffix ion
        prmSpec->setTolerance(i, locms2Tol + ms2Spec->parentMassTol);
      }
    }
  }

  void ExecPrmScoring::boostSilacPairs(SpecSet* inOutSpecs,
                                       SpecSet* ms2Spectra,
                                       vector<vector<unsigned int> >& inClusteredScans,
                                       vector<vector<unsigned int> >& outClusteredScans)
  {
    float modMassR = 10.008269;
    float modMassK = 8.014199;
    const unsigned int scanRange = 100;
    unsigned int minNumMP = 4 + 2;

    unsigned int numSilacSpecs = 0;

    float minRatioClust = m_params.getValueFloat("PRM_CLUSTER_RATIO", 0.72);

    bool filterNonBoost = m_params.getValueInt("FILTER_NONSILAC_PRMS", 1) > 0;

    DEBUG_VAR(filterNonBoost);

    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");

    DEBUG_VAR(inOutSpecs->size());

    vector<bool> usedSpecs(inOutSpecs->size(), false);

    PairedSpecSet cluster;

    SpecSet outSpecs(inOutSpecs->size());
    unsigned int idxUse = 0;
    unsigned int numPairs = 0;

    outClusteredScans = inClusteredScans;

    unsigned int specIdx = 0;
    unsigned int pairedIdx;

    list<unsigned int> lightIdxs;
    list<unsigned int> heavyIdxs;
    list<unsigned int> startIdxs;

    SpecSet locMS2Specs(scanRange);
    vector<unsigned int> locIdxMap(scanRange);

    SpecSet emptySpecs(scanRange);
    vector<TwoValues<float> > ratios;
    vector<TwoValues<float> > means;
    vector<float> varTerms;
    list<vector<float> > alignStats;
    vector<vector<float> > specStats;
    std::vector<unsigned int> idxKept;
    std::vector<TwoValues<float> > pvalues;

    ParameterList filterParams;
    filterParams.setValue("MAX_SHIFT", "12.0");
    filterParams.setValue("TOLERANCE_PEAK", m_params.getValue("TOLERANCE_PEAK"));
    filterParams.setValue("TOLERANCE_PM", m_params.getValue("TOLERANCE_PM"));
    filterParams.setValue("PAIRS_MIN_COSINE",
                          m_params.getValue("MIN_SILAC_COSINE", "0.4"));
    filterParams.setValue("MIN_NUM_MATCHED_PEAKS", "4");
    filterParams.setValue("SPEC_TYPE_MSMS", "1");
    filterParams.setValue("PAIRS_MATCH_MODE", "cosine");

    float srmOffset, curPMass, nextPMass, curPMtol, nextPMtol, lightPMass,
        heavyPMass;

    while (specIdx < inOutSpecs->size())
    {
      if (usedSpecs[specIdx] || (*inOutSpecs)[specIdx].size() == 0)
      {
        specIdx++;
        continue;
      }

      unsigned int specIdxLoc = specIdx;
      for (unsigned int i = 0; i < scanRange && specIdxLoc < inOutSpecs->size(); i++)
      {
        locIdxMap[i] = specIdxLoc;
        if (usedSpecs[specIdxLoc] || (*inOutSpecs)[specIdxLoc].size() == 0)
        {
          locMS2Specs[i].resize(0);
          specIdxLoc++;
          continue;
        }
        //DEBUG_VAR((*ms2Spectra)[specIdxLoc].parentMass);
        //DEBUG_VAR((*inOutSpecs)[specIdxLoc].parentMass);

        locMS2Specs[i] = (*ms2Spectra)[specIdxLoc];
        for (unsigned int j = 0; j < locMS2Specs[i].size(); j++)
        {
          locMS2Specs[i][j][1] = sqrt(locMS2Specs[i][j][1]);
        }
        locMS2Specs[i].normalize2();
        specIdxLoc++;
      }
      if (specIdxLoc == inOutSpecs->size())
      {
        locMS2Specs.resize(specIdxLoc - specIdx + 1);
        locIdxMap.resize(specIdxLoc - specIdx + 1);
      }

      // set of matched heavy/light pairs
      SpectrumPairSet filteredPairs;

      //DEBUG_VAR(specIdx);

      // Let ExecFilterPairs do cosine matching to find heavy/light pairs
      ExecFilterPairs moduleAlign(filterParams,
                                  &emptySpecs,
                                  &locMS2Specs,
                                  &filteredPairs,
                                  &ratios,
                                  &means,
                                  &varTerms,
                                  &alignStats,
                                  &specStats,
                                  &idxKept,
                                  &pvalues);

      if (!moduleAlign.invoke())
      {
        ERROR_MSG("Failed to invoke ExecFilterPairs!!!");
        abort();
      }

      //DEBUG_VAR(filteredPairs.size());

      curPMass = (*inOutSpecs)[specIdx].parentMass;
      curPMtol = (*inOutSpecs)[specIdx].parentMassTol;

      startIdxs.clear();
      lightIdxs.clear();
      heavyIdxs.clear();

      startIdxs.push_back(specIdx);

      srmOffset = -1.0;

      //DEBUG_VAR(specIdx);
      //DEBUG_VAR(curPMass);
      //DEBUG_VAR(curPMtol);
      filteredPairs.sort_pairs();

      for (unsigned int i = 0; i < filteredPairs.size(); i++)
      {
        unsigned int specIdx1 = locIdxMap[filteredPairs[i].spec1];
        unsigned int specIdx2 = locIdxMap[filteredPairs[i].spec2];

        if (specIdx1 == specIdx)
        {
          nextPMass = (*inOutSpecs)[specIdx2].parentMass;
          nextPMtol = (*inOutSpecs)[specIdx2].parentMassTol;
          pairedIdx = specIdx2;
        }
        else if (specIdx2 == specIdx)
        {
          nextPMass = (*inOutSpecs)[specIdx1].parentMass;
          nextPMtol = (*inOutSpecs)[specIdx1].parentMassTol;
          pairedIdx = specIdx1;
        }
        else
        {
          continue;
        }

        bool samePM = MZRange::EqualWithinRange(nextPMass, curPMass, curPMtol
            + nextPMtol);
        bool nextHR = MZRange::EqualWithinRange(nextPMass,
                                                curPMass + modMassR,
                                                curPMtol + nextPMtol);
        bool nextHK = MZRange::EqualWithinRange(nextPMass,
                                                curPMass + modMassK,
                                                curPMtol + nextPMtol);
        bool nextLR = MZRange::EqualWithinRange(nextPMass + modMassR,
                                                curPMass,
                                                curPMtol + nextPMtol);
        bool nextLK = MZRange::EqualWithinRange(nextPMass + modMassK,
                                                curPMass,
                                                curPMtol + nextPMtol);

        bool foundSim = samePM || nextHR || nextHK || nextLR || nextLK;

        if (!foundSim)
        {
          continue;
        }

        /*
         if (pairedIdx == 12395)
         {
         DEBUG_VAR(specIdx);
         DEBUG_VAR(pairedIdx);
         DEBUG_VAR(filteredPairs[i].score1);
         DEBUG_VAR((*inOutSpecs)[specIdx].parentMass);
         DEBUG_VAR((*inOutSpecs)[pairedIdx].parentMass);
         }
         */

        //DEBUG_VAR(filteredPairs[i].score1);

        if (srmOffset < 0)
        {
          if (samePM)
          {
            startIdxs.push_back(i);
          }
          else if (nextHR)
          {
            srmOffset = modMassR;
            lightIdxs = startIdxs;
            heavyIdxs.push_back(pairedIdx);
            lightPMass = curPMass;
            heavyPMass = curPMass + modMassR;
          }
          else if (nextHK)
          {
            srmOffset = modMassK;
            lightIdxs = startIdxs;
            heavyIdxs.push_back(pairedIdx);
            lightPMass = curPMass;
            heavyPMass = curPMass + modMassK;
          }
          else if (nextLR)
          {
            srmOffset = modMassR;
            heavyIdxs = startIdxs;
            lightIdxs.push_back(pairedIdx);
            heavyPMass = curPMass;
            lightPMass = curPMass - modMassR;
          }
          else
          {
            srmOffset = modMassK;
            heavyIdxs = startIdxs;
            lightIdxs.push_back(pairedIdx);
            heavyPMass = curPMass;
            lightPMass = curPMass - modMassK;
          }
        }
        else if (MZRange::EqualWithinRange(nextPMass, heavyPMass, curPMtol
            + nextPMtol))
        {
          heavyIdxs.push_back(pairedIdx);
        }
        else if (MZRange::EqualWithinRange(nextPMass, lightPMass, curPMtol
            + nextPMtol))
        {
          lightIdxs.push_back(pairedIdx);
        }
        else
        {
          continue;
        }

      }

      if (lightIdxs.size() == 0 || heavyIdxs.size() == 0)
      {
        if (!filterNonBoost)
        {
          outSpecs[specIdx] = (*inOutSpecs)[specIdx];
          outSpecs[specIdx].msFragType = Spectrum::FragType_CID;
        }

        specIdx++;
        continue;
      }

      SpecSet pairedSpecs(lightIdxs.size() + heavyIdxs.size());
      unsigned int locIdx = 0;
      vector<unsigned int> boostedIdxs(lightIdxs.size() + heavyIdxs.size());
      for (list<unsigned int>::const_iterator idxIt = lightIdxs.begin(); idxIt
          != lightIdxs.end(); idxIt++)
      {
        pairedSpecs[locIdx] = (*inOutSpecs)[*idxIt];
        pairedSpecs[locIdx].msFragType = Spectrum::FragType_CID;
        boostedIdxs[locIdx] = *idxIt;
        usedSpecs[*idxIt] = true;
        locIdx++;
      }
      for (list<unsigned int>::const_iterator idxIt = heavyIdxs.begin(); idxIt
          != heavyIdxs.end(); idxIt++)
      {
        pairedSpecs[locIdx] = (*inOutSpecs)[*idxIt];
        pairedSpecs[locIdx].msFragType = Spectrum::FragType_ETD;
        boostedIdxs[locIdx] = *idxIt;
        usedSpecs[*idxIt] = true;
        locIdx++;
      }

      numSilacSpecs += lightIdxs.size() + heavyIdxs.size();

      //DEBUG_VAR(specIdx);
      //DEBUG_VAR(lightIdxs.size());
      //DEBUG_VAR(heavyIdxs.size());

      cluster.initialize(&pairedSpecs, true, 0, peakTol, srmOffset);

      SpecSet boostedSpecs;

      cluster.boostPRMs(boostedSpecs);

      //DEBUG_VAR(boostedSpecs.size());

      for (locIdx = 0; locIdx < boostedSpecs.size(); locIdx++)
      {
        if (filterNonBoost)
        {
          outSpecs[idxUse] = boostedSpecs[locIdx];
          outSpecs[idxUse].msFragType = Spectrum::FragType_CID;
          outClusteredScans[idxUse] = inClusteredScans[boostedIdxs[locIdx]];
          idxUse++;
        }
        else
        {
          outSpecs[boostedIdxs[locIdx]] = boostedSpecs[locIdx];
          outSpecs[boostedIdxs[locIdx]].msFragType = Spectrum::FragType_CID;
        }
      }

      specIdx++;
      numPairs++;
    }

    if (filterNonBoost)
    {
      DEBUG_VAR(idxUse);
      outSpecs.resize(idxUse);
      outClusteredScans.resize(idxUse);
    }
    inOutSpecs->operator =(outSpecs);

    DEBUG_VAR(inOutSpecs->size());

    DEBUG_MSG("Found " << numSilacSpecs << " silac spectra and " << numPairs << " components");
  }

  struct SortPairs : public std::binary_function<pair<float, TwoValues<
      unsigned int> > , pair<float, TwoValues<unsigned int> > , bool>
  {
    bool operator()(pair<float, TwoValues<unsigned int> > left, pair<float,
        TwoValues<unsigned int> > right) const
    {
      return left.first > right.first;
    }
    ;
  };

  void ExecPrmScoring::mergeSamePrecursorStaged(SpecSet* inOutSpecs,
                                                vector<vector<unsigned int> >& inPairedScans,
                                                vector<vector<unsigned int> >& outClusteredScans,
                                                vector<vector<bool> >* reversedPeaks)
  {

    DEBUG_MSG("Merging spectra from the same precursor with hierarchical method");

    // Minimum allowable overlapping score between clustered PRM spectra
    float minRatioClust = m_params.getValueFloat("PRM_CLUSTER_RATIO", 0.72);

    // How many peaks to keep in a +/- 56 Da radius in clustered PRM spectra
    unsigned int rankFiltK = m_params.getValueInt("PRM_RANK_FILTER", 3);

    // Minimum allowable number of PRM spectra within a cluster
    unsigned int minClustSize = m_params.getValueInt("CLUSTER_MIN_SIZE", 0);

    DEBUG_VAR(rankFiltK);

    DEBUG_VAR(minRatioClust);

    // Keep a reference of scan # to lookup index
    map<unsigned int, unsigned int> scanToIdx;
    for (unsigned int i = 0; i < inOutSpecs->size(); i++)
    {
      scanToIdx[(*inOutSpecs)[i].scan] = i;
    }

    // Compute initial set of clusters (pairs or triples)

    // initial clustered spectra
    SpecSet* initialClusters = new SpecSet(inPairedScans.size());

    // keep track of clustered MS/MS indices
    vector<list<unsigned int> > clusterIdxs(initialClusters->size());

    // used to merge PRM spectra
    PairedSpecSet nextCluster;

    for (unsigned int i = 0; i < inPairedScans.size(); i++)
    {
      SpecSet pairedSpecs;
      vector<vector<bool> > pairedRev(0);
      //DEBUG_VAR(i);

      // Put triplet or paired spectra in a container
      for (unsigned int j = 0; j < inPairedScans[i].size(); j++)
      {
        //DEBUG_VAR(j);
        unsigned int specIdx = scanToIdx[inPairedScans[i][j]];
        pairedSpecs.push_back((*inOutSpecs)[specIdx]);
        //DEBUG_VAR((*inOutSpecs)[specIdx].parentMass);
        //DEBUG_VAR((*inOutSpecs)[specIdx].msFragType);
        //DEBUG_VAR((*inOutSpecs)[specIdx].size());

        if (reversedPeaks)
        {
          pairedRev.push_back((*reversedPeaks)[specIdx]);
        }
      }

      // also get the indices of MS/MS spectra
      unsigned int firstScan = inPairedScans[i][0];
      clusterIdxs[i].clear();
      clusterIdxs[i].push_back(scanToIdx[firstScan]);
      for (unsigned int j = 1; j < inPairedScans[i].size(); j++)
      {
        clusterIdxs[i].push_back(scanToIdx[inPairedScans[i][j]]);
      }

      if (reversedPeaks)
      {
        nextCluster.initialize(&pairedSpecs,
                               enforceDaTolerance,
                               &pairedRev,
                               m_params.getValueFloat("TOLERANCE_PEAK"));
      }
      else
      {
        nextCluster.initialize(&pairedSpecs,
                               enforceDaTolerance,
                               0,
                               m_params.getValueFloat("TOLERANCE_PEAK"));
      }

      // merge PRM spectra
      nextCluster.mergePRMs();
      nextCluster.getMergedSpectrum(&(*initialClusters)[i]);

      // rank filter peaks
      (*initialClusters)[i].rankFilterPeaks(rankFiltK);
    }

    DEBUG_VAR(initialClusters->size());
    DEBUG_VAR(inOutSpecs->size());

    // The final set of clusters after hierarchical clustering
    SpecSet finalClusters(initialClusters->size());

    finalClusters = *initialClusters;

    bool foundPair = true;
    PRMAlignment nextPair;
    vector<bool> mergedSpecs(finalClusters.size(), false);
    list<pair<float, TwoValues<unsigned int> > > candidatePairs;
    pair<float, TwoValues<unsigned int> > tempPair;
    unsigned int stage = 1;
    unsigned int numClusters = finalClusters.size();
    Spectrum sortedPrecursors;
    list<int> idxCheck;
    //map<MZRange, list<unsigned int> > precursorMap;

    if (minClustSize > 0)
    {

      while (foundPair)
      {
        DEBUG_MSG("Clustering stage " << stage << " ...");
        stage++;

        foundPair = false;
        candidatePairs.clear();

        DEBUG_MSG("Indexing precursors ...");

        // Use a sorted list (ie a spectrum) to index precursor masses
        sortedPrecursors.resize(finalClusters.size());
        unsigned int tempIdxUse = 0;

        for (unsigned int i = 0; i < finalClusters.size(); i++)
        {
          Spectrum* spec1 = &finalClusters[i];

          if (spec1->size() == 0)
          {
            continue;
          }

          sortedPrecursors[tempIdxUse][0] = spec1->parentMass;
          sortedPrecursors[tempIdxUse][1] = (float)i;
          sortedPrecursors.setTolerance(tempIdxUse, spec1->parentMassTol);
          tempIdxUse++;
        }
        sortedPrecursors.resize(tempIdxUse);

        DEBUG_MSG("Sorting ...");
        sortedPrecursors.sortPeaks();

        DEBUG_MSG("Finding pairs ...");

        // find all pairs of spectra with the same precursors
        for (unsigned int i = 0; i < finalClusters.size(); i++)
        {
          Spectrum* spec1 = &finalClusters[i];
          MZRange nextRange(spec1->parentMass, 0, spec1->parentMassTol);
          if (spec1->size() == 0)
          {
            continue;
          }
          nextPair.setSpec1(spec1);

          // lookup matching precursor masses
          sortedPrecursors.findPeaks(nextRange, &idxCheck);

          for (list<int>::const_iterator idxIt = idxCheck.begin(); idxIt
              != idxCheck.end(); idxIt++)
          {
            unsigned int j = floatToInt(sortedPrecursors[*idxIt][1]);

            if (j <= i)
            {
              continue;
            }

            Spectrum* spec2 = &finalClusters[j];

            if (spec2->size() == 0)
            {
              continue;
            }

            if (!MZRange::EqualWithinRange(spec1->parentMass,
                                           spec2->parentMass,
                                           spec1->parentMassTol
                                               + spec2->parentMassTol))
            {
              continue;
            }

            nextPair.setSpec2(spec2);

            //DEBUG_MSG("Considering " << spec1->parentMass << " and " << spec2->parentMass << " w/ tolerance " << spec1->parentMassTol + spec2->parentMassTol);

            pair<int, pair<float, float> > alignScore =
                nextPair.getShiftScore(0, 0, 1);

            float minR = min(alignScore.second.first, alignScore.second.second);

            //DEBUG_VAR(minR);

            if (minR >= minRatioClust)
            {
              tempPair.first = minR;
              tempPair.second[0] = i;
              tempPair.second[1] = j;
              candidatePairs.push_back(tempPair);
            }
          }
        }

        DEBUG_VAR(candidatePairs.size());

        if (candidatePairs.size() == 0)
        {
          break;
        }

        // order pairs in by decreasing matched score ratio
        candidatePairs.sort(SortPairs());

        for (unsigned int i = 0; i < finalClusters.size(); i++)
        {
          mergedSpecs[i] = false;
        }

        // merge pairs one-by-one, until all are merged
        for (list<pair<float, TwoValues<unsigned int> > >::const_iterator pIt =
            candidatePairs.begin(); pIt != candidatePairs.end(); pIt++)
        {
          unsigned int i = pIt->second[0];
          unsigned int j = pIt->second[1];

          if (mergedSpecs[i] || mergedSpecs[j])
          {
            // this pair has already been merged by transitive pairs
            continue;
          }

          foundPair = true;

          // keep track of clustered indices
          clusterIdxs[i].insert(clusterIdxs[i].end(),
                                clusterIdxs[j].begin(),
                                clusterIdxs[j].end());
          clusterIdxs[j].clear();
          SpecSet pairedSpecs;

          for (list<unsigned int>::iterator cIt = clusterIdxs[i].begin(); cIt
              != clusterIdxs[i].end(); cIt++)
          {
            pairedSpecs.push_back((*inOutSpecs)[*cIt]);
          }

          nextCluster.initialize(&pairedSpecs,
                                 enforceDaTolerance,
                                 0,
                                 m_params.getValueFloat("TOLERANCE_PEAK"));

          nextCluster.mergePRMs();
          nextCluster.getMergedSpectrum(&finalClusters[i]);
          finalClusters[i].rankFilterPeaks(rankFiltK);

          finalClusters[j].resize(0);
          numClusters--;
          mergedSpecs[i] = true;
          mergedSpecs[j] = true;
        }
        DEBUG_VAR(numClusters);
      }
    }

    outClusteredScans.resize(clusterIdxs.size());

    unsigned int idxUse = 0;
    DEBUG_VAR(minClustSize);

    SpecSet outSpecs(inOutSpecs->size());
    DEBUG_VAR(inOutSpecs->size());

    // Compute final set of condensed clusters
    for (unsigned int i = 0; i < clusterIdxs.size(); i++)
    {
      if (clusterIdxs[i].size() < minClustSize || finalClusters[i].size() == 0)
      {
        continue;
      }

      //DEBUG_VAR(idxUse);
      outSpecs[idxUse] = finalClusters[i];

      outClusteredScans[idxUse].resize(0);
      for (list<unsigned int>::iterator cIt = clusterIdxs[i].begin(); cIt
          != clusterIdxs[i].end(); cIt++)
      {
        unsigned int specIdx = *cIt;

        outClusteredScans[idxUse].push_back((*inOutSpecs)[specIdx].scan);
      }

      idxUse++;
    }
    DEBUG_TRACE;

    outSpecs.resize(idxUse);
    inOutSpecs->operator =(outSpecs);
    outClusteredScans.resize(idxUse);

    DEBUG_VAR(inOutSpecs->size());

    delete initialClusters;
  }

  void ExecPrmScoring::checkSpecset(SpecSet &inOutSpecs)
  {
    DEBUG_MSG("Checking sepecset. Size: " << inOutSpecs.size());

    for (int i = 0; i < inOutSpecs.size(); i++)
    {
      //DEBUG_VAR(inOutSpecs[i].size());
      if (inOutSpecs[i].size() <= 8)
      {
        DEBUG_MSG("Spectrum #" << i << " has " << inOutSpecs[i].size() << " peaks. Resizing to 0.");
        inOutSpecs[i].resize(0);
      }
    }

    DEBUG_MSG("Checking sepecset done.");
  }

} // namespace specnets
