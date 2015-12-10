/*
 * ExecFilterContigPairs.cpp
 *
 *  Created on: Dec 13, 2010
 *      Author: aguthals
 */

// Header Includes
#include "ExecFilterContigPairs.h"

using namespace std;

namespace specnets
{

  ExecFilterContigPairs::ExecFilterContigPairs(void) :
    m_inputContigs(0x0), m_inputPairs(0x0), ownInput(true),
        m_outputAlignments(0x0), ownOutput(true)
  {
    m_name = "ExecFilterContigPairs";
    m_type = "ExecFilterContigPairs";
  }

  ExecFilterContigPairs::ExecFilterContigPairs(const ParameterList & inputParams) :

    ExecBase(inputParams), m_inputContigs(0x0), ownInput(true),
        m_outputAlignments(0x0), ownOutput(true), computeAllPairs(true)
  {

    m_name = "ExecFilterContigPairs";
    m_type = "ExecFilterContigPairs";
  }

  ExecFilterContigPairs::ExecFilterContigPairs(const ParameterList & inputParams,
                                               SpecSet * inputContigs,
                                               SpectrumPairSet * inputPairs,
                                               SpectrumPairSet * outputAlignments) :

    ExecBase(inputParams), m_inputContigs(inputContigs), ownInput(false),
        m_outputAlignments(outputAlignments), ownOutput(false),
        computeAllPairs(false)
  {

    m_name = "ExecFilterContigPairs";
    m_type = "ExecFilterContigPairs";
  }

  ExecFilterContigPairs::ExecFilterContigPairs(const ParameterList & inputParams,
                                               SpecSet * inputContigs,
                                               SpectrumPairSet * inputPairs) :

    ExecBase(inputParams), m_inputContigs(inputContigs),
        m_inputPairs(inputPairs), ownInput(false), m_outputAlignments(0x0),
        ownOutput(true), computeAllPairs(false)
  {

    m_name = "ExecFilterContigPairs";
    m_type = "ExecFilterContigPairs";
  }

  ExecFilterContigPairs::~ExecFilterContigPairs(void)
  {
    if (ownInput) {
      delete m_inputContigs;
      delete m_inputPairs;
    }
    if (ownOutput) {
      delete m_outputAlignments;
    }
  }

  ExecBase * ExecFilterContigPairs::clone(const ParameterList & inputParams) const
  {
    return new ExecFilterContigPairs(inputParams);
  }

  bool ExecFilterContigPairs::invoke(void)
  {
    int num_contigs = m_inputContigs->size();
    int num_pairs = ((num_contigs * num_contigs) - num_contigs) / 2;

    float minScore = m_params.getValueFloat("MIN_METACONTIG_SCORE");

    int minNumMatchedPeaks = m_params.getValueInt("MIN_MATCHED_PEAKS");

    if (m_params.exists("TOLERANCE_PEAK")) {
    	float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");
    	m_inputContigs->setPeakTolerance(peakTol);
    }

    int start_align_idx = m_params.getValueInt("IDX_START", 0);

    int
        end_align_idx =
            (computeAllPairs) ? m_params.getValueInt("IDX_END", num_pairs - 1) : m_params.getValueInt("IDX_END",
                                                                                                      m_inputPairs->size()
                                                                                                          - 1);

    int idx1, idx2;
    int prevIdx1 = -1, prevIdx2 = -1;

    PRMAlignment nextPair;
    //nextPair.peak_tol = peakTol;
    //nextPair.tol_type = 0;

    m_outputAlignments->resize(end_align_idx - start_align_idx + 1);
    int outAlignIdx = 0;

    DEBUG_MSG("Aligning with at least " << minNumMatchedPeaks << " matching peaks and " << minScore << " alignment score");
    DEBUG_MSG("Computing contig/contig alignments for " << end_align_idx - start_align_idx + 1 << " input pairs");

    if (computeAllPairs) {
      int pairIdx = 0;

      for (int i = 0; i < num_contigs; i++) {

        if ((*m_inputContigs)[i].size() < minNumMatchedPeaks)
          continue;

        for (int j = i + 1; j < num_contigs; j++) {

          if ((*m_inputContigs)[j].size() < minNumMatchedPeaks)
            continue;

          if (pairIdx < start_align_idx) {
            ++pairIdx;
            continue;
          }
          else if (pairIdx > end_align_idx)
            break;

          idx1 = i;
          idx2 = j;

          if (idx1 != prevIdx1) {
            nextPair.setSpec1(&(*m_inputContigs)[idx1]);
            prevIdx1 = idx1;
          }

          if (idx2 != prevIdx2) {
            nextPair.setSpec2(&(*m_inputContigs)[idx2]);
            prevIdx2 = idx2;
          }

          bool alignRes = nextPair.computeShiftNoModFR(minNumMatchedPeaks,
                                                       minScore,
                                                       0);

          if (alignRes) {
            nextPair.spec1 = idx1;
            nextPair.spec2 = idx2;

            //DEBUG_MSG("(" << idx1 << "," << idx2 << "): scores = (" << nextPair.score1 << "," << nextPair.score1 << "); shift = " << nextPair.shift1 << "; spec2Rev = " << nextPair.spec2rev);

            (*m_outputAlignments)[outAlignIdx] = (SpectrumPair) nextPair;
            ++outAlignIdx;
          }

          ++pairIdx;
        }
      }
    }
    else {
      for (int i = start_align_idx; i <= end_align_idx; i++) {
        idx1 = (*m_inputPairs)[i].spec1;
        idx2 = (*m_inputPairs)[i].spec2;

        //cout << "Aligning " << idx1 << " and " << idx2 << "\n";

        if (idx1 != prevIdx1) {
          nextPair.setSpec1(&(*m_inputContigs)[idx1]);
          prevIdx1 = idx1;
        }

        if (idx2 != prevIdx2) {
          nextPair.setSpec2(&(*m_inputContigs)[idx2]);
          prevIdx2 = idx2;
        }

        bool alignRes = nextPair.computeShiftNoModFR(minNumMatchedPeaks,
                                                     minScore,
                                                     0);

        if (alignRes) {
          nextPair.spec1 = idx1;
          nextPair.spec2 = idx2;

          //DEBUG_MSG("(" << idx1 << "," << idx2 << "): scores = (" << nextPair.score1 << "," << nextPair.score1 << "); shift = " << nextPair.shift1 << "; spec2Rev = " << nextPair.spec2rev);

          (*m_outputAlignments)[outAlignIdx] = (SpectrumPair) nextPair;
          ++outAlignIdx;
        }
      }
    }

    m_outputAlignments->resize(outAlignIdx);

    DEBUG_MSG("Found " << m_outputAlignments->size() << " optimal alignments");

    return true;
  }

  bool ExecFilterContigPairs::loadInputData(void)
  {
    ownInput = true;
    m_inputContigs = new SpecSet();
    m_inputPairs = new SpectrumPairSet();
    ownOutput = true;
    m_outputAlignments = new SpectrumPairSet();

    //---------------------------------
    // Load spectrum data
    //---------------------------------
    if (m_params.exists("INPUT_CONTIGS")) {
      DEBUG_MSG("Loading contig spectra from <" << m_params.getValue("INPUT_CONTIGS") << "> ...");
      if (!m_inputContigs->LoadSpecSet_pkl(m_params.getValue("INPUT_CONTIGS").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_CONTIGS"));
        return false;
      }
    }
    else if (m_params.exists("INPUT_CONTIGS_PKLBIN")) {
      DEBUG_MSG("Loading contig spectra from <" << m_params.getValue("INPUT_CONTIGS_PKLBIN") << "> ...");
      if (!m_inputContigs->LoadSpecSet_pklbin(m_params.getValue("INPUT_CONTIGS_PKLBIN").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_CONTIGS_PKLBIN"));
        return false;
      }
    }
    if (m_inputContigs->size() == 0) {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_CONTIGS or INPUT_CONTIGS_PKLBIN?");
      return false;
    }

    DEBUG_MSG("Loading specs complete. Num specs [" << m_inputContigs->size() << "]");

    //---------------------------------
    // Load pairs
    //---------------------------------
    if (m_params.exists("INPUT_CONTIG_PAIRS")) {
      DEBUG_MSG("Loading contig pairs to compute from <" << m_params.getValue("INPUT_CONTIG_PAIRS") << "> ...");
      computeAllPairs = false;
      if (!m_inputPairs->loadFromBinaryFile(m_params.getValue("INPUT_CONTIG_PAIRS").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_CONTIG_PAIRS"));
        return false;
      }
      DEBUG_MSG("Loading alignments complete. Num pairs [" << m_inputPairs->size() << "]");
    }
    else {
      DEBUG_MSG("INPUT_CONTIG_PAIRS not specified, using the set of all unique contig pairs");
      computeAllPairs = true;
    }

    return true;
  }

  bool ExecFilterContigPairs::saveInputData(std::vector<std::string> & filenames)
  {
    std::string pairsFilename = getName() + "_pairs.bin";
    m_params.setValue("INPUT_CONTIG_PAIRS", pairsFilename);
    if (!fileExists(pairsFilename) && !computeAllPairs) {
      m_inputPairs->saveToBinaryFile(pairsFilename);
    }

    //SpecSet m_inputContigs; // the input spectra
    std::string spectraFilename = getName() + "_spectra.pklbin";
    m_params.setValue("INPUT_SPECS_PKLBIN", spectraFilename);
    if (!fileExists(spectraFilename)) {
      m_inputContigs->SaveSpecSet_pklbin(spectraFilename.c_str());
    }

    // Have to set up the output files also so the params will be correct on reload
    m_params.setValue("OUTPUT_CONTIG_ALIGNS", getName() + "_pair_aligns.bin");

    std::string paramFilename = getName() + ".params";
    m_params.writeToFile(paramFilename);

    filenames.push_back(paramFilename); // Parameter file MUST be first in vector
    filenames.push_back(pairsFilename);
    filenames.push_back(spectraFilename);

    return true;
  }

  bool ExecFilterContigPairs::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_CONTIG_ALIGNS")) {
      DEBUG_MSG("Saving " << m_outputAlignments->size() << " contig alignments to <" << m_params.getValue(
              "OUTPUT_CONTIG_ALIGNS") << ">");
      return m_outputAlignments->saveToBinaryFile(m_params.getValue("OUTPUT_CONTIG_ALIGNS").c_str());
    }
    return true;
  }

  bool ExecFilterContigPairs::loadOutputData(void)
  {
    if (m_params.exists("OUTPUT_CONTIG_ALIGNS")) {
      DEBUG_MSG("Loading contig alignments from <" << m_params.getValue(
              "OUTPUT_CONTIG_ALIGNS") << ">");
      return m_outputAlignments->loadFromBinaryFile(m_params.getValue("OUTPUT_CONTIG_ALIGNS").c_str());
    }
    return true;
  }

  std::vector<ExecBase *> const & ExecFilterContigPairs::split(int numSplit)
  {

    m_subModules.resize(0);

    if (numSplit < 2) {
      DEBUG_MSG("Number split [" << numSplit << "] must be at least 2");
      return m_subModules;
    }

    if (m_inputContigs->size() == 0) {
      DEBUG_MSG("Must have at least one spectrum");
      return m_subModules;
    }

    m_subModules.resize(numSplit);

    long num_ops = 0;
    int idx1, idx2;
    for (int pairIdx = 0; pairIdx < m_inputPairs->size(); pairIdx++) {
      idx1 = (*m_inputPairs)[pairIdx].spec1;
      idx2 = (*m_inputPairs)[pairIdx].spec2;

      num_ops += (*m_inputContigs)[idx1].size()
          * (*m_inputContigs)[idx2].size();
    }

    long numOpsPerChild = num_ops / ((long) numSplit);
    int globalPairIdx = 0, childIdx = 0;
    int startIdx, endIdx;

    DEBUG_MSG("Splitting into " << numSplit << " children");
    for (int i = 0; i < numSplit; i++) {

      num_ops = 0;
      startIdx = globalPairIdx;
      endIdx = globalPairIdx;

      while (globalPairIdx < m_inputPairs->size() && (num_ops <= numOpsPerChild
          || i == numSplit - 1)) {

        idx1 = (*m_inputPairs)[globalPairIdx].spec1;
        idx2 = (*m_inputPairs)[globalPairIdx].spec2;
        num_ops += (*m_inputContigs)[idx1].size()
            * (*m_inputContigs)[idx2].size();

        ++globalPairIdx;
      }
      endIdx = globalPairIdx - 1;

      ParameterList childParams(m_params);
      childParams.setValue("IDX_START", parseInt(startIdx));
      childParams.setValue("IDX_END", parseInt(endIdx));

      ExecBase * theClone = new ExecFilterContigPairs(childParams,
                                                      m_inputContigs,
                                                      m_inputPairs);

      theClone ->setName(makeName(m_name, i));
      m_subModules[i] = theClone;
    }

    DEBUG_MSG("Splitting success");
    DEBUG_TRACE;
    return m_subModules;
  }

  bool ExecFilterContigPairs::merge(void)
  {
    if (m_subModules.size() == 0) {
      DEBUG_MSG("No children found when merging");
      return false;
    }

    int num_contigs = m_inputContigs->size();
    int num_pairs = ((num_contigs * num_contigs) - num_contigs) / 2;

    m_outputAlignments->resize(num_pairs);
    int alignIdx = 0;

    DEBUG_MSG("Merging " << m_subModules.size() << " children");
    for (int child = 0; child < m_subModules.size(); child++) {
      ExecFilterContigPairs* theChild =
          (ExecFilterContigPairs*) m_subModules[child];
      SpectrumPairSet* childAligns = theChild->m_outputAlignments;

      for (int i = 0; i < childAligns->size(); i++) {
        (*m_outputAlignments)[alignIdx] = (*childAligns)[i];
        ++alignIdx;
      }
    }
    m_outputAlignments->resize(alignIdx);

    DEBUG_MSG("Merging success");
    return true;
  }

  bool ExecFilterContigPairs::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("MIN_METACONTIG_SCORE");
    VALIDATE_PARAM_EXIST("MIN_MATCHED_PEAKS");
    //VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");

    int startBaseIdx;
    if (m_params.exists("IDX_START")) {
      startBaseIdx = max(0, m_params.getValueInt("IDX_START"));
    }
    else {
      startBaseIdx = 0;
    }
    int endBaseIdx;
    if (m_params.exists("IDX_END")) {
      endBaseIdx = max(0, m_params.getValueInt("IDX_END"));
    }
    else {
      endBaseIdx = -1;
    }

    m_isValid = true;
    return true;
  }

}
