/*
 * ExecReportSPSStats.cpp
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

// Header Includes
#include "ExecReportSPSStats.h"

using namespace std;

namespace specnets
{

  ExecReportSPSStats::ExecReportSPSStats(void) :
    m_contigs(0x0), m_stars(0x0), m_abinfo(0x0), m_overlaps(0x0),
        m_protMatch(0x0), m_fasta(0x0), m_model(0x0), m_specIDs(0x0),
        m_mappedProj(0x0), m_targetProts(0x0), m_proteinSpectra(0x0),
        m_peptides(0x0), ownInput(true), ownOutput(true)
  {
    m_name = "ExecReportSPSStats";
    m_type = "ExecReportSPSStats";
  }

  ExecReportSPSStats::ExecReportSPSStats(const ParameterList & inputParams) :

    ExecBase(inputParams), m_contigs(0x0), m_stars(0x0), m_abinfo(0x0),
        m_overlaps(0x0), m_protMatch(0x0), m_fasta(0x0), m_model(0x0),
        m_specIDs(0x0), m_mappedProj(0x0), m_targetProts(0x0),
        m_proteinSpectra(0x0), m_peptides(0x0), ownInput(true), ownOutput(true)
  {

    m_name = "ExecReportSPSStats";
    m_type = "ExecReportSPSStats";
  }

  ExecReportSPSStats::~ExecReportSPSStats(void)
  {
    if (ownInput) {
      delete m_contigs;
      delete m_stars;
      delete m_abinfo;
      delete m_overlaps;
      delete m_protMatch;
      delete m_fasta;
      delete m_model;
      delete m_specIDs;
      delete m_targetProts;
      delete m_proteinSpectra;
      delete m_peptides;
    }
    if (ownOutput) {
      delete m_mappedProj;
    }
  }

  ExecBase * ExecReportSPSStats::clone(const ParameterList & inputParams) const
  {
    return new ExecReportSPSStats(inputParams);
  }

  bool ExecReportSPSStats::invoke(void)
  {
    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK", 0.5);
    /*
     SpecSet* sps_contigs,
     abinfo_t* sps_components,
     SpecSet* star_spectra,
     SpecSet* matchma_overlaps,
     vector<vector<int> >* matchma_prot_idx,
     vector<string>* _proteins,
     SpecSet* protein_spectra,
     vector<string>* spectrum_ids,
     MS2ScoringModel& input_ion_types,
     float peak_tol
     */
    //DEBUG_MSG("invoking");
    m_mappedProj->mapProt(m_contigs,
                          m_abinfo,
                          m_stars,
                          m_overlaps,
                          m_protMatch,
                          m_targetProts,
                          m_proteinSpectra,
                          m_peptides,
                          *m_model,
                          peakTol);
    //DEBUG_MSG("done invoking");
    return true;
  }

  bool ExecReportSPSStats::loadInputData(void)
  {
    m_contigs = new SpecSet;
    m_stars = new SpecSet;
    m_abinfo = new abinfo_t;
    m_overlaps = new SpecSet;
    m_protMatch = new vector<vector<int> > ;
    m_fasta = new DB_fasta;
    m_model = new MS2ScoringModel;
    m_specIDs = new PeptideSpectrumMatchSet;
    m_mappedProj = new MappedSpecnets;

    if (m_params.exists("INPUT_CONTIGS")) {
      DEBUG_MSG("Loading contig spectra from <" << m_params.getValue("INPUT_CONTIGS") << "> ...");
      if (!m_contigs->LoadSpecSet_pkl(m_params.getValue("INPUT_CONTIGS").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_CONTIGS"));
        return false;
      }
    }
    else if (m_params.exists("INPUT_CONTIGS_PKLBIN")) {
      DEBUG_MSG("Loading contig spectra from <" << m_params.getValue("INPUT_CONTIGS_PKLBIN") << "> ...");
      if (!m_contigs->LoadSpecSet_pklbin(m_params.getValue("INPUT_CONTIGS_PKLBIN").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_CONTIGS_PKLBIN"));
        return false;
      }
    }
    if (m_contigs->size() == 0) {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_CONTIGS or INPUT_CONTIGS_PKLBIN?");
      return false;
    }

    if (m_params.exists("INPUT_STARS")) {
      DEBUG_MSG("Loading star spectra from <" << m_params.getValue("INPUT_STARS") << "> ...");
      if (!m_stars->LoadSpecSet_pkl(m_params.getValue("INPUT_STARS").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_STARS"));
        return false;
      }
    }
    else if (m_params.exists("INPUT_STARS_PKLBIN")) {
      DEBUG_MSG("Loading star spectra from <" << m_params.getValue("INPUT_STARS_PKLBIN") << "> ...");
      if (!m_stars->LoadSpecSet_pklbin(m_params.getValue("INPUT_STARS_PKLBIN").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_STARS_PKLBIN"));
        return false;
      }
    }
    if (m_stars->size() == 0) {
      ERROR_MSG("Input spectra size is 0!, did you specify INPUT_STARS or INPUT_STARS_PKLBIN?");
      return false;
    }

    if (m_params.exists("INPUT_CONTIG_ABINFO")) {
      DEBUG_MSG("Loading abinfo from <" << m_params.getValue("INPUT_CONTIG_ABINFO") << "> ...");
      if (!Load_abinfo(m_params.getValue("INPUT_CONTIG_ABINFO").c_str(),
                       *m_abinfo)) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_CONTIG_ABINFO"));
        return false;
      }
    }
    if (m_abinfo->size() == 0) {
      ERROR_MSG("Input abinfo size is 0!, did you specify INPUT_CONTIG_ABINFO?");
      return false;
    }

    if (m_params.exists("INPUT_MIDX")) {
      DEBUG_MSG("Loading matched peak indicies from <" << m_params.getValue("INPUT_MIDX") << "> ...");
      if (!m_overlaps->LoadSpecSet_pklbin(m_params.getValue("INPUT_MIDX").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_MIDX"));
        return false;
      }
    }
    if (m_overlaps->size() == 0) {
      ERROR_MSG("Input matched peak indicies size is 0!, did you specify INPUT_MIDX?");
      return false;
    }

    if (m_params.exists("INPUT_MP")) {
      DEBUG_MSG("Loading matched protein indicies from <" << m_params.getValue("INPUT_MP") << "> ...");
      if (!Load_binArray(m_params.getValue("INPUT_MP").c_str(), *m_protMatch)) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_MP"));
        return false;
      }
    }
    if (m_protMatch->size() == 0) {
      ERROR_MSG("Input matched protein indicies size is 0!, did you specify INPUT_MP?");
      return false;
    }

    if (m_params.exists("INPUT_FASTA")) {
      DEBUG_MSG("Loading target proteins from <" << m_params.getValue("INPUT_FASTA") << "> ...");
      if (!m_fasta->Load(m_params.getValue("INPUT_FASTA").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_FASTA"));
        return false;
      }
    }
    if (m_fasta->size() == 0) {
      ERROR_MSG("Input fasta proteins size is 0!, did you specify INPUT_FASTA?");
      return false;
    }

    if (m_params.exists("INPUT_ION_TYPES")) {
      DEBUG_MSG("Loading ion types model from <" << m_params.getValue("INPUT_ION_TYPES") << "> ...");
      if (!m_model->LoadModel(m_params.getValue("INPUT_ION_TYPES").c_str())) {
        ERROR_MSG("Could not load " << m_params.getValue("INPUT_ION_TYPES"));
        return false;
      }
    }

    if (m_params.exists("INPUT_SPEC_IDS")) {
      string specIDFormat = m_params.getValue("SPEC_ID_FORMAT");
      bool loadSuccess = true;
      DEBUG_MSG("Loading " << specIDFormat << " results from <" << m_params.getValue("INPUT_SPEC_IDS") << "> ...");
      if (specIDFormat == "InspecT") {
        if (!m_specIDs->loadInspectResultsFile(m_params.getValue("INPUT_SPEC_IDS").c_str())) {
          loadSuccess = false;
        }
      }
      else if (specIDFormat == "MSGFDB") {
        if (!m_specIDs->loadMSGFDBResultsFile(m_params.getValue("INPUT_SPEC_IDS").c_str())) {
          loadSuccess = false;
        }
      }
      else {
        ERROR_MSG("Found unsupported spectrum ID format \'" << specIDFormat << "\'");
        return false;
      }

      if (!loadSuccess) {
        ERROR_MSG("Failed to load " << m_params.getValue("INPUT_SPEC_IDS") << " in " << specIDFormat << " format.");
        return false;
      }
    }

    m_targetProts = new vector<string> (m_fasta->size());
    m_proteinSpectra = new SpecSet(m_fasta->size());
    m_peptides = new vector<string> (m_stars->size());

    set<int> target_proteins;
    list<string> strIdxs;

    if (!splitText(m_params.getValue("TARGET_PROTEINS").c_str(), strIdxs, ";")) {
      return false;
    }

    for (list<string>::iterator strIt = strIdxs.begin(); strIt != strIdxs.end(); strIt++) {
      target_proteins.insert(atoi((*strIt).c_str()));
    }

    FilterFastaProteins(target_proteins, *m_targetProts);

    FilterSpecIds(*m_peptides);

    for (int p = 0; p < m_fasta->size(); p++) {
      (*m_proteinSpectra)[p] = m_fasta->getMassesSpec(p);
    }

    int minContigTag = m_params.getValueInt("MIN_CONTIG_AA_TAG", -1);
    float peakTol = m_params.getValueFloat("TOLERANCE_PEAK", 0.5);

    AAJumps jumps(1);

    int endsChop = m_params.getValueInt("ENDS_CHOP", 0);

    if (endsChop > 0) {
      DEBUG_MSG("Ignoring the first and last " << endsChop << " abruijn vertice(s) of every contig");
      ChopContigEnds(endsChop);
    }

    if (minContigTag > 0) {
      DEBUG_MSG("Ignoring any contigs with an AA tag of " << minContigTag << " residues or less");
      for (int i = 0; i < m_contigs->size(); i++) {
        if ((*m_contigs)[i].size() == 0) {
          continue;
        }
        vector<MZRange> c_masses((*m_contigs)[i].size());
        for (int j = 0; j < c_masses.size(); j++) {
          c_masses[j].set((*m_contigs)[i][j][0], (*m_contigs)[i][j][1], peakTol);
        }
        string tag = getLongestTag(c_masses, jumps);

        //DEBUG_MSG("Contig " << i << " has tag " << tag);
        int tagLen = tag.length();
        if (tagLen < minContigTag) {
          (*m_contigs)[i].resize(0);
        }
      }
    }

    return true;
  }

  bool ExecReportSPSStats::saveInputData(std::vector<std::string> & filenames)
  {
    return false;
  }

  bool ExecReportSPSStats::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_SPS_STATS_FILE")) {
      MappedSPSStatTable outTab(m_mappedProj);
      outTab.prepareTable();
      DEBUG_MSG("Saving cumulative statistics table to <" << m_params.getValue("OUTPUT_SPS_STATS_FILE") << "> ...");
      if (!outTab.printToCSV(m_params.getValue("OUTPUT_SPS_STATS_FILE").c_str())) {
        return false;
      }
    }
    if (m_params.exists("OUTPUT_CONTIG_STATS_FILEROOT")) {
      string my_fileroot = m_params.getValue("OUTPUT_CONTIG_STATS_FILEROOT");
      string::size_type loc = my_fileroot.find(".csv");
      if (loc != string::npos) {
        my_fileroot.erase(loc);
      }
      string filename;
      string max_idx = parseInt(m_contigs->size() - 1);
      int equalizeLength = max_idx.length();

      string place_holder(equalizeLength, '*');

      DEBUG_MSG("Saving contig statistics tables to <" << my_fileroot << place_holder << "> ...");

      MappedContigStatTable outTab(m_mappedProj);
      for (int i = 0; i < m_contigs->size(); i++) {
        if ((*m_contigs)[i].size() == 0)
          continue;

        filename = my_fileroot;
        filename.append(parseInt(i, equalizeLength));
        filename.append(".csv");

        outTab.prepareTable(i);
        if (!outTab.printToCSV(filename.c_str())) {
          return false;
        }
      }
    }
    if (m_params.exists("OUTPUT_SEQACC_POS_TABLE")) {
      OutputTable seqTab;
      vector<pair<int, int> > values;
      m_mappedProj->getGapAccVsPos(-1, values, &(seqTab.values));
      DEBUG_MSG("Saving Gap Accuracy vs. Position Table to <" << m_params.getValue("OUTPUT_SEQACC_POS_TABLE") << "> ...");
      if (!seqTab.printToCSV(m_params.getValue("OUTPUT_SEQACC_POS_TABLE").c_str())) {
        return false;
      }
    }
    if (m_params.exists("OUTPUT_PARAMS_FILECOPY")) {
      DEBUG_MSG("Saving a copy of the input parameters to <" << m_params.getValue("OUTPUT_PARAMS_FILECOPY") << "> ...");
      if (!m_params.writeToFile(m_params.getValue("OUTPUT_PARAMS_FILECOPY"))) {
        return false;
      }
    }
    return true;
  }

  bool ExecReportSPSStats::loadOutputData(void)
  {
    return false;
  }

  vector<ExecBase *> const & ExecReportSPSStats::split(int numSplit)
  {
    m_subModules.resize(0);
    return m_subModules;
  }

  // -------------------------------------------------------------------------
  bool ExecReportSPSStats::merge(void)
  {
    return false;
  }

  // -------------------------------------------------------------------------
  bool ExecReportSPSStats::validateParams(std::string & error)
  {
    m_isValid = false;

    VALIDATE_PARAM_EXIST("INPUT_CONTIGS_PKLBIN");
    VALIDATE_PARAM_EXIST("INPUT_STARS_PKLBIN");
    VALIDATE_PARAM_EXIST("INPUT_CONTIG_ABINFO");
    VALIDATE_PARAM_EXIST("INPUT_MIDX");
    VALIDATE_PARAM_EXIST("INPUT_MP");
    VALIDATE_PARAM_EXIST("INPUT_FASTA");
    VALIDATE_PARAM_EXIST("INPUT_ION_TYPES");
    VALIDATE_PARAM_EXIST("INPUT_SPEC_IDS");
    VALIDATE_PARAM_EXIST("SPEC_ID_FORMAT");
    VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");
    VALIDATE_PARAM_EXIST("TARGET_PROTEINS");

    m_isValid = true;
    return true;
  }

  void ExecReportSPSStats::FilterFastaProteins(set<int>& target_proteins,
                                               vector<string>& put_proteins)
  {
    put_proteins.resize(m_fasta->size());
    for (int i = 0; i < put_proteins.size(); i++) {
      if (target_proteins.count(i) > 0 || target_proteins.size() == 0) {
        put_proteins[i] = m_fasta->getSequence(i);
        //cout << "protein " << i << " has size " << put_proteins[i].length() << "\n";
      }
      else {
        put_proteins[i] = "";
      }
    }
  }

  void ExecReportSPSStats::FilterSpecIds(vector<string>& put_peptides)
  {
    put_peptides.resize(m_stars->size(), "");
    for (int i = 0; i < m_specIDs->size(); i++) {
      if ((*m_specIDs)[i]->m_annotation.length() == 0) {
        continue;
      }
      else {
        int specIdx = (*m_specIDs)[i]->m_scanNum - 1;

        // square brackets indicate gaps, not mods. So any modified residues need to be flanked by round brackets.
        put_peptides[specIdx] = makeBracketsMods((*m_specIDs)[i]->m_annotation);
      }
    }
  }

  void ExecReportSPSStats::ChopContigEnds(int endsChop)
  {
    vector<pair<vector<int> , vector<double> > > *abruijn_verts;
    vector<int>* spectrumIdxs;
    vector<int>* spectrumFlip;

    vector<int> specIdxOut;
    vector<int> specFlipOut;
    vector<pair<vector<int> , vector<double> > > abruijnVertsOut;
    set<int> specIdxUse;
    map<int, int> specFlipUse;

    Spectrum newContig;
    Spectrum newMidx;

    for (int i = 0; i < m_contigs->size(); i++) {
      if ((*m_contigs)[i].size() == 0) { // ignore empty contigs
        continue;
      }
      if ((*m_contigs)[i].size() <= 2 * endsChop) {
        (*m_contigs)[i].resize(0);
        continue;
      }
      abruijn_verts = &(*m_abinfo)[i].second;

      specIdxUse.clear();
      specFlipUse.clear();
      spectrumIdxs = &(*m_abinfo)[i].first.first;
      spectrumFlip = &(*m_abinfo)[i].first.second;

      for (int k = 0; k < spectrumIdxs->size(); k++) {
        specFlipUse[(*spectrumIdxs)[k]] = (*spectrumFlip)[k];
      }

      bool foundFirst = false;
      float firstMass = 0;

      bool foundFirstEnd = false;
      float firstEnd = 0;
      float lastEnd = 0;

      // create new contig, abinfo
      newContig = (*m_contigs)[i];
      newMidx = (*m_overlaps)[i];

      newContig.resize((*m_contigs)[i].size() - 2 * endsChop);
      newMidx.resize(newContig.size());
      abruijnVertsOut.resize(newContig.size());
      for (int j = 0; j < (*m_contigs)[i].size(); j++) {
        spectrumIdxs = &((*abruijn_verts)[j].first);
        if (j < endsChop) { // ignore beginning of contig
          continue;
        }
        else if ((*m_contigs)[i].size() - j <= endsChop) { // ignore end of contig
          // remember the cummulative mass removed from the end in order to update parent mass
          if (!foundFirstEnd) {
            firstEnd = (*m_contigs)[i][j][0];
            foundFirstEnd = true;
          }
          lastEnd = (*m_contigs)[i][j][0];
        }
        else {
          // remember the spectra idxs assembled in remaining abruijn vertices
          for (int k = 0; k < spectrumIdxs->size(); k++) {
            specIdxUse.insert((*spectrumIdxs)[k]);
          }
          // remember the cummulative mass removed from the beginning to update heavier peaks and parent mass
          if (!foundFirst) {
            firstMass = (*m_contigs)[i][j][0];
            foundFirst = true;
          }

          newContig[j - endsChop] = (*m_contigs)[i][j];
          newContig[j - endsChop][0] -= firstMass;

          // copy assembled peak masses
          abruijnVertsOut[j - endsChop] = (*abruijn_verts)[j];
        }
      }
      // subtract removed mass from parent mass
      newContig.parentMass -= (lastEnd - firstEnd + firstMass);

      specIdxOut.resize(specIdxUse.size());
      specFlipOut.resize(specIdxUse.size());
      int idxUse = 0;
      for (set<int>::iterator specIt = specIdxUse.begin(); specIt
          != specIdxUse.end(); specIt++) {
        specIdxOut[idxUse] = *specIt;
        specFlipOut[idxUse] = specFlipUse[*specIt];
        idxUse++;
      }

      idxUse = 0;
      for (int k = 0; k < (*m_overlaps)[i].size(); k++) {
        int contigPeakIdx = floatToInt((*m_overlaps)[i][k][0]);
        if (contigPeakIdx >= endsChop && contigPeakIdx < (*m_contigs)[i].size()
            - endsChop) {
          newMidx[idxUse] = (*m_overlaps)[i][k];
          newMidx[idxUse][0] -= (float) endsChop;
          idxUse++;
        }
      }
      newMidx.resize(idxUse);

      // update contig SpecSet, abinfo, and overlapped indicies
      (*m_contigs)[i] = newContig;
      (*m_abinfo)[i].second = abruijnVertsOut;
      (*m_abinfo)[i].first.first = specIdxOut;
      (*m_abinfo)[i].first.second = specFlipOut;
      (*m_overlaps)[i] = newMidx;
    }
  }
}
