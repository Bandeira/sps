/*
 * specnets_statistics.cpp
 *
 *  Created on: Sep 23, 2010
 *      Author: aguthals
 */

#include "specnets_statistics.h"

/**
 * Computes the sequence coverage of spectrum_ids on the protein
 * @param protein target protein sequence
 * @param spectrum_ids specnets annotation string for each spectrum
 *   identified
 * @return percent of target protein covered by annotation strings
 */
float SpecnetsStatistics::percentSpectrumCoverage(string& protein, vector<
    string>& spectrum_ids)
{
  vector<bool>* covered_bins = new vector<bool> (protein.length(), false);

  float num_covered_bins = 0;
  string next_id_str, peptideStr;
  string::size_type start_pos;

  for (int i = 0; i < spectrum_ids.size(); i++)
  {
    next_id_str = spectrum_ids[i];

    if (next_id_str.length() == 0)
      continue;

    const char* next_id = next_id_str.c_str();
    char just_letters[strlen(next_id) + 1];

    getPepSeq(next_id, &just_letters[0]);

    peptideStr = &just_letters[0];

    //look for matching sequence in protein
    start_pos = protein.find(peptideStr);

    while (start_pos != string::npos)
    {
      //these residues are now covered if not already
      for (int j = start_pos; j < start_pos + peptideStr.length(); j++)
      {
        if (!(*covered_bins)[j])
        {
          (*covered_bins)[j] = true;
          num_covered_bins += 1.0;
        }
      }
      start_pos = protein.find(peptideStr, start_pos + 1);
    }
  }

  float total_bins = (float)covered_bins->size();

  delete covered_bins;

  return (num_covered_bins / total_bins) * 100.0;
}

/**
 * Computes the sequence coverage of mapped contig sequences on the protein
 * @param protein target protein sequence
 * @param matchma_overlaps for each mapped contig (spectrum idx), this matches
 *   contig peak indicies (mass) to protein residue indicies (intensity).
 *   This must only include contigs mapped to this protein.
 * @return percent of target protein covered by mapped contig sequences
 */
float SpecnetsStatistics::percentSequencingCoverage(string& protein,
                                                    SpecSet& matchma_overlaps)
{
  //all residues are assumed to be uncovered
  vector<bool>* covered_bins = new vector<bool> (protein.length(), false);

  int first_residue_idx, last_residue_idx, temp_idx;

  float num_covered_bins = 0;

  for (int i = 0; i < matchma_overlaps.size(); i++)
  {
    if (matchma_overlaps[i].size() == 0)
      continue;

    first_residue_idx = floatToInt(matchma_overlaps[i].peakList.front()[1]);
    last_residue_idx = floatToInt(matchma_overlaps[i].peakList.back()[1]);
    temp_idx = first_residue_idx;

    first_residue_idx = min(first_residue_idx, last_residue_idx);
    last_residue_idx = max(temp_idx, last_residue_idx);

    //these residues are now covered if not already
    for (int j = first_residue_idx; j <= last_residue_idx; j++)
    {
      if (!(*covered_bins)[j])
      {
        (*covered_bins)[j] = true;
        num_covered_bins += 1.0;
      }
    }
  }

  float total_bins = (float)(covered_bins->size());

  delete covered_bins;

  return (num_covered_bins / total_bins) * 100.0;
}

/**
 * Computes the percentage of assembled protein sequence represented in at
 *   least 3 spectra
 * @param protein target protein sequence
 * @param sps_components contig component information detailing which
 *   spectra were included in each contig. Must only include contigs mapped
 *   to this protein.
 * @param spectrum_ids specnets annotation string for each spectrum
 *   identified. must have same indices as spectra in sps_components
 * @return percentage of assembled protein sequence represented in at
 *   least 1 spectrum
 */
float SpecnetsStatistics::percentContigCoverage(string& protein,
                                                abinfo_t& sps_components,
                                                vector<string>& spectrum_ids)
{
  //all residues are assumed to be uncovered
  vector<int>* covered_bins = new vector<int> (protein.length(), 0);

  float num_covered_bins = 0;
  string next_id_str, peptideStr;
  string::size_type start_pos;
  vector<int>* spectrum_idxs;

  for (abinfo_t::iterator compit = sps_components.begin(); compit
      != sps_components.end(); compit++)
  {
    if (compit->second.second.size() == 0)
      continue;

    spectrum_idxs = &compit->second.first.first;

    for (int i = 0; i < spectrum_idxs->size(); i++)
    {
      next_id_str = spectrum_ids[(*spectrum_idxs)[i]];

      if (next_id_str.length() == 0)
        continue;

      const char* next_id = next_id_str.c_str();
      char just_letters[strlen(next_id) + 1];

      getPepSeq(next_id, &just_letters[0]);

      peptideStr = &just_letters[0];

      //look for matching sequence in protein
      start_pos = protein.find(peptideStr);

      while (start_pos != string::npos)
      {
        //these residues are now covered if not already
        for (int j = start_pos; j < start_pos + peptideStr.length(); j++)
        {
          (*covered_bins)[j]++;
          if ((*covered_bins)[j] == 1)
            num_covered_bins += 1.0;
        }
        start_pos = protein.find(peptideStr, start_pos + 1);
      }
    }
  }

  float total_bins = (float)covered_bins->size();

  delete covered_bins;

  return (num_covered_bins / total_bins) * 100.0;
}

/**
 * Finds the protein index that a given contig has the most mapped spectra to.
 * @param contigIdx
 * @param proteins sequence of target proteins
 * @param sps_components contig component information detailing which
 *   spectra were included in each contig
 * @param peptides matched peptides of spectra included in contig
 * @return index of target protein with most mapped spectra. -1 if no spectra map
 */
int SpecnetsStatistics::getContigProtein(int contigIdx,
                                         vector<string>& proteins,
                                         abinfo_t& sps_components,
                                         vector<string>& peptides)
{

  if (sps_components.count(contigIdx) == 0
      || sps_components[contigIdx].second.size() == 0)
    return -1;

  vector<int>* spectrumIdxs = &sps_components[contigIdx].first.first;
  int specIdx;
  string::size_type start_pos;
  string next_id_str, peptideStr;

  vector<int> protIdxs(proteins.size(), 0);

  for (int i = 0; i < spectrumIdxs->size(); i++)
  {
    specIdx = (*spectrumIdxs)[i];
    next_id_str = peptides[specIdx];

    if (peptides[specIdx].length() == 0)
      continue;

    const char* next_id = next_id_str.c_str();
    char just_letters[strlen(next_id) + 1];

    //get rid of special characters and modifications so we can match to protein
    getPepSeq(next_id, &just_letters[0]);
    peptideStr = &just_letters[0];

    for (int j = 0; j < proteins.size(); j++)
    {
      //look for matching sequence in protein
      start_pos = proteins[j].find(peptideStr);
      if (start_pos == string::npos)
        continue;

      protIdxs[j]++;
    }
  }

  int maxIdxCount = 0;
  int maxIdx = -1;
  int count;

  for (int j = 0; j < protIdxs.size(); j++)
  {
    count = protIdxs[j];
    if (count > maxIdxCount)
    {
      maxIdxCount = count;
      maxIdx = j;
    }
  }

  return maxIdx;
}

/**
 * Maps peak indices of annotated spectra to protein indices.
 * @param proteins sequence of target proteins
 * @param spectra PRM spectra of target proteins
 * @param peptides specnets-style peptides matched indexed by corresponding
 *   spectrum indices. In format
 * @param inputIonTypes MS2ScoringModel with b and y ion types specified
 * @param peakTol peak tolerance in Da
 * @param mappedSpectra output data structure of spectrum index -> peak index
 *   -> [protein index, residue index in protein]. protein index is -1 if
 *   spectrum is not annotated, peak index is -1 if peak is not annotated
 *   as b or y.
 * @return
 */
void SpecnetsStatistics::mapSpectraToProteins(vector<string>& proteins,
                                              SpecSet& spectra,
                                              vector<string>& peptides,
                                              MS2ScoringModel& inputIonTypes,
                                              float peakTol,
                                              vector<vector<TwoValues<int> > >& mappedSpectra)
{
  string next_id_str, peptideStr;
  string::size_type start_pos;
  string ionNamesInclude("b,y");
  Spectrum workingSpec;

  //b/y offset for PRM spectra
  float offset = 0.0 - AAJumps::massHion;

  mappedSpectra.resize(spectra.size());

  for (int i = 0; i < spectra.size(); i++)
  {
    workingSpec = spectra[i];

    mappedSpectra[i].resize(workingSpec.size());
    for (int k = 0; k < workingSpec.size(); k++)
    {
      mappedSpectra[i][k][0] = -1;
      mappedSpectra[i][k][1] = -1;
    }

    if (peptides[i].length() == 0)
      continue;

    next_id_str = peptides[i];
    const char* next_id = next_id_str.c_str();
    char just_letters[strlen(next_id) + 1];

    //get rid of special characters and modifications so we can match to protein
    getPepSeq(next_id, &just_letters[0]);
    peptideStr = &just_letters[0];

    //annotate b and y ions
    workingSpec.annotate(next_id_str,
                         ionNamesInclude,
                         inputIonTypes,
                         offset,
                         offset,
                         peakTol);

    string by;
    int loc;

    for (int j = 0; j < proteins.size(); j++)
    {
      //look for matching sequence in protein
      start_pos = proteins[j].find(peptideStr);
      if (start_pos == string::npos)
        continue;

      for (int k = 0; k < mappedSpectra[i].size(); k++)
      {
        //map to peak index to protein index
        mappedSpectra[i][k][0] = j;

        //map peak index to residue index if b or y

        if (workingSpec.annotation[k].first != NULL)
        {
          by = workingSpec.annotation[k].first->name;
          loc = (int)workingSpec.annotation[k].second;
          loc = (by == "y") ? peptideStr.length() - loc : loc;
          loc--;
        }
        else if (MZRange::EqualWithinRange(workingSpec[k][0], 0, peakTol))
        {
          by = "b";
          loc = -1;
        }
        else if (MZRange::EqualWithinRange(workingSpec[k][0],
                                           workingSpec.parentMass
                                               - AAJumps::massHion,
                                           peakTol))
        {
          by = "y";
          loc = -1;
        }
        else if (MZRange::EqualWithinRange(workingSpec[k][0],
                                           AAJumps::massH2O,
                                           peakTol))
        {
          by = "y";
          loc = peptideStr.length() - 1;
        }
        else if (MZRange::EqualWithinRange(workingSpec[k][0],
                                           workingSpec.parentMass
                                               - AAJumps::massMH,
                                           peakTol))
        {
          by = "b";
          loc = peptideStr.length() - 1;
        }
        else
          continue;

        mappedSpectra[i][k][1] = ((int)start_pos) + 1 + loc;
      }
      break;
    }
  }
}

/**
 * Labels contig peaks as un-annotated (0), incorrect (1), chimeric (2), or correct (3).
 *   un-annotated: All derived peaks come from un-annotated spectra
 *   incorrect: All derived peaks from annotated spectra were not b or y ions
 *   chimeric: All derived peaks from annotated spectra were mapped to different
 *     proteins or different residues
 *   correct: All derived peaks from annotated spectra were mapped to the same protein
 *     and the same residue
 * @param contigIdx
 * @param sps_components contig component information detailing which spectra were
 *   included in each contig peak
 * @param star_spectra spectra that were assembled into contigs
 * @param mappedSpectra output by mapSpectraToProteins. Maps spectrum peaks to protein
 *   residues
 * @param labels output data structure where contig peak index -> label
 * @return true if one peak is chimeric, false otherwise
 */
bool SpecnetsStatistics::labelContigPeaks(int contigIdx,
                                          abinfo_t& sps_components,
                                          SpecSet& star_spectra,
                                          vector<vector<TwoValues<int> > >& mappedSpectra,
                                          vector<short>& labels)
{
  vector < pair<vector<int> , vector<double> > > *abruijn_verts
      = &sps_components[contigIdx].second;
  vector<int>* spectrumIdxs;
  vector<double>* peakMasses;

  int spectrumIdx, peakIdx, proteinIdx, residueIdx;
  float peakMass;

  labels.resize(abruijn_verts->size());
  set<int> mappedProteinIndicies;
  set<int> mappedResidueIndicies;
  bool chimeric = false;

  //check every contig peak
  for (int i = 0; i < abruijn_verts->size(); i++)
  {
    spectrumIdxs = &(*abruijn_verts)[i].first;
    peakMasses = &(*abruijn_verts)[i].second;

    mappedProteinIndicies.clear();
    mappedResidueIndicies.clear();

    //check all derived spectrum peaks
    for (int j = 0; j < spectrumIdxs->size(); j++)
    {
      spectrumIdx = (*spectrumIdxs)[j];
      peakMass = (*peakMasses)[j];
      peakIdx = star_spectra[spectrumIdx].findClosest(peakMass);
      proteinIdx = mappedSpectra[spectrumIdx][peakIdx][0];
      residueIdx = mappedSpectra[spectrumIdx][peakIdx][1];
      if (proteinIdx >= 0)
        mappedProteinIndicies.insert(proteinIdx);
      if (residueIdx >= 0)
        mappedResidueIndicies.insert(residueIdx);
    }

    if (mappedProteinIndicies.size() == 0)
      labels[i] = 0;
    else if (mappedResidueIndicies.size() == 0)
      labels[i] = 1;
    else if (mappedProteinIndicies.size() > 1 || mappedResidueIndicies.size()
        > 1)
      labels[i] = 2;
    else
      labels[i] = 3;

    if (labels[i] == 2)
      chimeric = true;
  }

  return chimeric;
}

/**
 * Counts the number of chimeric abruijn vertices that are derived from spectrum endpoints
 * @param contigIdx contig index
 * @param sps_components contig component information detailing which spectra were
 *   included in each contig peak
 * @param star_spectra spectra that were assembled into contigs
 * @param mappedSpectra mapping spectrum peaks to protein residues
 * @param labels labeling each chimeric contig peak as 2
 * @param peakTol peak mass tolerance in Da
 * @return number of chimeric peaks deriving endpoints in this contig
 */
int SpecnetsStatistics::countChimericEndpts(int contigIdx,
                                            abinfo_t& sps_components,
                                            SpecSet& star_spectra,
                                            vector<vector<TwoValues<int> > >& mappedSpectra,
                                            vector<short>& labels,
                                            float peakTol)
{
  vector < pair<vector<int> , vector<double> > > *abruijn_verts
      = &sps_components[contigIdx].second;

  vector<int>* spectrumIdxs;
  vector<double>* peakMasses;

  int spectrumIdx, peakIdx, proteinIdx, residueIdx;
  float peakMass;

  int count = 0;

  for (int i = 0; i < abruijn_verts->size(); i++)
  {
    if (labels[i] != 2)
      continue;

    spectrumIdxs = &(*abruijn_verts)[i].first;
    peakMasses = &(*abruijn_verts)[i].second;

    //check all derived spectrum peaks
    for (int j = 0; j < spectrumIdxs->size(); j++)
    {
      spectrumIdx = (*spectrumIdxs)[j];
      peakMass = (*peakMasses)[j];
      peakIdx = star_spectra[spectrumIdx].findClosest(peakMass);
      proteinIdx = mappedSpectra[spectrumIdx][peakIdx][0];
      residueIdx = mappedSpectra[spectrumIdx][peakIdx][1];

      if (proteinIdx >= 0 && residueIdx >= 0
          && (peakMass < 57.0 || peakMass
              >= star_spectra[spectrumIdx].parentMass - AAJumps::massHion
                  - peakTol))
      {
        count++;
        break;
      }
    }
  }

  return count;
}

/**
 * Calculates the percent of a contig's peaks that are derived from b and y
 *   ions in annotated spectra
 * @param contigIdx index of contig to consider
 * @param sps_components contig component information detailing which
 *   spectra were included in each contig.
 * @param star_spectra PRM spectra assembled into contigs
 * @param spectrum_ids specnets annotation string for each spectrum
 *   identified. must have same indices as spectra in sps_components
 * @param inputIonTypes loaded MS2ScoringModel with b and y types
 * @param peakTol peak mass tolerance
 * @param pyPeaks for every spectrum peak included in every abruijn
 *   vertex in the contig, this maps abriujn idx -> spectrum idx ->
 *   "b", "y", or ""
 * @return percent of peaks derived from b ions [0], percent of peaks
 *   derived from y ions [1].
 */
TwoValues<float> SpecnetsStatistics::contigPercentBY(int contigIdx,
                                                     abinfo_t& sps_components,
                                                     SpecSet& star_spectra,
                                                     vector<string>& peptides,
                                                     MS2ScoringModel& inputIonTypes,
                                                     float peakTol,
                                                     vector<map<int, string> >& byPeaks)
{
  string next_id_str;
  vector<int>* spectrum_idxs;
  vector<int>::iterator idxIt;
  vector < pair<vector<int> , vector<double> > > *abruijn_verts;

  //make sure contig has peaks
  if (sps_components.count(contigIdx) == 0
      || sps_components[contigIdx].second.size() == 0)
    return TwoValues<float> (0, 0);

  spectrum_idxs = &sps_components[contigIdx].first.first;
  abruijn_verts = &sps_components[contigIdx].second;

  byPeaks.resize(abruijn_verts->size());
  vector<bool> labeled(abruijn_verts->size(), false);

  string ionNamesInclude("b");
  float offset = 0.0 - AAJumps::massHion;
  float num_b = 0, num_y = 0;

  //annotate only b ions
  for (idxIt = spectrum_idxs->begin(); idxIt != spectrum_idxs->end(); idxIt++)
  {
    next_id_str = peptides[*idxIt];
    if (next_id_str.length() == 0)
      continue;
    star_spectra[*idxIt].annotate(next_id_str,
                                  ionNamesInclude,
                                  inputIonTypes,
                                  offset,
                                  offset,
                                  peakTol);
  }

  bool found_b;
  int peakIdx, specIdx;

  //count only b ions
  for (int i = 0; i < abruijn_verts->size(); i++)
  {
    found_b = false;
    for (int j = 0; j < (*abruijn_verts)[i].first.size(); j++)
    {
      specIdx = (*abruijn_verts)[i].first[j];

      //peaks are assumed to be un-annotated
      byPeaks[i][specIdx] = "";

      if (peptides[specIdx].length() == 0)
        continue;

      peakIdx
          = star_spectra[specIdx].findClosest((*abruijn_verts)[i].second[j]);
      if (star_spectra[specIdx].annotation[peakIdx].first != NULL
          || MZRange::EqualWithinRange(star_spectra[specIdx][peakIdx][0],
                                       0,
                                       peakTol)
          || MZRange::EqualWithinRange(star_spectra[specIdx][peakIdx][0],
                                       star_spectra[specIdx].parentMass
                                           - AAJumps::massMH,
                                       peakTol))
      {
        byPeaks[i][specIdx] = "b";
        found_b = true;
      }
    }
    if (found_b)
    {
      num_b += 1.0;
      labeled[i] = true;
    }
  }

  ionNamesInclude = "y";

  //annotate only y ions
  for (idxIt = spectrum_idxs->begin(); idxIt != spectrum_idxs->end(); idxIt++)
  {
    next_id_str = peptides[*idxIt];
    if (next_id_str.length() == 0)
      continue;
    star_spectra[*idxIt].annotate(next_id_str,
                                  ionNamesInclude,
                                  inputIonTypes,
                                  offset,
                                  offset,
                                  peakTol);
  }

  bool found_y;

  //fill in the rest with y ions
  for (int i = 0; i < abruijn_verts->size(); i++)
  {
    found_y = false;
    for (int j = 0; j < (*abruijn_verts)[i].first.size(); j++)
    {
      specIdx = (*abruijn_verts)[i].first[j];

      if (peptides[specIdx].length() == 0)
        continue;

      peakIdx
          = star_spectra[specIdx].findClosest((*abruijn_verts)[i].second[j]);
      if (star_spectra[specIdx].annotation[peakIdx].first != NULL
          || MZRange::EqualWithinRange(star_spectra[specIdx][peakIdx][0],
                                       AAJumps::massH2O,
                                       peakTol)
          || MZRange::EqualWithinRange(star_spectra[specIdx][peakIdx][0],
                                       star_spectra[specIdx].parentMass
                                           - AAJumps::massHion,
                                       peakTol))
      {
        byPeaks[i][specIdx] = "y";
        found_y = true;
      }
    }
    if (found_y && !labeled[i])
    {
      num_y += 1.0;
    }
  }

  float num_peaks = (float)byPeaks.size();
  return TwoValues<float> (100.0 * (num_b / num_peaks), 100.0 * (num_y
      / num_peaks));
}

/**
 * Computes the sequence length of a contig in AA count
 * @param contigIdx contig index
 * @param matchma_overlaps for each mapped contig (spectrum idx), this matches
 *   contig peak indices (mass) to protein residue indices (intensity).
 * @return AA sequence length covered by contig
 */
int SpecnetsStatistics::sequenceAALength(int contigIdx,
                                         SpecSet& matchma_overlaps)
{

  Spectrum* overlap = &matchma_overlaps[contigIdx];

  if (overlap->size() == 0)
    return 0;

  return floatToInt(abs(overlap->peakList.back()[1]
      - overlap->peakList.front()[1]) + 1.0);

}

/**
 * Computes the sequence length of a contig
 * @param contigIdx
 * @param contigs
 * @return sequence length of a contig
 */
float SpecnetsStatistics::sequenceDaLength(int contigIdx, SpecSet& contigs)
{
  return contigs[contigIdx].peakList.back()[0]
      - contigs[contigIdx].peakList.front()[0];
}

/**
 * Counts the number of spectra used to assemble a contig
 * @param contigIdx
 * @param sps_components contig component information detailing which
 *   spectra were included in each contig.
 * @return number of spectra assembled into contig
 */
int SpecnetsStatistics::countSpectraInContig(int contigIdx,
                                             abinfo_t& sps_components)
{
  if (sps_components[contigIdx].second.size() == 0)
    return 0;

  return sps_components[contigIdx].first.first.size();
}

/**
 * Counts the number of peptides assembled into a contig
 * @param contigIdx
 * @param sps_components contig component information detailing which
 *   spectra were included in each contig.
 * @param star_spectra contains all spectra assembled into contig
 * @param parent_mass_tol parent mass tolerance in Da
 * @return number of spectra with different parent masses assembled into contig
 */
int SpecnetsStatistics::countPeptidesInContig(int contigIdx,
                                              abinfo_t& sps_components,
                                              SpecSet& star_spectra,
                                              float parent_mass_tol)
{
  vector<int>* spectrum_idxs;
  map<MZRange, vector<MZRange> > parentMasses;

  if (sps_components[contigIdx].second.size() == 0)
    return 0;

  spectrum_idxs = &sps_components[contigIdx].first.first;

  for (int i = 0; i < spectrum_idxs->size(); i++)
  {
    int idx = (*spectrum_idxs)[i];
    MZRange parentMassRange(star_spectra[idx].parentMass,
                            parent_mass_tol,
                            false);
    MZRange newAverage = MZRange::InsertMZRange(parentMasses, parentMassRange);
  }

  return (int)parentMasses.size();
}

//NON-STATIC METHODS

/**
 * Default constructor
 */
SpecnetsStatistics::SpecnetsStatistics()
{
}

/**
 * Constructor. Just calls reset.
 * @param in_contigs specnets contigs
 * @param in_sps_components contig component information detailing which
 *   spectra and which peaks were included in each contig
 * @param in_star_spectra PRM spectra assembled into contigs
 * @param in_matchma_overlaps matchma output detailing where mapped
 *   contigs overlap on their protein
 * @param in_matchma_prot_idx matchma output detailing which contigs are
 *   mapped to which proteins
 * @param in_proteins sequence of target proteins in same order as matchma
 *   saw them
 * @param in_spectrum_ids sequence of PRM spectra annotations in same order
 *   as in_star_spectra. Must be in specnets format
 * @param in_inputIonTypes loaded MS2ScoringModel with b and y ion types
 *   specified
 * @param in_peak_tol peak tolerance in Da
 * @param in_parent_mass_tol parent mass tolerance in Da
 * @return
 */
SpecnetsStatistics::SpecnetsStatistics(SpecSet& in_contigs,
                                       abinfo_t& in_sps_components,
                                       SpecSet& in_star_spectra,
                                       SpecSet& in_matchma_overlaps,
                                       vector<vector<int> >& in_matchma_prot_idx,
                                       vector<string>& in_proteins,
                                       vector<string>& in_spectrum_ids,
                                       MS2ScoringModel& in_inputIonTypes,
                                       float in_peak_tol,
                                       float in_parent_mass_tol)
{
  reset(in_contigs,
        in_sps_components,
        in_star_spectra,
        in_matchma_overlaps,
        in_matchma_prot_idx,
        in_proteins,
        in_spectrum_ids,
        in_inputIonTypes,
        in_peak_tol,
        in_parent_mass_tol);
}

/**
 * Resets SpecnetsStatistics instance with new data fields. Equivalent
 *   to reconstructing the object.
 * @param in_contigs specnets contigs
 * @param in_sps_components contig component information detailing which
 *   spectra and which peaks were included in each contig
 * @param in_star_spectra PRM spectra assembled into contigs
 * @param in_matchma_overlaps matchma output detailing where mapped
 *   contigs overlap on their protein
 * @param in_matchma_prot_idx matchma output detailing which contigs are
 *   mapped to which proteins
 * @param in_proteins sequence of target proteins in same order as matchma
 *   saw them
 * @param in_spectrum_ids sequence of PRM spectra annotations in same order
 *   as in_star_spectra. Must be in specnets format
 * @param in_inputIonTypes loaded MS2ScoringModel with b and y ion types
 *   specified
 * @param in_peak_tol peak tolerance in Da
 * @param in_parent_mass_tol parent mass tolerance in Da.
 * @return
 */
void SpecnetsStatistics::reset(SpecSet& in_contigs,
                               abinfo_t& in_sps_components,
                               SpecSet& in_star_spectra,
                               SpecSet& in_matchma_overlaps,
                               vector<vector<int> >& in_matchma_prot_idx,
                               vector<string>& in_proteins,
                               vector<string>& in_spectrum_ids,
                               MS2ScoringModel& in_inputIonTypes,
                               float in_peak_tol,
                               float in_parent_mass_tol)
{
  contigs = in_contigs;
  sps_components = in_sps_components;
  star_spectra = in_star_spectra;
  matchma_overlaps = in_matchma_overlaps;
  matchma_prot_idx = in_matchma_prot_idx;
  proteins = in_proteins;
  spectrum_ids = in_spectrum_ids;
  input_ion_types = in_inputIonTypes;
  peak_tol = in_peak_tol;
  parent_mass_tol = in_parent_mass_tol;

  vector<bool> stars_reversed(star_spectra.size(), false);
  for (abinfo_t::iterator seqIt = sps_components.begin(); seqIt
      != sps_components.end(); seqIt++)
  {
    if (seqIt->second.second.size() == 0)
      continue;
    for (int j = 0; j < seqIt->second.first.first.size(); j++)
    {
      stars_reversed[seqIt->second.first.first[j]]
          = seqIt->second.first.second[j] == 1;
    }
  }

  star_spectra.addZPMpeaks(peak_tol, 0, true);

  for (unsigned int i = 0; i < star_spectra.size(); i++)
  {
    if (stars_reversed[i])
      star_spectra[i].reverse(0);
  }

  mapSpectraToProteins(proteins,
                       star_spectra,
                       spectrum_ids,
                       input_ion_types,
                       peak_tol,
                       mappedSpectra);
}

/**
 * Outputs statistics bundle on all Specnets contigs organized per protein. Includes
 *   any numbers computed by class static methods.
 * @param filename
 * @return true if file was written successfully, false if not
 */
bool SpecnetsStatistics::outputCummulativeCSVReport(string& filename)
{
  float tot_residues = 0;
  for (int p = 0; p < proteins.size(); p++)
    tot_residues += (float)proteins[p].length();

  float norm_factor;
  int stats_vector_size = 1;
  vector<int> protIdxs(1, -1);
  for (int i = 0; i < proteins.size(); i++)
  {
    if (proteins[i].length() > 0)
    {
      protIdxs.push_back(i);
      stats_vector_size++;
    }
  }

  vector<int> contig_size(stats_vector_size, 0);
  int tot_contig_size = 0;

  vector<int> contig_peaks(stats_vector_size, 0);
  int tot_contig_peaks = 0;

  vector<float> peaks_annot(stats_vector_size, 0);
  float tot_peaks_annot = 0;

  vector<float> contig_annot_peaks(stats_vector_size, 0);
  float tot_contig_annot_peaks = 0;

  vector<float> spec_cov(stats_vector_size, 0);
  float tot_spec_cov = 0;

  vector<float> contig_cov(stats_vector_size, 0);
  float tot_contig_cov = 0;

  vector<float> seq_cov(stats_vector_size, 0);
  float tot_seq_cov = 0;

  vector<float> assembled_spectra(stats_vector_size, 0);
  float tot_assembled_spectra = 0;

  vector<float> assembled_spectra_annot(stats_vector_size, 0);
  float tot_assembled_spectra_annot = 0;

  vector<float> assembled_peptides(stats_vector_size, 0);
  float tot_assembled_peptides = 0;

  vector<float> assembled_peptides_annot(stats_vector_size, 0);
  float tot_assembled_peptides_annot = 0;

  vector<float> aa_length(stats_vector_size, 0);
  float tot_aa_length = 0;

  vector<float> da_length(stats_vector_size, 0);
  float tot_da_length = 0;

  vector<float> da_length_annot(stats_vector_size, 0);
  float tot_da_length_annot = 0;

  vector<float> perc_b_ions(stats_vector_size, 0);
  float tot_b_ions = 0;

  vector<float> perc_y_ions(stats_vector_size, 0);
  float tot_y_ions = 0;

  vector<float> peaks_chimeric(stats_vector_size, 0);
  float tot_peaks_chimeric = 0;

  vector<float> perc_chimeric_endpts(stats_vector_size, 0);
  float tot_perc_chimeric_endpts = 0;

  vector<float> contigs_chimeric(stats_vector_size, 0);
  float tot_contigs_chimeric = 0;

  vector<float> contigs_chimeric_annot(stats_vector_size, 0);
  float tot_contigs_chimeric_annot = 0;

  vector<float> peaks_correct(stats_vector_size, 0);
  float tot_peaks_correct = 0;

  vector<float> peaks_incorrect(stats_vector_size, 0);
  float tot_peaks_incorrect = 0;

  vector<float> peaks_notannot(stats_vector_size, 0);
  float tot_peaks_notannot = 0;

  vector<float> peaks_notannot_annot(stats_vector_size, 0);
  float tot_peaks_notannot_annot = 0;

  vector<int> contig_matchma(stats_vector_size, 0);
  int tot_contig_matchma = 0;

  vector<int> contig_inspect(stats_vector_size, 0);
  int tot_contig_inspect = 0;

  vector<int> contig_inspect_matchma(stats_vector_size, 0);
  int tot_contig_inspect_matchma = 0;

  vector<TwoValues<float> > long_aa_coverage(stats_vector_size,
                                             TwoValues<float> (0, 0));
  float AA_coverage;
  TwoValues<float> tot_long_aa_coverage(0, 0);

  SpecSet loc_overlaps(matchma_overlaps.size());
  abinfo_t loc_sps_components;

  TwoValues<float> percBY;
  vector < map<int, string> > byPeaks;
  vector<short> labels;
  bool chimeric;
  int p;

  for (int q = 1; q < stats_vector_size; q++)
  {
    p = protIdxs[q];
    string protein = proteins[p];
    norm_factor = ((float)protein.length()) / tot_residues;

    //vector<bool>* test_alloc = new vector<bool>(protein.length());
    //  free(test_alloc);
    spec_cov[q] = percentSpectrumCoverage(protein, spectrum_ids);
    tot_spec_cov += spec_cov[q] * norm_factor;

    for (int s = 0; s < loc_overlaps.size(); s++)
      loc_overlaps[s] = matchma_overlaps[s];

    Copy_abinfo(sps_components, loc_sps_components);

    for (int i = 0; i < contigs.size(); i++)
    {
      if (contigs[i].size() == 0)
      {
        if (loc_sps_components.count(i) > 0)
          loc_sps_components.erase(i);
        loc_overlaps[i].peakList.resize(0);
        continue;
      }

      if (matchma_prot_idx[i][0] == p)
      {
        contig_matchma[q]++;

        AA_coverage = (float)sequenceAALength(i, matchma_overlaps);
        aa_length[q] += AA_coverage;

        if (AA_coverage > long_aa_coverage[q][0])
        {
          long_aa_coverage[q][0] = AA_coverage;
          long_aa_coverage[q][1] = (float)i;
        }
        if (AA_coverage > tot_long_aa_coverage[0])
        {
          tot_long_aa_coverage[0] = AA_coverage;
          tot_long_aa_coverage[1] = (float)i;
        }
      }
      else
      {
        if (loc_sps_components.count(i) > 0)
          loc_sps_components.erase(i);
        loc_overlaps[i].peakList.resize(0);
      }

      int inspect_prot_idx = getContigProtein(i,
                                              proteins,
                                              sps_components,
                                              spectrum_ids);

      if (matchma_prot_idx[i][0] == p && inspect_prot_idx == p)
        contig_inspect_matchma[q]++;

      if (inspect_prot_idx == p)
      {
        contig_annot_peaks[q] += contigs[i].size();

        assembled_spectra_annot[q]
            += (float)countSpectraInContig(i, sps_components);
        assembled_peptides_annot[q]
            += (float)countPeptidesInContig(i,
                                            sps_components,
                                            star_spectra,
                                            parent_mass_tol);

        contig_inspect[q]++;

        da_length_annot[q] += sequenceDaLength(i, contigs);

        chimeric = labelContigPeaks(i,
                                    sps_components,
                                    star_spectra,
                                    mappedSpectra,
                                    labels);

        if (chimeric)
          contigs_chimeric_annot[q] += 1.0;

        for (int j = 0; j < labels.size(); j++)
        {
          if (labels[j] == 0)
            peaks_notannot_annot[q] += 1.0;
        }
      }

      if (matchma_prot_idx[i][0] == p || inspect_prot_idx == p)
      {
        assembled_spectra[q] += (float)countSpectraInContig(i, sps_components);
        assembled_peptides[q] += (float)countPeptidesInContig(i,
                                                              sps_components,
                                                              star_spectra,
                                                              parent_mass_tol);

        contig_size[q]++;
        da_length[q] += sequenceDaLength(i, contigs);

        percBY = contigPercentBY(i,
                                 sps_components,
                                 star_spectra,
                                 spectrum_ids,
                                 input_ion_types,
                                 peak_tol,
                                 byPeaks);

        perc_b_ions[q] += (percBY[0] / 100.0) * (float)byPeaks.size();
        perc_y_ions[q] += (percBY[1] / 100.0) * (float)byPeaks.size();
        contig_peaks[q] += byPeaks.size();

        chimeric = labelContigPeaks(i,
                                    sps_components,
                                    star_spectra,
                                    mappedSpectra,
                                    labels);

        perc_chimeric_endpts[q] += countChimericEndpts(i,
                                                       sps_components,
                                                       star_spectra,
                                                       mappedSpectra,
                                                       labels,
                                                       peak_tol);

        if (chimeric)
          contigs_chimeric[q] += 1.0;

        for (int j = 0; j < labels.size(); j++)
        {
          if (labels[j] == 0)
            peaks_notannot[q] += 1.0;
          else if (labels[j] == 1)
          {
            peaks_incorrect[q] += 1.0;
            peaks_annot[q] += 1.0;
          }
          else if (labels[j] == 2)
          {
            peaks_chimeric[q] += 1.0;
            peaks_annot[q] += 1.0;
          }
          else
          {
            peaks_correct[q] += 1.0;
            peaks_annot[q] += 1.0;
          }
        }
      }
    }

    tot_perc_chimeric_endpts += perc_chimeric_endpts[q];
    perc_chimeric_endpts[q] = 100.0 * perc_chimeric_endpts[q]
        / peaks_chimeric[q];

    tot_assembled_spectra += assembled_spectra[q];
    tot_assembled_spectra_annot += assembled_spectra_annot[q];
    tot_assembled_peptides += assembled_peptides[q];
    tot_assembled_peptides_annot += assembled_peptides_annot[q];
    tot_aa_length += aa_length[q];
    tot_da_length += da_length[q];
    tot_da_length_annot += da_length_annot[q];

    tot_b_ions += perc_b_ions[q];
    tot_y_ions += perc_y_ions[q];
    tot_contigs_chimeric += contigs_chimeric[q];
    tot_contigs_chimeric_annot += contigs_chimeric_annot[q];
    tot_peaks_notannot += peaks_notannot[q];
    tot_contig_annot_peaks += contig_annot_peaks[q];
    tot_peaks_notannot_annot += peaks_notannot_annot[q];
    tot_peaks_incorrect += peaks_incorrect[q];
    tot_peaks_chimeric += peaks_chimeric[q];
    tot_peaks_correct += peaks_correct[q];
    tot_peaks_annot += peaks_annot[q];

    assembled_spectra[q] /= (float)contig_size[q];
    assembled_peptides[q] /= (float)contig_size[q];
    aa_length[q] /= (float)contig_matchma[q];
    da_length[q] /= (float)contig_size[q];
    contigs_chimeric[q] = 100.0 * contigs_chimeric[q] / (float)contig_size[q];
    contigs_chimeric_annot[q] = 100.0 * contigs_chimeric_annot[q]
        / (float)contig_inspect[q];

    perc_b_ions[q] = 100.0 * perc_b_ions[q] / peaks_annot[q];
    perc_y_ions[q] = 100.0 * perc_y_ions[q] / peaks_annot[q];
    peaks_notannot[q] = 100.0 * peaks_notannot[q] / (float)contig_peaks[q];
    peaks_notannot_annot[q] = 100.0 * peaks_notannot_annot[q]
        / (float)contig_annot_peaks[q];
    peaks_incorrect[q] = 100.0 * peaks_incorrect[q] / peaks_annot[q];
    peaks_chimeric[q] = 100.0 * peaks_chimeric[q] / peaks_annot[q];
    peaks_correct[q] = 100.0 * peaks_correct[q] / peaks_annot[q];

    tot_contig_size += contig_size[q];
    tot_contig_matchma += contig_matchma[q];
    tot_contig_inspect += contig_inspect[q];
    tot_contig_inspect_matchma += contig_inspect_matchma[q];
    tot_contig_peaks += contig_peaks[q];

    contig_cov[q] = percentContigCoverage(protein,
                                          loc_sps_components,
                                          spectrum_ids);
    tot_contig_cov += contig_cov[q] * norm_factor;

    seq_cov[q] = percentSequencingCoverage(protein, loc_overlaps);
    tot_seq_cov += seq_cov[q] * norm_factor;
  }

  perc_chimeric_endpts[0] = 100.0 * tot_perc_chimeric_endpts
      / tot_peaks_chimeric;

  contig_cov[0] = tot_contig_cov;
  seq_cov[0] = tot_seq_cov;
  spec_cov[0] = tot_spec_cov;

  contig_size[0] = tot_contig_size;
  contig_matchma[0] = tot_contig_matchma;
  contig_inspect[0] = tot_contig_inspect;
  contig_inspect_matchma[0] = tot_contig_inspect_matchma;

  contig_peaks[0] = tot_contig_peaks;
  contig_annot_peaks[0] = tot_contig_annot_peaks;

  assembled_spectra[0] = tot_assembled_spectra / (float)tot_contig_size;
  assembled_spectra_annot[0] = tot_assembled_spectra_annot
      / (float)tot_contig_inspect;
  assembled_peptides[0] = tot_assembled_peptides / (float)tot_contig_size;
  assembled_peptides_annot[0] = tot_assembled_peptides_annot
      / (float)tot_contig_inspect;
  aa_length[0] = tot_aa_length / (float)tot_contig_matchma;
  da_length[0] = tot_da_length / (float)tot_contig_size;
  da_length_annot[0] = tot_da_length_annot / (float)tot_contig_inspect;
  contigs_chimeric[0] = 100.0 * tot_contigs_chimeric / (float)tot_contig_size;
  contigs_chimeric_annot[0] = 100.0 * tot_contigs_chimeric_annot
      / (float)tot_contig_inspect;

  perc_b_ions[0] = 100.0 * tot_b_ions / tot_peaks_annot;
  perc_y_ions[0] = 100.0 * tot_y_ions / tot_peaks_annot;
  peaks_notannot[0] = 100.0 * tot_peaks_notannot / (float)tot_contig_peaks;
  peaks_notannot_annot[0] = 100.0 * tot_peaks_notannot_annot
      / (float)tot_contig_annot_peaks;
  peaks_incorrect[0] = 100.0 * tot_peaks_incorrect / tot_peaks_annot;
  peaks_chimeric[0] = 100.0 * tot_peaks_chimeric / tot_peaks_annot;
  peaks_correct[0] = 100.0 * tot_peaks_correct / tot_peaks_annot;

  long_aa_coverage[0] = tot_long_aa_coverage;

  int num_contigs_found = 0;
  for (int i = 0; i < contigs.size(); i++)
  {
    if (contigs[i].size() > 0)
      num_contigs_found++;
  }

  printf("Saving cummulative statistics to %s ... ", filename.c_str());
  fflush( stdout);
  FILE* output = fopen(filename.c_str(), "w");
  if (output == NULL)
  {
    cerr << "ERROR: Could not write to file " << filename << endl;
    return false;
  }

  fprintf(output, "Contigs%s%d\n\n", CSV_SEP, num_contigs_found);

  fprintf(output, "Protein Index%s", CSV_SEP);
  fprintf(output, "# Mapped Contigs%s", CSV_SEP);
  fprintf(output, "# Matchma Mapped Contigs%s", CSV_SEP);
  fprintf(output, "# Inspect Mapped Contigs%s", CSV_SEP);
  fprintf(output, "# Matchma & Inspect Mapped Contigs%s", CSV_SEP);

  fprintf(output, "%s Spectrum Coverage%s", "%", CSV_SEP);
  fprintf(output, "%s Contig Coverage%s", "%", CSV_SEP);
  fprintf(output, "%s Sequencing Coverage%s", "%", CSV_SEP);

  fprintf(output, "Spectra Assembled Per Contig%s", CSV_SEP);
  fprintf(output, "Peptides Assembled Per Contig%s", CSV_SEP);
  fprintf(output, "Da Covered Per Contig%s", CSV_SEP);
  fprintf(output, "AA Covered Per Contig%s", CSV_SEP);
  fprintf(output, "Longest AA Covered By Contig%s", CSV_SEP);
  fprintf(output, "Longest AA Contig Index%s", CSV_SEP);

  fprintf(output, "%s Chimeric Contigs%s", "%", CSV_SEP);
  fprintf(output, "%s Annotated Chimeric Contigs%s", "%", CSV_SEP);
  fprintf(output, "%s Annotated Chimeric Peaks%s", "%", CSV_SEP);
  fprintf(output, "%s Chimeric Endpoints%s", "%", CSV_SEP);
  fprintf(output, "%s Annotated Incorrect Peaks%s", "%", CSV_SEP);
  fprintf(output, "%s Un-Annotated Peaks%s", "%", CSV_SEP);
  fprintf(output, "%s Un-Annotated Peaks in Annotated Contigs%s", "%", CSV_SEP);
  fprintf(output, "%s Annotated Correct Peaks%s", "%", CSV_SEP);
  fprintf(output, "%s Annotated Peaks From B ions%s", "%", CSV_SEP);
  fprintf(output, "%s Annotated Peaks From Y ions\n", "%");

  for (int q = 0; q < stats_vector_size; q++)
  {
    if (q == 0)
      fprintf(output, "all");
    else
      fprintf(output, "%d", protIdxs[q]);

    fprintf(output, "%s%d", CSV_SEP, contig_size[q]);
    fprintf(output, "%s%d", CSV_SEP, contig_matchma[q]);
    fprintf(output, "%s%d", CSV_SEP, contig_inspect[q]);
    fprintf(output, "%s%d", CSV_SEP, contig_inspect_matchma[q]);

    fprintf(output, "%s%.1f", CSV_SEP, spec_cov[q]);
    fprintf(output, "%s%.1f", CSV_SEP, contig_cov[q]);
    fprintf(output, "%s%.1f", CSV_SEP, seq_cov[q]);

    fprintf(output, "%s%.1f", CSV_SEP, assembled_spectra[q]);
    fprintf(output, "%s%.1f", CSV_SEP, assembled_peptides[q]);
    fprintf(output, "%s%.1f", CSV_SEP, da_length[q]);
    fprintf(output, "%s%.1f", CSV_SEP, aa_length[q]);
    fprintf(output, "%s%d", CSV_SEP, floatToInt(long_aa_coverage[q][0]));
    fprintf(output, "%s%d", CSV_SEP, floatToInt(long_aa_coverage[q][1]));

    fprintf(output, "%s%.1f", CSV_SEP, contigs_chimeric[q]);
    fprintf(output, "%s%.1f", CSV_SEP, contigs_chimeric_annot[q]);
    fprintf(output, "%s%.1f", CSV_SEP, peaks_chimeric[q]);
    fprintf(output, "%s%.1f", CSV_SEP, perc_chimeric_endpts[q]);
    fprintf(output, "%s%.1f", CSV_SEP, peaks_incorrect[q]);
    fprintf(output, "%s%.1f", CSV_SEP, peaks_notannot[q]);
    fprintf(output, "%s%.1f", CSV_SEP, peaks_notannot_annot[q]);
    fprintf(output, "%s%.1f", CSV_SEP, peaks_correct[q]);
    fprintf(output, "%s%.1f", CSV_SEP, perc_b_ions[q]);
    fprintf(output, "%s%.1f", CSV_SEP, perc_y_ions[q]);

    fprintf(output, "\n");
  }

  fclose(output);

  printf("finished\n");
  fflush(stdout);

  return true;
}

/**
 * Calls outputContigCSVReport on every contig with > 0 peaks. The
 *   output filename for each contig is <fileroot><0-padded index>.csv
 * @param fileroot root of filename for each contig report.
 * @return true if every report was written successfully, false if one
 *  failed
 */
bool SpecnetsStatistics::outputAllContigsCSVReport(string& fileroot)
{
  string my_fileroot = fileroot;
  string::size_type loc = my_fileroot.find(".csv");
  if (loc != string::npos)
  {
    my_fileroot.erase(loc);
  }
  string filename;
  string max_idx = parseInt(contigs.size() - 1);
  int equalizeLength = max_idx.length();

  string place_holder(equalizeLength, '*');
  printf("Saving contig statistics to %s%s.csv ... ", my_fileroot.c_str(), place_holder.c_str());
  fflush( stdout);
  for (int i = 0; i < contigs.size(); i++)
  {
    if (contigs[i].size() == 0)
      continue;

    filename = my_fileroot;
    filename.append(parseInt(i, equalizeLength));
    filename.append(".csv");

    if (!outputContigCSVReport(filename, i))
      return false;
  }
  printf("finished\n");
  fflush(stdout);
  return true;
}

/**
 * Outputs statistics bundle on a Specnets contig. Includes statistics computed by class
 *   static methods that operate on a per-contig basis.
 * @param filename
 * @param contigIdx
 * @return true if file was written successfully, false if not
 */
bool SpecnetsStatistics::outputContigCSVReport(string& filename, int contigIdx)
{
  //border between abruijn vertex rows
  const char* border = "-----";

  const char* labelRef[] =
      { "un-annotated", "incorrect", "chimeric", "correct" };

  //Organize numbers from static methods
  int num_assembled_spec = countSpectraInContig(contigIdx, sps_components);

  int num_assembled_pep = countPeptidesInContig(contigIdx,
                                                sps_components,
                                                star_spectra,
                                                parent_mass_tol);

  int prot_idx = matchma_prot_idx[contigIdx][0];

  int matchma_beg_idx = (prot_idx < 0) ? -1
      : floatToInt(matchma_overlaps[contigIdx].peakList.front()[1]);

  int matchma_end_idx = (prot_idx < 0) ? -1
      : floatToInt(matchma_overlaps[contigIdx].peakList.back()[1]);

  int matchma_cov = (prot_idx < 0) ? -1 : sequenceAALength(contigIdx,
                                                           matchma_overlaps);

  float da_length = sequenceDaLength(contigIdx, contigs);

  int inspect_prot_idx = getContigProtein(contigIdx,
                                          proteins,
                                          sps_components,
                                          spectrum_ids);

  bool has_annotation = inspect_prot_idx >= 0;

  vector<short> labels;
  bool chimeric = labelContigPeaks(contigIdx,
                                   sps_components,
                                   star_spectra,
                                   mappedSpectra,
                                   labels);

  vector < map<int, string> > byPeaks;
  TwoValues<float> percBY = contigPercentBY(contigIdx,
                                            sps_components,
                                            star_spectra,
                                            spectrum_ids,
                                            input_ion_types,
                                            peak_tol,
                                            byPeaks);

  float perc_not_annot = 100.0 - percBY[0] - percBY[1];
  float seq_accuracy = percBY[0] + percBY[1];

  vector<int> label_totals(4, 0);
  int countB = 0, countY = 0, countUnAnnot = 0;
  vector < string > byVerts(byPeaks.size(), "");

  bool foundB, foundY;

  //check if each abruijn vertex includes a b or y ion.
  for (int i = 0; i < labels.size(); i++)
  {
    label_totals[(int)labels[i]]++;
    foundB = false;
    foundY = false;

    for (map<int, string>::iterator labelIt = byPeaks[i].begin(); labelIt
        != byPeaks[i].end(); labelIt++)
    {
      if (labelIt->second == "b")
        foundB = true;
      if (labelIt->second == "y")
        foundY = true;
    }

    if (foundB)
    {
      countB++;
      byVerts[i] = "b";
    }
    else if (foundY)
    {
      countY++;
      byVerts[i] = "y";
    }
    else
      countUnAnnot++;
  }

  vector < list<int> > spectrum_idxs;
  vector < list<int> > spectrum_rev;
  vector < list<int> > peak_idxs;
  vector < list<float> > peak_masses;
  vector < list<int> > peak_endpts;
  vector < list<int> > protein_idxs;
  vector < list<int> > residue_idxs;
  vector < list<string> > by_labels;

  list<int>::iterator intIt;
  list<float>::iterator floatIt;
  list<string>::iterator labelIt;

  prepareAbruijnVertexOutput(contigIdx,
                             byPeaks,
                             spectrum_idxs,
                             spectrum_rev,
                             peak_idxs,
                             peak_masses,
                             peak_endpts,
                             protein_idxs,
                             residue_idxs,
                             by_labels);

  float perc_0 = 100.0 * ((float)label_totals[0] / (float)labels.size());
  float perc_1 = 100.0 * ((float)label_totals[1] / (float)labels.size());
  float perc_2 = 100.0 * ((float)label_totals[2] / (float)labels.size());
  float perc_3 = 100.0 * ((float)label_totals[3] / (float)labels.size());

  //Print statistics
  FILE* output = fopen(filename.c_str(), "w");
  if (output == NULL)
  {
    cerr << "ERROR: Could not write to file " << filename << endl;
    return false;
  }

  //General contig statistics
  fprintf(output, "Contig Index%s%d\n", CSV_SEP, contigIdx);
  fprintf(output, "# of Vertices%s%lu\n", CSV_SEP, (unsigned long)labels.size());
  fprintf(output, "# Assembled Spectra%s%d\n", CSV_SEP, num_assembled_spec);
  fprintf(output, "# Assembled Peptides%s%d\n", CSV_SEP, num_assembled_pep);

  int has_annot = (has_annotation) ? 1 : 0;
  fprintf(output, "# Assembled Annotated Spectra%s%d\n", CSV_SEP, has_annot);

  fprintf(output, "%s Vertices from b or y%s%.1f\n", "%", CSV_SEP, seq_accuracy);
  fprintf(output, "%s Vertices from b%s%.1f\n", "%", CSV_SEP, percBY[0]);
  fprintf(output, "%s Vertices from y%s%.1f\n", "%", CSV_SEP, percBY[1]);

  fprintf(output,
          "%s Vertices not from b or y%s%.1f\n",
          "%",
          CSV_SEP,
          perc_not_annot);

  int chim = (chimeric) ? 1 : 0;
  fprintf(output, "Chimeric%s%d\n", CSV_SEP, chim);

  fprintf(output, "%s Vertices %s %s%.1f\n", "%", labelRef[0], CSV_SEP, perc_0);
  fprintf(output, "%s Vertices %s %s%.1f\n", "%", labelRef[1], CSV_SEP, perc_1);
  fprintf(output, "%s Vertices %s %s%.1f\n", "%", labelRef[2], CSV_SEP, perc_2);
  fprintf(output, "%s Vertices %s %s%.1f\n", "%", labelRef[3], CSV_SEP, perc_3);
  fprintf(output, "Inspect Protein Index%s%d\n", CSV_SEP, inspect_prot_idx);
  fprintf(output, "Matchma Protein Index%s%d\n", CSV_SEP, prot_idx);
  fprintf(output, "Matchma AA Coverage%s%d\n", CSV_SEP, matchma_cov);

  fprintf(output,
          "Matchma Beginning Residue Index%s%d\n",
          CSV_SEP,
          matchma_beg_idx);

  fprintf(output,
          "Matchma Ending Residue Index%s%d\n",
          CSV_SEP,
          matchma_end_idx);

  fprintf(output, "\n");

  fprintf(output,
          "Abruijn Vertex Index%sMass%sIntensity%sLabel%sAnnotation%sAssembled Ions\n",
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP,
          CSV_SEP);

  //Abruijn vertex statistics
  for (int i = 0; i < labels.size(); i++)
  {
    fprintf(output,
            "%d%s%.1f%s%.1f%s%s%s%s",
            i,
            CSV_SEP,
            contigs[contigIdx][i][0],
            CSV_SEP,
            contigs[contigIdx][i][1],
            CSV_SEP,
            labelRef[(int)labels[i]],
            CSV_SEP,
            byVerts[i].c_str());

    fprintf(output, "%sSpectrum Index", CSV_SEP);
    for (intIt = spectrum_idxs[i].begin(); intIt != spectrum_idxs[i].end(); intIt++)
      fprintf(output, "%s%d", CSV_SEP, *intIt);

    fprintf(output,
            "\n%s%s%s%s%sSpectrum Reversed",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
    for (intIt = spectrum_rev[i].begin(); intIt != spectrum_rev[i].end(); intIt++)
      fprintf(output, "%s%d", CSV_SEP, *intIt);

    fprintf(output,
            "\n%s%s%s%s%sPeak Index",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
    for (intIt = peak_idxs[i].begin(); intIt != peak_idxs[i].end(); intIt++)
      fprintf(output, "%s%d", CSV_SEP, *intIt);

    fprintf(output,
            "\n%s%s%s%s%sPeak Mass",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
    for (floatIt = peak_masses[i].begin(); floatIt != peak_masses[i].end(); floatIt++)
      fprintf(output, "%s%.1f", CSV_SEP, *floatIt);

    fprintf(output,
            "\n%s%s%s%s%sPeak Endpoint",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
    for (intIt = peak_endpts[i].begin(); intIt != peak_endpts[i].end(); intIt++)
      fprintf(output, "%s%d", CSV_SEP, *intIt);

    fprintf(output,
            "\n%s%s%s%s%sProtein Index",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
    for (intIt = protein_idxs[i].begin(); intIt != protein_idxs[i].end(); intIt++)
    {
      if (*intIt >= 0)
        fprintf(output, "%s%d", CSV_SEP, *intIt);
      else
        fprintf(output, "%s", CSV_SEP);
    }

    fprintf(output,
            "\n%s%s%s%s%sResidue Index",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
    for (intIt = residue_idxs[i].begin(); intIt != residue_idxs[i].end(); intIt++)
    {
      if (*intIt >= 0)
        fprintf(output, "%s%d", CSV_SEP, *intIt);
      else
        fprintf(output, "%s", CSV_SEP);
    }

    fprintf(output,
            "\n%s%s%s%s%sAnnotation",
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP,
            CSV_SEP);
    for (labelIt = by_labels[i].begin(); labelIt != by_labels[i].end(); labelIt++)
      fprintf(output, "%s%s", CSV_SEP, (*labelIt).c_str());

    fprintf(output,
            "\n%s%s%s%s%s%s%s%s%s%s%s",
            border,
            CSV_SEP,
            border,
            CSV_SEP,
            border,
            CSV_SEP,
            border,
            CSV_SEP,
            border,
            CSV_SEP,
            border);

    for (int j = 0; j < by_labels[i].size(); j++)
      fprintf(output, "%s%s", CSV_SEP, border);

    fprintf(output, "\n");
  }

  fprintf(output, "\nVertex Labels%sTotal\n", CSV_SEP);
  for (int i = 0; i < label_totals.size(); i++)
    fprintf(output, "%s%s%d\n", labelRef[i], CSV_SEP, label_totals[i]);

  fprintf(output, "\nVertex Annotation%sTotal\n", CSV_SEP);
  fprintf(output, "b%s%d\n", CSV_SEP, countB);
  fprintf(output, "y%s%d\n", CSV_SEP, countY);
  fprintf(output, "none%s%d\n", CSV_SEP, countUnAnnot);

  fclose(output);
  return true;
}

/**
 * Helper function for outputContigCSVReport. Prepares abruijn vertex info
 *   into parallel lists for easy printing
 * @param contigIdx
 * @param by_peaks from contigPercentBY
 * @param spectrum_idxs abruijn vert idx -> list of included spectrum indices
 * @param spectrum_rev abruijn vert idx -> list of reversed spectrum idx labels
 * @param peak_idxs abruijn vert idx -> list of included peak indices
 * @param peak_masses abruijn vert idx -> list of included peak masses
 * @param peak_endpts abruijn vert idx -> list of endpoint labels
 * @param protein_idxs abruijn vert idx -> list of matched protein indices
 * @param residue_idxs abruijn vert idx -> list of matched residue indices
 * @param by_labels abruijn vert idx -> list of matched annotations
 * @return
 */
void SpecnetsStatistics::prepareAbruijnVertexOutput(int contigIdx,
                                                    vector<map<int, string> >& by_peaks,
                                                    vector<list<int> >& spectrum_idxs,
                                                    vector<list<int> >& spectrum_rev,
                                                    vector<list<int> >& peak_idxs,
                                                    vector<list<float> >& peak_masses,
                                                    vector<list<int> >& peak_endpts,
                                                    vector<list<int> >& protein_idxs,
                                                    vector<list<int> >& residue_idxs,
                                                    vector<list<string> >& by_labels)
{
  vector < pair<vector<int> , vector<double> > > *abruijn_verts
      = &sps_components[contigIdx].second;

  vector<int>* spectrumIdxs;
  vector<double>* peakMasses;

  map<int, int> specRev;
  vector<int>* Idxs = &sps_components[contigIdx].first.first;
  vector<int>* Rev = &sps_components[contigIdx].first.second;
  for (int i = 0; i < Idxs->size(); i++)
  {
    specRev[(*Idxs)[i]] = (*Rev)[i];
  }

  int spectrum_idx, peak_idx, protein_idx, residue_idx, peak_endpt;
  float peak_mass;
  string by_label;

  spectrum_idxs.resize(abruijn_verts->size());
  spectrum_rev.resize(abruijn_verts->size());
  peak_idxs.resize(abruijn_verts->size());
  peak_masses.resize(abruijn_verts->size());
  peak_endpts.resize(abruijn_verts->size());
  protein_idxs.resize(abruijn_verts->size());
  residue_idxs.resize(abruijn_verts->size());
  by_labels.resize(abruijn_verts->size());

  for (int i = 0; i < abruijn_verts->size(); i++)
  {
    spectrumIdxs = &(*abruijn_verts)[i].first;
    peakMasses = &(*abruijn_verts)[i].second;

    spectrum_idxs[i].clear();
    spectrum_rev[i].clear();
    peak_idxs[i].clear();
    peak_masses[i].clear();
    peak_endpts[i].clear();
    protein_idxs[i].clear();
    residue_idxs[i].clear();
    by_labels[i].clear();

    for (int j = 0; j < spectrumIdxs->size(); j++)
    {
      spectrum_idx = (*spectrumIdxs)[j];
      peak_idx = star_spectra[spectrum_idx].findClosest((*peakMasses)[j]);
      protein_idx = mappedSpectra[spectrum_idx][peak_idx][0];
      residue_idx = mappedSpectra[spectrum_idx][peak_idx][1];
      by_label = by_peaks[i][spectrum_idx];
      peak_mass = star_spectra[spectrum_idx][peak_idx][0];
      peak_endpt = (peak_mass < 57 || peak_mass
          >= star_spectra[spectrum_idx].parentMass - AAJumps::massHion
              - peak_tol) ? 1 : 0;

      spectrum_idxs[i].push_back(spectrum_idx);
      spectrum_rev[i].push_back(specRev[spectrum_idx]);
      peak_idxs[i].push_back(peak_idx);
      peak_masses[i].push_back(peak_mass);
      peak_endpts[i].push_back(peak_endpt);
      protein_idxs[i].push_back(protein_idx);
      residue_idxs[i].push_back(residue_idx);
      by_labels[i].push_back(by_label);
    }
  }
}

