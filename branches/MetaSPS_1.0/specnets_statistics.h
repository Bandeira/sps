/**
 @file specnets_statistics.h

 @note
 Copyright 2006, The Regents of the University of California
 All Rights Reserved

 Permission to use, copy, modify and distribute any part of this
 program for educational, research and non-profit purposes, without fee,
 and without a written agreement is hereby granted, provided that the
 above copyright notice, this paragraph and the following three paragraphs
 appear in all copies.

 Those desiring to incorporate this work into commercial
 products or use for commercial purposes should contact the Technology
 Transfer & Intellectual Property Services, University of California,
 San Diego, 9500 Gilman Drive, Mail Code 0910, La Jolla, CA 92093-0910,
 Ph: (858) 534-5815, FAX: (858) 534-7345, E-MAIL:invent@ucsd.edu.

 IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
 FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
 INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN
 IF THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY
 OF SUCH DAMAGE.

 THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY
 OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
 ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF CALIFORNIA MAKES NO
 REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
 EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF
 THE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
 */

#ifndef SPECNETS_STATISTICS_H_
#define SPECNETS_STATISTICS_H_

#include <cstring>
#include <string>
#include <list>
#include <map>
#include <set>
#include <vector>

#include "abruijn.h"
#include "aminoacid.h"
#include "mzrange.h"
#include "spectrum_scoring.h"
#include "spectrum.h"
#include "twovalues.h"
#include "utils.h"

class SpecnetsStatistics
{
public:

  /**
   * Computes the sequence coverage of spectrum_ids on the protein
   * @param protein target protein sequence
   * @param spectrum_ids specnets annotation string for each spectrum
   *   identified
   * @return percent of target protein covered by annotation strings
   */
  static float percentSpectrumCoverage(string& protein,
                                       vector<string>& spectrum_ids);

  /**
   * Computes the sequence coverage of mapped contig sequences on the protein
   * @param protein target protein sequence
   * @param matchma_overlaps for each mapped contig (spectrum idx), this matches
   *   contig peak indices (mass) to protein residue indices (intensity).
   *   This must only include contigs mapped to this protein.
   * @return percent of target protein covered by mapped contig sequences
   */
  static float percentSequencingCoverage(string& protein,
                                         SpecSet& matchma_overlaps);

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
  static float percentContigCoverage(string& protein,
                                     abinfo_t& sps_components,
                                     vector<string>& spectrum_ids);

  /**
   * Finds the protein index that a given contig has the most mapped spectra to.
   * @param contigIdx
   * @param proteins sequence of target proteins
   * @param sps_components contig component information detailing which
   *   spectra were included in each contig
   * @param peptides matched peptides of spectra included in contig
   * @return index of target protein with most mapped spectra. -1 if no spectra map
   */
  static int getContigProtein(int contigIdx,
                              vector<string>& proteins,
                              abinfo_t& sps_components,
                              vector<string>& peptides);

  /**
   * Maps peak indices of annotated spectra to protein indices.
   * @param proteins sequence of target proteins
   * @param spectra PRM spectra of target proteins
   * @param peptides specnets-style peptides matched indexed by corresponding
   *   spectrum indices
   * @param inputIonTypes MS2ScoringModel with b and y ion types specified
   * @param peakTol peak tolerance in Da
   * @param mappedSpectra output data structure of spectrum index -> peak index
   *   -> [protein index, residue index in protein]. protein index is -1 if
   *   spectrum is not annotated, peak index is -1 if peak is not annotated
   *   as b or y.
   * @return
   */
  static void
  mapSpectraToProteins(vector<string>& proteins, SpecSet& spectra, vector<
      string>& peptides, MS2ScoringModel& inputIonTypes, float peakTol, vector<
      vector<TwoValues<int> > >& mappedSpectra);

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
  static bool labelContigPeaks(int contigIdx,
                               abinfo_t& sps_components,
                               SpecSet& star_spectra,
                               vector<vector<TwoValues<int> > >& mappedSpectra,
                               vector<short>& labels);

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
  static int
  countChimericEndpts(int contigIdx,
                      abinfo_t& sps_components,
                      SpecSet& star_spectra,
                      vector<vector<TwoValues<int> > >& mappedSpectra,
                      vector<short>& labels,
                      float peakTol);

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
  static TwoValues<float> contigPercentBY(int contigIdx,
                                          abinfo_t& sps_components,
                                          SpecSet& star_spectra,
                                          vector<string>& spectrum_ids,
                                          MS2ScoringModel& inputIonTypes,
                                          float peakTol,
                                          vector<map<int, string> >& byPeaks);

  /**
   * Calculates the percent of a contig's peaks that are derived from b and y
   *   ions in annotated spectra. A contig peak is labeled as a b or y only if
   *   it is derived from overlapping spectra
   * @param contigIdx index of contig to consider
   * @param sps_components contig component information detailing which
   *   spectra were included in each contig.
   * @param star_spectra PRM spectra assembled into contigs
   * @param spectrum_ids specnets annotation string for each spectrum
   *   identified. must have same indices as spectra in sps_components
   * @param inputIonTypes loaded MS2ScoringModel with b and y types
   * @param peakTol peak mass tolerance
   * @param pyPeaks for every abruijn vertex in contig, this holds a string
   *   detailing if the peak was derived from a b ion ("b"), a y ion ("y"), or
   *   no ion ("").
   * @param byPerc percent of peaks derived from b ions [0], percent of peaks
   *   derived from y ions [1].
   * @return true if contig is derived from non-overlapping annotated spectra,
   *   false otherwise
   */
  static bool
  contigPercentBYCheckChimeric(int contigIdx,
                               abinfo_t& sps_components,
                               SpecSet& star_spectra,
                               vector<string>& spectrum_ids,
                               MS2ScoringModel& inputIonTypes,
                               float peakTol,
                               vector<string>& byPeaks,
                               TwoValues<float>& byPerc);

  /**
   * Computes the sequence length of a contig in AA count
   * @param contigIdx contig index
   * @param matchma_overlaps for each mapped contig (spectrum idx), this matches
   *   contig peak indices (mass) to protein residue indices (intensity).
   * @return AA sequence length covered by contig
   */
  static int sequenceAALength(int contigIdx, SpecSet& matchma_overlaps);

  /**
   * Computes the sequence length of a contig
   * @param contigIdx
   * @param contigs
   * @return sequence length of a contig
   */
  static float sequenceDaLength(int contigIdx, SpecSet& contigs);

  /**
   * Counts the number of spectra used to assemble a contig
   * @param contigIdx
   * @param sps_components contig component information detailing which
   *   spectra were included in each contig.
   * @return number of spectra assembled into contig
   */
  static int countSpectraInContig(int contigIdx, abinfo_t& sps_components);

  /**
   * Counts the number of peptides assembled into a contig
   * @param contigIdx
   * @param sps_components contig component information detailing which
   *   spectra were included in each contig.
   * @param star_spectra contains all spectra assembled into contig
   * @param parent_mass_tol parent mass tolerance in Da
   * @return number of spectra with different parent masses assembled into contig
   */
  static int countPeptidesInContig(int contigIdx,
                                   abinfo_t& sps_components,
                                   SpecSet& star_spectra,
                                   float parent_mass_tol);

  //NON-STATIC METHODS

  SpecSet contigs;
  abinfo_t sps_components;
  SpecSet star_spectra;
  SpecSet matchma_overlaps;
  vector<vector<int> > matchma_prot_idx;
  vector<string> proteins;
  vector<string> spectrum_ids;
  MS2ScoringModel input_ion_types;
  float peak_tol;
  float parent_mass_tol;

  vector<vector<TwoValues<int> > > mappedSpectra;

  /**
   * Default constructor
   */
  SpecnetsStatistics();

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
  SpecnetsStatistics(SpecSet& in_contigs,
                     abinfo_t& in_sps_components,
                     SpecSet& in_star_spectra,
                     SpecSet& in_matchma_overlaps,
                     vector<vector<int> >& in_matchma_prot_idx,
                     vector<string>& in_proteins,
                     vector<string>& in_spectrum_ids,
                     MS2ScoringModel& in_inputIonTypes,
                     float in_peak_tol,
                     float in_parent_mass_tol);

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
  void reset(SpecSet& in_contigs,
             abinfo_t& in_sps_components,
             SpecSet& in_star_spectra,
             SpecSet& in_matchma_overlaps,
             vector<vector<int> >& in_matchma_prot_idx,
             vector<string>& in_proteins,
             vector<string>& in_spectrum_ids,
             MS2ScoringModel& in_inputIonTypes,
             float in_peak_tol,
             float in_parent_mass_tol);

  /**
   * Outputs statistics bundle on all Specnets contigs organized per protein. Includes
   *   any numbers computed by class static methods.
   * @param filename
   * @return true if file was written successfully, false if not
   */
  bool outputCummulativeCSVReport(string& filename);

  /**
   * Calls outputContigCSVReport on every contig with > 0 peaks. The
   *   output filename for each contig is <fileroot><0-padded index>.csv
   * @param fileroot root of filename for each contig report.
   * @return true if every report was written successfully, false if one
   *  failed
   */
  bool outputAllContigsCSVReport(string& fileroot);

  /**
   * Outputs statistics bundle on a Specnets contig. Includes statistics computed by class
   *   static methods that operate on a per-contig basis.
   * @param filename
   * @param contigIdx
   * @return true if file was written successfully, false if not
   */
  bool outputContigCSVReport(string& filename, int contigIdx);

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
  void prepareAbruijnVertexOutput(int contigIdx,
                                  vector<map<int, string> >& by_peaks,
                                  vector<list<int> >& spectrum_idxs,
                                  vector<list<int> >& spectrum_rev,
                                  vector<list<int> >& peak_idxs,
                                  vector<list<float> >& peak_masses,
                                  vector<list<int> >& peak_endpts,
                                  vector<list<int> >& protein_idxs,
                                  vector<list<int> >& residue_idxs,
                                  vector<list<string> >& by_labels);

};

#endif /* SPECNETS_STATISTICS_H_ */
