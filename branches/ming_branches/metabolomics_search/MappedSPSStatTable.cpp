/*
 * MappedSpecnetsStatTable.cpp
 *
 *  Created on: Mar 4, 2011
 *      Author: aguthals
 */

#include "MappedSPSStatTable.h"

namespace specnets
{

  MappedSPSStatTable::MappedSPSStatTable(MappedSpecnets* _mapped_sps_proj) :
    OutputTable()
  {
    mapped_sps_proj = _mapped_sps_proj;
  }

  /**
   * Prepares output table with all necessary statistics
   * @return
   */
  void MappedSPSStatTable::prepareTable()
  {
    pair<string, bool> deflt("", false);
    int numRows = 31;

    values.resize(numRows);

    values[0].resize(2, deflt);
    values[0][0].first = "Contigs";
    values[0][0].second = true;
    values[0][1].first = parseInt(mapped_sps_proj->getNumContigs());

    values[1].resize(2, deflt);
    values[1][0].first = "Spectra";
    values[1][0].second = true;
    values[1][1].first = parseInt(mapped_sps_proj->getNumSpectra());

    values[2].resize(2, deflt);
    values[2][0].first = "Identified Spectra";
    values[2][0].second = true;
    values[2][1].first = parseInt(mapped_sps_proj->getNumSpecIdent());

    values[3].resize(0);

    for (int i = 4; i < numRows; i++) {
      values[i].resize(mapped_sps_proj->proteins->size() + 2, deflt);
      //cout << "resizing " << i << " to " << mapped_sps_proj->proteins->size() + 2 << ", actual size is " << values[i].size() << "\n"; cout.flush();
    }

    values[4][0].first = "Protein Index";
    values[4][0].second = true;
    values[4][1].first = "all";
    values[4][1].second = true;

    values[5][0].first = "Identified or Mapped Contigs";
    values[5][0].second = true;

    values[6][0].first = "Identified Contigs";
    values[6][0].second = true;

    values[7][0].first = "Mapped Contigs";
    values[7][0].second = true;

    values[8][0].first = "Identified and Mapped Contigs";
    values[8][0].second = true;

    values[9][0].first = "Identified Spectra";
    values[9][0].second = true;

    values[10][0].first = "Assembled Spectra";
    values[10][0].second = true;

    values[11][0].first = "Assembled Identified Spectra (%)";
    values[11][0].second = true;

    values[12][0].first = "Correctly Mapped Verts (%)";
    values[12][0].second = true;

    values[13][0].first = "Spectrum Coverage (%)";
    values[13][0].second = true;

    values[14][0].first = "Sequencing Coverage (%)";
    values[14][0].second = true;

    values[15][0].first = "Sequencing Coverage Redundancy";
    values[15][0].second = true;

    values[16][0].first = "Spectra Assembled Per Contig";
    values[16][0].second = true;

    values[17][0].first = "Peptides Assembled Per Contig";
    values[17][0].second = true;

    values[18][0].first = "Da Length Per Contig";
    values[18][0].second = true;

    values[19][0].first = "AA Covered Per Contig";
    values[19][0].second = true;

    values[20][0].first = "Longest AA Region Covered by Contig";
    values[20][0].second = true;

    values[21][0].first = "Contig Covering Longest AA Region";
    values[21][0].second = true;

    values[22][0].first = "Annotated Vertices Correct (%)";
    values[22][0].second = true;

    values[23][0].first = "Annotated Vertices Chimeric (%)";
    values[23][0].second = true;

    values[24][0].first = "Annotated Vertices Incorrect (%)";
    values[24][0].second = true;

    values[25][0].first = "Un-annotated Vertices (%)";
    values[25][0].second = true;

    values[26][0].first = "B Vertices (%)";
    values[26][0].second = true;

    values[27][0].first = "Y Vertices (%)";
    values[27][0].second = true;

    values[28][0].first = "Annotated Gaps Correct (%)";
    values[28][0].second = true;

    values[29][0].first = "Annotated Gaps Incorrect (%)";
    values[29][0].second = true;

    values[30][0].first = "Un-annotated Gaps (%)";
    values[30][0].second = true;

    int numProts = mapped_sps_proj->proteins->size();
    //cout << "proteins size=" << numProts << "\n"; cout.flush();

    for (int prot_idx = -1; prot_idx < numProts; prot_idx++) {

      int tp = prot_idx + 2;
      //cout << "on " << prot_idx << "\n";
      if (prot_idx >= 0 && (*mapped_sps_proj->proteins)[prot_idx].length() == 0) {
        //        cout << "skipping " << prot_idx << "\n";
        continue;
      }

      if (prot_idx >= 0) {
        values[4][tp].first = parseInt(prot_idx);
        //        cout << "setting " << 4 << ", " << tp << " to " << parseInt(prot_idx) << ", actually got " << values[4][tp].first << "\n";
      }
      values[5][tp].first
          = parseInt(mapped_sps_proj->getNumContigsMapped(prot_idx));
      values[6][tp].first
          = parseInt(mapped_sps_proj->getNumContigsSpecMapped(prot_idx));
      values[7][tp].first
          = parseInt(mapped_sps_proj->getNumContigsVertMapped(prot_idx));
      values[8][tp].first
          = parseInt(mapped_sps_proj->getNumContigsVertSpecMapped(prot_idx));
      values[9][tp].first
          = parseInt(mapped_sps_proj->getNumSpecMapped(prot_idx));
      values[10][tp].first
          = parseInt(mapped_sps_proj->getNumAssembledSpec(prot_idx));
      values[11][tp].first
          = parseFloat(mapped_sps_proj->getPercAssemSpecIdent(prot_idx), 1);
      values[12][tp].first
          = parseFloat(mapped_sps_proj->getMatchmaAccuracy(prot_idx), 1);
      values[13][tp].first
          = parseFloat(mapped_sps_proj->getPercSpecCov(prot_idx), 1);
      values[14][tp].first
          = parseFloat(mapped_sps_proj->getPercSeqCov(prot_idx), 1);
      values[15][tp].first
          = parseFloat(mapped_sps_proj->getCovRedundancy(prot_idx), 1);
      values[16][tp].first
          = parseFloat(mapped_sps_proj->getSpecPerContig(prot_idx), 1);
      values[17][tp].first
          = parseFloat(mapped_sps_proj->getPepPerContig(prot_idx), 1);
      values[18][tp].first
          = parseFloat(mapped_sps_proj->getDaLengthPerContig(prot_idx), 1);
      values[19][tp].first
          = parseFloat(mapped_sps_proj->getAALengthPerContig(prot_idx), 1);
      pair<int, int> aaRes = mapped_sps_proj->getLongestAAContig(prot_idx);
      values[20][tp].first = parseInt(aaRes.first);
      values[21][tp].first = parseInt(aaRes.second);
      values[22][tp].first = parseFloat(mapped_sps_proj->getPercVerts(prot_idx,
                                                                      3), 1);
      values[23][tp].first = parseFloat(mapped_sps_proj->getPercVerts(prot_idx,
                                                                      2), 1);
      values[24][tp].first = parseFloat(mapped_sps_proj->getPercVerts(prot_idx,
                                                                      1), 1);
      values[25][tp].first = parseFloat(mapped_sps_proj->getPercVerts(prot_idx,
                                                                      0), 1);
      pair<float, float> byRes = mapped_sps_proj->getPercBYVerts(prot_idx);
      values[26][tp].first = parseFloat(byRes.first, 1);
      values[27][tp].first = parseFloat(byRes.second, 1);
      values[28][tp].first = parseFloat(mapped_sps_proj->getPercGaps(prot_idx,
                                                                     2), 1);
      values[29][tp].first = parseFloat(mapped_sps_proj->getPercGaps(prot_idx,
                                                                     1), 1);
      values[30][tp].first = parseFloat(mapped_sps_proj->getPercGaps(prot_idx,
                                                                     0), 1);
    }
  }
}
