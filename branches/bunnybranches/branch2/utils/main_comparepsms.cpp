//
//  main_specdump - stand alone executable for dumping spectrum data
//
#include "abruijn.h"
#include "ClusterData.h"
#include "CommandLineParser.h"
#include "db_fasta.h"
#include "Logger.h"
#include "ReportTableClusterConsensus.h"
#include "spectrum.h"
#include "SpecSet.h"
#include "SpectrumPairSet.h"
#include "tuple.h"

#include <stdlib.h>

#define MATCH_THRESHOLD 0.8

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;

// -------------------------------------------------------------------------
string stripAnnotation(string & annotation)
{
  string cleanAnnotation; 
  static string aminoAcids("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
  for (int iChar = 0; iChar < annotation.length(); iChar++)
  {
    if (aminoAcids.find_first_of(annotation[iChar]) != string::npos)
    {
      cleanAnnotation += annotation[iChar];
    }
  }
  return cleanAnnotation;
}

// -------------------------------------------------------------------------
float getPercentMatch(string & string1, string & string2)
{
  int i = 0;
  int j = 0;
  int nMatchCount = 0;
  while (i < string1.size() && j < string2.size()) {
    if (string1[i] == string2[j]) {
      nMatchCount++;
      i++;
    }
    j++;
  }

  return  (float)nMatchCount / (float)string1.size();
}


// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 4) {
    cerr << "Usage: main_scorecontig msgf_psm_file spec_psm_file abinfo_file" << endl;
    return -1;
  }

  PeptideSpectrumMatchSet	psmSetMsgf;
  if (!psmSetMsgf.loadFromFile(argv[1])) {
    ERROR_MSG("Loading MSGF PSM file [" << argv[1] << "]");
    return -2;
  }
  PeptideSpectrumMatchSet	psmSetSpec;
  if (!psmSetSpec.loadFromFile(argv[2])) {
    ERROR_MSG("Loading Spectra PSM file [" << argv[2] << "]");
    return -3;
  }

  abinfo_t contigAbinfo;
  if (!Load_abinfo(argv[3], contigAbinfo)) {
    ERROR_MSG("Loading Abinfo file [" << argv[3] << "]");
    return -3;
  }
  //DEBUG_VAR(contigAbinfo.size());

  map<int, int> specToContig;  
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > >::iterator itr = contigAbinfo.begin();
  std::map<unsigned, // contig index
      std::pair<std::pair<vector<int> , vector<int> >, // spectrum index, flipped(1)/not-flipped(0)
          vector<std::pair< // ABruijn vertices
              vector<int> , vector<double> > // Spectrum index, peak mass
          > > >::iterator itr_end = contigAbinfo.end();

  int nUnkonwnContigs = contigAbinfo.size();

  for ( ; itr != itr_end; itr++) {
    int contigIndex = itr->first + 1;  // We prefer 1-based contigs
    //DEBUG_MSG("C = " << contigIndex);
    vector<int> specs = itr->second.first.first;
    for (int i = 0; i < specs.size(); i++) {
      int specIndex = specs[i] + 1;  // We prefer 1-based spectra
      //DEBUG_MSG("  S = " << specIndex);
      specToContig[specIndex] = contigIndex;
    }
  }

  map<int, psmPtr> mapScanMsgf;
  for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ ) {

    psmPtr psmSpec = psmSetMsgf[iPSM];
    mapScanMsgf[psmSpec->m_scanNum] = psmSpec;
  } // for (int iPSM = 0; iPSM < psmSetMsgf.size(); iPSM++ )


  int nNoMatch = 0;
  int nDiffProtein = 0;
  int nExactMatch = 0;
  int nMatchG80 = 0;
  int nDiffPeptide = 0;

  set<int> setNoMsgf;  
  set<int> setExactMatch;  
  set<int> setInexactMatch;  
  set<int> setAssociatedPeptide;
  set<int> setAssociatedBadPeptide;
  set<int> setDiffPeptide;
  map<int, float> mapSpecToMatchPercent;
  map<int, int> mapSpecToMatchType;

  for (int iPSM = 0; iPSM < psmSetSpec.size(); iPSM++ ) {

    const psmPtr & psmSpec = psmSetSpec[iPSM];

    int specIndex = psmSpec->m_scanNum;
    int contigIndex = specToContig[specIndex];

    string cleanAnnotationSpec = stripAnnotation(psmSpec ->m_annotation);

    map<int, psmPtr>::iterator itr = mapScanMsgf.find(specIndex);
    if (itr == mapScanMsgf.end()) {
      setNoMsgf.insert(specIndex);
      continue;
    }

    psmPtr psmMsgf = itr->second;
    string cleanAnnotationMsgf = stripAnnotation(psmMsgf->m_annotation);

    if (cleanAnnotationMsgf == cleanAnnotationSpec) {
      nExactMatch++;
      setExactMatch.insert(specIndex);
      mapSpecToMatchPercent[specIndex] = 1.0;
      mapSpecToMatchType[specIndex] = 1;
      continue;
    }

    //cout << cleanAnnotationSpec << "  " << cleanAnnotationMsgf << endl;
    float percentMatch1 = getPercentMatch(cleanAnnotationSpec, cleanAnnotationMsgf);
    //cout << percentMatch1 << endl;
    float percentMatch2 = getPercentMatch(cleanAnnotationMsgf, cleanAnnotationSpec);
    //cout << percentMatch2 << endl;

    float maxMatch = max(percentMatch1, percentMatch2);
    mapSpecToMatchPercent[specIndex] = maxMatch;

    if (maxMatch >= MATCH_THRESHOLD) {
      //cout << psmSpec->m_origAnnotation << "  " << psmMsgf->m_annotation << endl;
      //cout << cleanAnnotationSpec << "  " << cleanAnnotationMsgf << endl;
      //cout << psmSpec->m_protein << "  " << psmMsgf->m_protein << endl;
      //cout << psmSpec->m_matchOrientation << endl << endl;
      nMatchG80++;

      setInexactMatch.insert(specIndex);
      mapSpecToMatchType[specIndex] = 2;
      continue;
    }

    setDiffPeptide.insert(specIndex);
    mapSpecToMatchType[specIndex] = 3;

  } // for (int iPSM = 0; iPSM < psmSetSpec.size(); iPSM++ )

  cout << endl;
  cout << "Total Msgf scans       : " << psmSetMsgf.size() << endl;
  cout << "Total Spec scans       : " << psmSetSpec.size() << endl;
  cout << "Scan not found in MSGF : " << setNoMsgf.size() << endl;
  cout << "Peptides match exactly : " << setExactMatch.size() << endl;
  cout << "Peptides match > " << MATCH_THRESHOLD << "   : " << setInexactMatch.size() << endl;
  cout << "Mismatched Peptides    : " << setDiffPeptide.size() << endl;
  cout << endl;

  map<int, vector<int> > mapContigVecType;
  for (int iPSM = 0; iPSM < psmSetSpec.size(); iPSM++ ) {
    const psmPtr & psmSpec = psmSetSpec[iPSM];
    int specIndex = psmSpec->m_scanNum;
    int contigIndex = specToContig[specIndex];
    int matchType = mapSpecToMatchType[specIndex];
    if (matchType == 2) matchType = 1;
    if (matchType == 4) matchType = 1;
    if (matchType == 5) matchType = 3;
    mapContigVecType[contigIndex].push_back(matchType);
  }

  int nGoodContig = 0;
  int nBadContig = 0;
  int nHybridContig = 0;
  map<int, int> mapContigType;
  map<int, vector<int> >::iterator itrm = mapContigVecType.begin();
  map<int, vector<int> >::iterator itrm_end = mapContigVecType.end();
  for (; itrm != itrm_end; itrm++) {
    int contigIndex = itrm->first;
    vector<int> & vecType = itrm->second;
    for (int i = 0; i < vecType.size(); i++) {
      if (vecType[i] != 0) {
        if (mapContigType.find(contigIndex) == mapContigType.end()) {
          mapContigType[contigIndex] = vecType[i];
        } else if (mapContigType[contigIndex] != vecType[i]) {
          mapContigType[contigIndex] = 8;
        }
      }
    }
    if (mapContigType[contigIndex] == 1) nGoodContig++;
    if (mapContigType[contigIndex] == 3) nBadContig++;
    if (mapContigType[contigIndex] == 8) nHybridContig++;
  }

  cout << "corr\tcontig\tscan\ttype\tspec_anno\tspec_anno_clean\tmsgf_anno\tmsgf_anno_clean\tspec_prot\tmsgf_prot\tspec_prob\tmsgf_prob\tspec_pvalue\tmsgf_pvalue\torient\tmatch%" << endl;
  for (int iPSM = 0; iPSM < psmSetSpec.size(); iPSM++ ) {

    const psmPtr & psmSpec = psmSetSpec[iPSM];

    int specIndex = psmSpec->m_scanNum;
    int contigIndex = specToContig[specIndex];
    string cleanAnnotationSpec = stripAnnotation(psmSpec->m_origAnnotation);

    cout << mapContigType[contigIndex] << "\t";
    cout << contigIndex << "\t" << specIndex << "\t" << mapSpecToMatchType[specIndex] << "\t";
    cout << psmSpec->m_origAnnotation << "\t" << cleanAnnotationSpec << "\t";

    map<int, psmPtr>::iterator itr = mapScanMsgf.find(specIndex);
    if (itr == mapScanMsgf.end()) {
      cout << "---" << "\t" << "---" << "\t" ;
      cout << psmSpec->m_protein << "\t" << "---" << "\t";
      cout << psmSpec->m_score <<  "\t0.0" << "\t";
      cout << psmSpec->m_pValue <<  "\t0.0" << "\t";
      cout << psmSpec->m_matchOrientation  << "\t0.0" << endl;
      continue;
    }

    psmPtr psmMsgf = itr->second;
    string cleanAnnotationMsgf = stripAnnotation(psmMsgf->m_annotation);

    cout << psmMsgf->m_annotation << "\t" << cleanAnnotationMsgf << "\t";
    cout << psmSpec->m_protein << "\t" << psmMsgf->m_protein << "\t";
    cout << psmSpec->m_score <<  "\t" << -log10(psmMsgf->m_score)  << "\t";
    cout << psmSpec->m_pValue <<  "\t" << psmMsgf->m_pValue  << "\t";
    cout << psmSpec->m_matchOrientation << "\t";
    cout << mapSpecToMatchPercent[specIndex] << endl;
  }

  cout << "Total Contigs          : " << contigAbinfo.size() << endl;
  cout << "Total Contigs 2        : " << mapContigVecType.size() << endl;
  cout << "Correct Contigs        : " << nGoodContig << endl;
  cout << "Bad Contigs            : " << nBadContig << endl;
  cout << "Hybrid Contigs         : " << nHybridContig << endl;
  cout << "Unknown Contigs        : " << contigAbinfo.size() - nGoodContig - nBadContig - nHybridContig << endl;
  cout << endl;

  return 0;
}
