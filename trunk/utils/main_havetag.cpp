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

using namespace std;
using namespace sps;
using namespace specnets;
using namespace spsReports;

// -------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  Logger::setDefaultLogger(Logger::getLogger(0));

  if (argc != 3) {
    cerr << "Usage: main_havetag spectrum_comp_file psm_file" << endl;
    return -1;
  }

  DEBUG_TRACE;
  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(argv[1], "\t", "", lines)) {
      ERROR_MSG("Unable to open spectrum comparison score file! " << argv[1]);
      return -3;
  }
  DEBUG_VAR(lines.size());

  PeptideSpectrumMatchSet	psmSetSpec;
  if (!psmSetSpec.loadFromFile(argv[2])) {
    ERROR_MSG("Loading Spectra PSM file [" << argv[2] << "]");
    return -3;
  }
  DEBUG_VAR(psmSetSpec.size());

  for (int iResult = 1; iResult < lines.size(); iResult++ ) {
    if (lines[iResult].size() < 8) {
      continue;
    }
    //DEBUG_VAR(iResult);
    //DEBUG_VAR(lines[iResult].size());
    int contigNum;
    sscanf(lines[iResult][1].c_str(), "%d", &contigNum);
    //DEBUG_VAR(contigNum);
    int scanNum;
    sscanf(lines[iResult][2].c_str(), "%d", &scanNum);
    //DEBUG_VAR(scanNum);
    string msgfProtein = lines[iResult][9];
    if (msgfProtein == "---") {
      continue;
    }
    //DEBUG_VAR(msgfProtein );
    cout << contigNum << "\t" << scanNum << "\t" << msgfProtein << "\t";

    for (int iPSM = 0; iPSM < psmSetSpec.size(); iPSM++ ) {
      if (contigNum == psmSetSpec[iPSM]->m_scanNum) {
        if (msgfProtein == psmSetSpec[iPSM]->m_protein) {
          cout << psmSetSpec[iPSM]->m_origAnnotation << "\t";
        }
      }
    }
    cout << endl;
  }


  return 0;
}
