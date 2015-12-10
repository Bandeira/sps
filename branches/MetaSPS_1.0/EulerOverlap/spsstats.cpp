/*
 * spsstats.cpp
 *
 *  Created on: Oct 1, 2010
 *      Author: aguthals
 */

#include <cstring>
#include <list>
#include <set>
#include <string>
#include <vector>

#include "SpsStatsHandler.h"
#include "parameter_names.h"

#include "../inputParams.h"
#include "../utils.h"

using namespace std;

int main(int argc, char *argv[], char **envp)
{
  const char* paramFile = (argc == 1) ? "spsmerge.params" : argv[1];

  InputParams params;
  if (!params.readParams(paramFile))
  {
    cerr << "ERROR: Cannot open parameters file " << paramFile << endl;
    cout << endl << SPSSTATS_USAGE;
    return -1;
  }

  vector<const char *> paramStrs(13);

  paramStrs[0] = INPUT_CONTIG_SPECTRA;
  paramStrs[1] = INPUT_STAR_SPECTRA;
  paramStrs[2] = INPUT_CONTIG_ABINFO;
  paramStrs[3] = INPUT_ION_TYPES;
  paramStrs[4] = TOLERANCE_PEAK;
  paramStrs[5] = TOLERANCE_PM;
  paramStrs[6] = INPUT_MATCHED_PEAKS_IDX;
  paramStrs[7] = INPUT_MATCHED_PROTS;
  paramStrs[8] = TARGET_PROTEINS;
  paramStrs[9] = INPUT_FASTA;
  paramStrs[10] = INPUT_INSPECT_RESULTS;
  paramStrs[11] = OUTPUT_STATS_FILE;
  paramStrs[12] = OUTPUT_CONTIG_STATS_FILEROOT;

  if (!params.confirmParams(paramStrs))
  {
    cerr << "ERROR: Parameters file " << paramFile
        << " is incomplete. One of the following is missing: " << paramStrs[0];
    for (int i = 1; i < paramStrs.size(); i++)
      cerr << ", " << paramStrs[i];
    cerr << endl;
    return -1;
  }

  set<int> target_proteins;
  list < string > strIdxs;

  if (!splitText(params.getValue(TARGET_PROTEINS), strIdxs, ":"))
    return -1;

  for (list<string>::iterator strIt = strIdxs.begin(); strIt != strIdxs.end(); strIt++)
    target_proteins.insert(atoi((*strIt).c_str()));

  if (!SpsStatsHandler::OutputCummulativeStatistics(params.getValue(INPUT_CONTIG_SPECTRA),
                                                    params.getValue(INPUT_CONTIG_ABINFO),
                                                    params.getValue(INPUT_STAR_SPECTRA),
                                                    params.getValue(INPUT_MATCHED_PEAKS_IDX),
                                                    params.getValue(INPUT_MATCHED_PROTS),
                                                    params.getValue(INPUT_FASTA),
                                                    params.getValue(INPUT_INSPECT_RESULTS),
                                                    params.getValue(INPUT_ION_TYPES),
                                                    params.getValue(OUTPUT_STATS_FILE),
                                                    params.getValue(OUTPUT_CONTIG_STATS_FILEROOT),
                                                    target_proteins,
                                                    (float)params.getValueDouble(TOLERANCE_PEAK),
                                                    (float)params.getValueDouble(TOLERANCE_PM)))
    return -1;

}
