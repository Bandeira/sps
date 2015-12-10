///////////////////////////////////////////////////////////////////////////////
#ifndef __CLUSTER_INFO_INTERFACE_H__
#define __CLUSTER_INFO_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "ClusterData.h"
#include "SpecSet.h"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace spsReports;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
class ClusterInfoInterface {


  // cluster data object
  ClusterData clusterData;

  // input files filenames (pklbin)
  vector<string> inputFiles;

  // original files filenames
  vector<string> originalFiles;

  // scan numbers filenames (.bin files)
  vector<string> scanFileNames;

  // the scan numbers
  vector<vector<vector<int> > > scanNumbers;

  //
  string binFiles, projectdir;

  // input filenames of filename containers
  string binFilesNames, inputFilesNames, originalFilesNames;

  // output file and directory
  string outdir, outFileName;

  // Input data (specsets)
  vector<SpecSet> inputData;


  int  loadStringVector(string &inputFilesNames, vector<std::string> &inputFiles);



 public:


  // Constructors and destructor
  ClusterInfoInterface();
  ~ClusterInfoInterface();

  // Option parsing
  int processOptions(int argc, char **argv);

  // Output help
  int help(ostream &);

  // Output version information
  int version(ostream &);

  // Output error messages
  int error(const string &);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
