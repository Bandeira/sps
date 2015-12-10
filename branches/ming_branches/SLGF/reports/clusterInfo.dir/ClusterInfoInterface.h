///////////////////////////////////////////////////////////////////////////////
#ifndef __CLUSTER_INFO_INTERFACE_H__
#define __CLUSTER_INFO_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "ClusterData.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace spsReports;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
class ClusterInfoInterface {



  ClusterData clusterData;

  vector<string> inputFiles;

  vector<string> scanFileNames;

  vector<vector<vector<int> > > scanNumbers;

  string binFiles, projectdir;

  string binFilesNames, inputFilesNames;

  string outdir, outFileName;

  void  loadStringVector(string &inputFilesNames, vector<std::string> &inputFiles);



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
