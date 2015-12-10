////////////////////////////////////////////////////////////////////////////////
#ifndef __CLUSTER_DATA_H__
#define __CLUSTER_DATA_H__
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
#include <list>
#include <map>

#include "SpecSet.h"


////////////////////////////////////////////////////////////////////////////////

#define MS_CLUSTER_MGF_DIR "spectra/out/mgf"
#define MS_CLUSTER_OUT_DIR "spectra/out"
#define MS_CLUSTER_TEMP_DIR "spectra/tmp"
#define MS_CLUSTER_DATA_SUBPATH "spectra/out/clust"
#define MS_CLUSTER_FILE_PREFIX  "clusters_"
#define CLUSTER_DATA_FILENAME   "spectra/clusterData.bin"


using namespace std;
using namespace specnets;

namespace spsReports {

////////////////////////////////////////////////////////////////////////////////
// type ClusterData
//
// Defines the data structure for storing clusters data
//
// Where to load cluster information to
// vector index 1 is cluster #
// Unsigned 1 is spectrum file index
// Unsigned 2 is spectrum index
//
typedef map<unsigned,list<pair<unsigned,unsigned> > > ClusterData_t;

////////////////////////////////////////////////////////////////////////////////
// file clusterData.bin structure
//
// 1 unsinged int -> number of records
// <Records>
//
// Per record:
// 1 unsigned int -> record number
// 1 unsigned int -> number of entries for record
// <Entries>
//
// per entry:
// 1 unsigned int: file index
// 1 unsigned int: spectrum index
//
////////////////////////////////////////////////////////////////////////////////
class ClusterData {

  // Helper method to sort MsCluster file names
  static bool clusterFilenameSort(const string &a, const string &b);
  // Read MsCluster file names, given a directory and a filename prefix
  int  getClusterFileNames(const string &location, const string &prefix, bool sort);

 public:

  // MsCluster clustering data structure
   ClusterData_t data;

   // file names of input MsCluster files
   vector<string> fileNames;

  // Constructors and destructor
  ClusterData()  {fileNames.resize(0);};
  ~ClusterData() {};

  // Load data form MsCluster output files.
  int loadMsClusterData(const string &projectDir);
  // Load data stored in binary format, but load it
  // as if the data appeared twice
  int loadDataDouble(const string &projectDir);
  // Load data stored in binary format
  int loadData(const string &projectDir);
  // Output file in binary (.bin) format
  int saveData(const string &projectDir);


  // General load method
  int load(const string &projectDir, bool saveBinary = true, bool rebuild = false);

  // check if the file is already present. Returns 0 if the file doesn't exist, 1 if it does
  int check(const string &projectDir);

  // cluster count method
  int getClusterCount(void)  {return data.size();};
  // spectra count methods
  int getSpectraCount(void);

  // mapping methods
  int  getConsensusFromInputSpectra(int fileIndex, int inputSpectra);
  list<pair<unsigned,unsigned> > *getInputSpectraFromConsensus(int consensus);

  // output cluster data to cout
  void dump(ostream &sout, bool web);

  // write file in csv format
  int writeCsv(vector<string> *inputFiles, vector<vector<vector<int> > > *scanNumbers, string &outFileName, vector<SpecSet> &specs);

};
////////////////////////////////////////////////////////////////////////////////
};  // namespace sps
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
