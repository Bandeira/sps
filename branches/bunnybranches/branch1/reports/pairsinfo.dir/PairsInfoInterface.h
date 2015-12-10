///////////////////////////////////////////////////////////////////////////////
#ifndef __PAIRS_INFO_INTERFACE_H__
#define __PAIRS_INFO_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"
#include "SpectrumPairSet.h"


///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
class PairsInfoInterface {


  // cluster data object
  SpectrumPairSet pairSet;

  // Input data
  string inputFilename, projectdir;

  // output file and directory
  string outdir, outFileName;


  float m_minScore1, m_minScore2, m_edgeTopKBoth, m_edgeTopKOne;
  bool  m_exists_minScore1, m_exists_minScore2, m_exists_edgeTopKBoth, m_exists_edgeTopKOne;


  int writeOutFile(void);
  
  


 public:


  // Constructors and destructor
  PairsInfoInterface();
  ~PairsInfoInterface();

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
