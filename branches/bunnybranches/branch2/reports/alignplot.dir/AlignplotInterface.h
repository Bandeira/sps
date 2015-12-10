///////////////////////////////////////////////////////////////////////////////
#ifndef __ALIGNPLOT_INTERFACE_H__
#define __ALIGNPLOT_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "spectrum.h"
#include "SpecSet.h"
#include "PlotContig.h"
///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
class AlignplotInterface {

  // abruijn graph
  abinfo_t  abinfo;
  // star spectra
  SpecSet   star;
  // seqs file
  SpecSet   seqs;


  // Spectrum draw object -- used to draw image
  PlotContig  plotContig;

  // Spectrum index.
  int m_spectrumIndex;

  // spectrum scan
  string m_spectrumScan;


 public:


  // Constructors and destructor
  AlignplotInterface();
  ~AlignplotInterface();

  // Option parsing
  int processOptions(int argc, char **argv);

  // Execution based on options
  int processOptions(ParameterList & commandLineParams);


  // Output help
  int help(ostream &);

  // Output version information
  int version(ostream &);

  // Output error messages
  int error(const string &);

  // Spectrum drawing method
  int plot(void);


  void getData(string & data) {data = plotContig.getImage();};


};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////