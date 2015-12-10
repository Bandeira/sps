///////////////////////////////////////////////////////////////////////////////
#ifndef __SPECPLOT_INTERFACE_H__
#define __SPECPLOT_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "spectrum.h"
#include "PlotSpectrum.h"
///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
class SpecplotInterface {


  // Spectrum draw object -- used to draw image
  PlotSpectrum  plotSpectrum;

  // Specset object - used to load specset from file
  SpecSet specSet;

  // Spectrum index.
  int m_spectrumIndex;

  // spectrum scan
  string m_spectrumScan;


 public:


  // Constructors and destructor
  SpecplotInterface();
  ~SpecplotInterface();

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


  void getData(string & data) {data = plotSpectrum.getImage();};


};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
