///////////////////////////////////////////////////////////////////////////////
#ifndef __SPECPLOT_INTERFACE_H__
#define __SPECPLOT_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "CommandLineParser.h"
#include "ParameterList.h"

#include "spectrum.h"
#include "SpecSet.h"
#include "PlotSpectrum.h"
///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;


#define STR(s) #s
#define XSTR(s) STR(s)


///////////////////////////////////////////////////////////////////////////////
class SpecplotInterface {

 protected:

  // Spectrum draw object -- used to draw image
  PlotSpectrum  plotSpectrum;

  // Specset object - used to load specset from file
  SpecSet *specSet;

  // Is the specset object owned by specplot?
  bool specSet_own;

  // Spectrum index.
  int m_spectrumIndex;

  // spectrum scan
  string m_spectrumScan;


 public:


  // Constructors and destructor
  //! \name CONSTRUCTORS
  //@{
  SpecplotInterface();
  //@}

  //! \name DESTRUCTOR
  //@{
  ~SpecplotInterface();
  //@}

  // Option parsing
  int processOptions(int argc, char **argv);

  // Execution based on options
  int processOptions(ParameterList & commandLineParams);

  // load files
  virtual int load(string &fn, string &ext, bool spectrumPresent);

  // Output help
  int help(ostream &);

  // Output version information
  int version(ostream &);

  // Output error messages
  int error(const string &);

  // Spectrum drawing method
  int plot(void);


  void getData(string & data) {data = plotSpectrum.getImage();};

  void setData(int type, void *data);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
