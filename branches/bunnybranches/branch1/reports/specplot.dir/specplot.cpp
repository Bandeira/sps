///////////////////////////////////////////////////////////////////////////////
#include "specplot.h"

#include "SpecplotInterface.h"
#include "PlotSpectrum.h"

#include "Timer.h"
///////////////////////////////////////////////////////////////////////////////
using namespace specnets;
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //Timer_c t0;
  // Build PlotSpec object
  SpecplotInterface  specplotInterface;
  // Process options
  int ret = specplotInterface.processOptions(argc, argv);
  // test return value
  if(ret != 0)
    return 0;
  // generate the image
  //cout << "Time to process parameters: " << t0.stop() << endl;
  //Timer_c t1;

  specplotInterface.plot();
  //cout << "Time to generate the image: " << t1.stop() << endl;
  // exit
  return 0;

/*

  // Spectrum draw object -- used to draw image
  PlotSpectrum  plotSpectrum;

  // Specset object - used to load specset from file
  SpecSet specSet;

  // Load specta data
  specSet.LoadSpecSet_pklbin("specs_ms.pklbin");

  // Spectrum index to draw
  plotSpectrum.m_spectrumIndex = 3197;

  // Sets draw object in DrawSpectrum as PlotGnu
  plotSpectrum.setDrawMethod(RENDERER_TYPE_GNU);

  plotSpectrum.m_peptide.push_back("[33.9](IV,-18.0)MSQSPSSLAVSAGEK");
  plotSpectrum.m_peptide.push_back("GGGGGGGGGGGIH");
  plotSpectrum.m_peptide.push_back("AMDQEPOGFMH");
  plotSpectrum.m_peptide.push_back("LNGFGIHEAERG");
  plotSpectrum.m_peptide.push_back("AGHAETATHAFADG");


  plotSpectrum.setZoom(1);

  plotSpectrum.m_spectrum = &specSet.specs[3197];

  // Draws the image
  plotSpectrum.draw("Test.png"); */
}
///////////////////////////////////////////////////////////////////////////////
