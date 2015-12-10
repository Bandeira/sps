///////////////////////////////////////////////////////////////////////////////
#include "alignplot.h"

#include "AlignplotInterface.h"
#include "PlotContig.h"
#include "abruijn.h"

///////////////////////////////////////////////////////////////////////////////
using namespace specnets;
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  // Build PlotCont object -- used to draw image
  AlignplotInterface  alignCont;

  // Process options
  alignCont.processOptions(argc, argv);

  return 0;
}
///////////////////////////////////////////////////////////////////////////////
