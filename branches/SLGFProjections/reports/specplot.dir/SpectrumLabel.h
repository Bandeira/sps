///////////////////////////////////////////////////////////////////////////////
#ifndef __LABEL_PLACER_H__
#define __LABEL_PLACER_H__
///////////////////////////////////////////////////////////////////////////////

#include "PlotBase.h"

#include "spectrum.h"
#include "aminoacid.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Data structure used to simplify and speed up parameter passing
struct ParamsData2 {
  vector< vector<int> >&matrix;
  int matX;
  int matY;
  int x1;
  int y1;
  int i;
  int currentX;
  int currentY;
  int labelSizeX;
  int labelSizeY;

  ParamsData2(vector< vector<int> >&m) : matrix(m) {};
};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class SpectrumLabel : public PlotBase  {


  // vector to hold peak information (for peak label generation)
  vector<PeakLabelItem>           peakLabelList;



  void placePeakLabels(void);

  // Checks if a specific position is already taken (peak label placement)
  bool checkPlacementPosition(ParamsData2 &params);
  // Checks all the possible positions for a peak label
  void checkPlacementPositions(ParamsData2 &params);


  // mz lower and upper limits
  double m_mzLowerLimit;
  double m_mzUpperLimit;

  
 protected:

  // Constructors and Destructor
  SpectrumLabel(RendererBase *renderer, double ll, double ul) :
    m_mzLowerLimit(ll), m_mzUpperLimit(ul)
    {m_rendererObject = renderer; m_ownRenderer = false;};

  ~SpectrumLabel() {};


  // Default draw object entry point
  virtual void drawExec(void);


};
///////////////////////////////////////////////////////////////////////////////
}; //namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
