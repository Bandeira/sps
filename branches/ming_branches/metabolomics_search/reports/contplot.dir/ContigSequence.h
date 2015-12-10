///////////////////////////////////////////////////////////////////////////////
#ifndef __CONTIG_SEQUENCE_H__
#define __CONTIG_SEQUENCE_H__
///////////////////////////////////////////////////////////////////////////////
#include "spectrum.h"
#include "abruijn.h"
#include "PlotBase.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
// Auxiliary classes and methods

///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
class ContigSequence : public PlotBase {


  // sequence
  string *m_sequence;

  // mass intervals
  vector<double> *m_breaks;

  // previous sequence intervals. Null if no vertical dashed lines are to be drawn
  vector<double> *m_previousBreaks;

  // the peptide label
  string m_label;

  // y graph position to draw the sequence
  double m_yPosition;

  // y graph position of the previous sequence
  double m_yPositionPrevious;

  // sequence mass offset
  double m_offset;

  // specifies if the label is to be drawn
  bool m_drawLabel;

  // specifies the color used to draw the arrow and the dashed lines
  string m_colour;

  // specifies if the conneting dashed line is bend
  bool m_bendLine;

  // mz lower and upper limits
  double m_mzLowerLimit;
  double m_mzUpperLimit;


  // clean the sequence, bu eliminating non-existing characters
  void cleanSequence(string &iSequence);

  // draws peptide labels
  void drawPeptideLabel(void);

  // Draws the upper peptides (consensus and homolog)
  void drawPeptide(void);

  // draw a red arrow
  void addRedArrow(double x1, double y1, double x2, double y2);

  // draws the top red arrows
  void drawTopArrows(void);

  // break the sequence into a vector of strings
  void breakSequence(vector<string> &oSequenceMapping, string &iSequence, bool omitValues, bool useParantesis = false);

  // draw a single dashed line
  void addDashedSegment(double x1, double y1, double x2, double y2);

  // draw a single dashed line bending to the right
  void addDashedSegmentBent(double x1, double y1, double x2, double y2);

  // draw the dashed lines conneting the current sequence with the previous
  void drawDashedLines(void);
  

 public:


  // Constructors and Destructor
  ContigSequence(RendererBase *renderer, vector<double> *b, vector<double> *pb, string *s, double o, char *l, double y, double yp, bool dl, double ll, double ul, string c, bool bl = false) :
    m_breaks(b), m_previousBreaks(pb), m_sequence(s), m_label(l), m_yPosition(y), m_yPositionPrevious(yp), m_offset(o), m_drawLabel(dl), m_colour(c), m_bendLine(bl), m_mzLowerLimit(ll), m_mzUpperLimit(ul)
    {m_rendererObject = renderer; m_ownRenderer = false;};

  ~ContigSequence() {};


  // Default draw object entry point
  virtual void drawExec(void);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
