///////////////////////////////////////////////////////////////////////////////
#ifndef __RENDERER_GNU_H__
#define __RENDERER_GNU_H__
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <sstream>


#include "RendererBase.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////


using namespace std;

///////////////////////////////////////////////////////////////////////////////
class RendererGnu : public RendererBase {

  // Where to hold GnuRenderer command stream (command stack)
  stringstream   m_commandStack, m_gnuplotQueue;

  // Stores gnuRenderer command file name
  string          m_gnuplotCommandFileName;

  // Keeps track of terminator
  bool  m_terminatorIssued;

  // Output helpers section
  virtual void outputVector(stringstream &out, vector<string> &);
  virtual void outputDoubleVector(stringstream &out, vector<pair<string,string> > &);

  // Sets GnuRenderer Command File name
  virtual void setGnuplotCommandFileName(void);
  // Generate gnuplot command queue
  virtual void createGnuplotCommandList(void);
  // Writes command file for GnuRenderer
  virtual void writeGnuplotFile(void);
  // invokes GnuRenderer
  virtual void callGnuplot(void);
  virtual void callGnuplot2(void);
  // Removes gnuRenderer command file
  virtual void deleteGnuplotCommandFile(void);

  // Axis range tail
  virtual void setAxisRange(stringstream &out, const double &RangeStart, const double &RangeEnd);
  virtual void setAxisRange(stringstream &out, const double &RangeStart);

  // Method for partial arrow and line drawing methods
  virtual void drawArrowInternal(stringstream &out, const RendererLine &line);

  // output a single coordinate
  virtual void drawCoordinate(stringstream &out, const RendererCoordinate &coordinate);

  // output item order string
  virtual void outputOrder(stringstream &out, const RendererDataOrder &order);


  // output a single point
  virtual void drawPoint(stringstream &out, const RendererPoint &RendererPoint);

  // get file format command for gnuplot
  virtual string getFileFormat(void);



  // sets the title
  virtual void setTitle(stringstream &out, const char *title);
  virtual void setTitle(stringstream &out, const RendererTitle &title);

  // Clear command stack
  virtual void clear(stringstream &out);

  // Draw a single line
  virtual void drawLine(stringstream &out, const RendererLine &line);
  // Draw an arrow
  virtual void drawArrow(stringstream &out, const RendererLine &line);
  // Draw a label
  virtual void drawLabel(stringstream &out);
  // Draw a curve
  virtual void drawCurve(stringstream &out);
  // Draw graph axis
  virtual void drawAxis(stringstream &out);

  // Draw data
  virtual void drawData(stringstream &out, vector<RendererData> &data);

  // Set margins
  virtual void setMarginLeft(stringstream &out, const double &margin, bool src = false);
  virtual void setMarginRight(stringstream &out, const double &margin, bool src = false);
  virtual void setMarginBottom(stringstream &out, const double &margin, bool src = false);
  virtual void setMarginTop(stringstream &out, const double &margin, bool src = false);
  virtual void setMarginAll(stringstream &out, const double &margin, bool src = false);

  // Set origin and size
  virtual void setOrigin(stringstream &out, double &x, double &y);
  virtual void setSize(stringstream &out, double &x, double &y);

  // Set draw area
  virtual void writeDataSection(stringstream &out, RendererArea &area);

  // Set Border
  virtual void setBorder(stringstream &out, const int &border);

  // Set axis range
  virtual void setXrange(stringstream &out, const double &, const double &);
  virtual void setXrange(stringstream &out, const double &);
  virtual void setYrange(stringstream &out, const double &, const double &);
  virtual void setYrange(stringstream &out, const double &);

  // Set axis labels
  virtual void setXlabel(stringstream &out, const string &);
  virtual void setXlabel(stringstream &out, const char *);
  virtual void setYlabel(stringstream &out, const string &);
  virtual void setYlabel(stringstream &out, const char *);

  virtual void setLabel(stringstream &out, const RendererLabel &label);


  // set axis intervals
  virtual void setXinterval(stringstream &out);
  virtual void setYinterval(stringstream &out);

  virtual void setXinterval(stringstream &out, const double &);
  virtual void setYinterval(stringstream &out, const double &);

  virtual void setXinterval(stringstream &out, vector<pair<string,string> > &);
  virtual void setYinterval(stringstream &out, vector<pair<string,string> > &);

  virtual void setXintervalOut(stringstream &out);
  virtual void setYintervalOut(stringstream &out);



 public:


  // Constructor
  RendererGnu();
  // Desctructor
  ~RendererGnu();

  // initialize object
  virtual void initialize(void);

  // Execute commands in the stack
  virtual int execute(void);

  // Clear command stack
  virtual void clear(void)       {m_commandStack.str("");m_gnuplotQueue.str("");RendererBase::clear();};

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
