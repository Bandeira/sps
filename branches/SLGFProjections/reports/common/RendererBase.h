///////////////////////////////////////////////////////////////////////////////
#ifndef __RENDERER_BASE_H__
#define __RENDERER_BASE_H__
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
using namespace std;

typedef enum {RENDERER_TYPE_NONE, RENDERER_TYPE_GNU} RendererType;

typedef enum { RENDERER_DATA_ORDER_NONE, RENDERER_DATA_ORDER_FRONT, RENDERER_DATA_ORDER_BACK } RendererDataOrder;


typedef enum { RENDERER_DATA_TYPE_INT, RENDERER_DATA_TYPE_DOUBLE } RendererDataType;

typedef enum { RENDERER_COORD_NONE, RENDERER_COORD_FIRST, RENDERER_COORD_SECOND, RENDERER_COORD_SCREEN, RENDERER_COORD_GRAPH } RendererCoordType;

typedef enum { RENDERER_OFFSET_RIGHT, RENDERER_OFFSET_CENTER, RENDERER_OFFSET_LEFT } RendererLabelOffsetType;


///////////////////////////////////////////////////////////////////////////////
// Renderer Coordinate
//
// Stores information about a coordinate: value and coordinate system used
//
struct RendererCoordinate {
  double        value;
  RendererCoordType type;
};
///////////////////////////////////////////////////////////////////////////////
// Point class
//
// Stores a 2D point coordinates
//
struct RendererPoint {
  RendererCoordinate   x, y;
};
///////////////////////////////////////////////////////////////////////////////
// Line class
//
// Stores line info
//
struct RendererLine {
  // Start point coordinates
  RendererPoint     start;
  // End point coordinates
  RendererPoint     end;
  // Line type
  int               type;
  // Line width
  int               width;
  // Line color
  string            color;
  // drawing order
  RendererDataOrder order;

  RendererLine() : order(RENDERER_DATA_ORDER_NONE) {};

};
///////////////////////////////////////////////////////////////////////////////
// RendererLabel class
//
// Stores label info
//
struct RendererLabel {
  string                    label;
  RendererPoint             position;
  double                    offsetX;
  double                    offsetY;
  RendererLabelOffsetType   location;
  string                    color;
  // drawing order
  RendererDataOrder         order;

  RendererLabel() : order(RENDERER_DATA_ORDER_NONE) {};

};
///////////////////////////////////////////////////////////////////////////////
// RendererTitle class
//
// Stores title info
//
struct RendererTitle {
  string    label;
  double    offsetX;
  double    offsetY;
  string    color;
  string    fontName;
  double    fontSize;
};

///////////////////////////////////////////////////////////////////////////////
// RendererData class
//
// Plot data container. Contains data to be ploted.
//
struct RendererData {
  string                              title;        // Data title.
  RendererDataType                    dataType;     // Data type: may be int or double
  int                                 dataRecords;  // # of records
  std::vector<pair<double,double> >   dataDouble;   // Data records for double data type
  std::vector<pair<int,int> >         dataInt;      // Data records for int data type
  int                                 lineWidth;    // width of data lines
  int                                 lineType;     // type of data lines
  string                              lineColor;    // color of data lines
};
///////////////////////////////////////////////////////////////////////////////
// Renderer Area - specifies a drawing area
///////////////////////////////////////////////////////////////////////////////
struct RendererArea {

  /////////////////////////////////////////////////////////////////////////////
  // grapth specific margins and location

  // borders
  int    m_border;

  // image margins
  double m_marginLeft, m_marginRight, m_marginBottom, m_marginTop;

  // area size
  bool m_useSize;
  double m_sizeX, m_sizeY;

  // area location
  bool m_useOrigin;
  double m_originX, m_originY;


  //////////////////////////////////////////////////////////////////////////////
  // axis data

  // axis labels
  string m_xLabel, m_yLabel;

  // x and y start and end values
  bool   m_useXstart, m_useXend, m_useYstart, m_useYend;
  double m_startX, m_endX;
  double m_startY, m_endY;

  // data for x tics
  bool   m_useTicsX;
  double m_ticsX;
  vector<pair<string,string> > m_ticsExtendedX;
  // data for y tics
  bool   m_useTicsY;
  double m_ticsY;
  vector<pair<string,string> > m_ticsExtendedY;


  //////////////////////////////////////////////////////////////////////////////
  // data on area

  // general labels
  vector<RendererLabel> m_labels;

  // arbitrary lines
  vector<RendererLine>  m_lines;

  // arbitrary arrow
  vector<RendererLine> m_arrows;

  // data to renderer
  vector<vector<RendererData> > m_data;


  //////////////////////////////////////////////////////////////////////////////
  // Methods

  // constructor
  RendererArea() {reset();};

  void reset(void)
  {
    m_useXstart   = m_useXend     = false;
    m_useYstart   = m_useYend     = false;
    m_border      = 0;
    m_marginLeft  = m_marginRight = m_marginBottom = m_marginTop = 0.0;
    m_useSize     = false;
    m_useOrigin   = false;
    m_xLabel      = ""; m_yLabel  = "";
    m_useTicsX    = m_useTicsY    = false;
    m_ticsExtendedX.clear(); m_ticsExtendedY.clear();
    m_labels.clear();
    m_lines.clear();
    m_arrows.clear();
    m_data.clear();
  };

};

///////////////////////////////////////////////////////////////////////////////
// RendererImage class
//
// Image data container
//
struct RendererImage {

  // Image dimensions
  int     imageSizeX;
  int     imageSizeY;

  //Image stretch limits. Used to limit automatic image stretching
  int     imageStretchLimitX;
  int     imageStretchLimitY;

  // Zoom value for image
  double  zoom;

  // Pixel size. Typically 1.0
  double  pixelWidth;
  double  pixelHeight;

  // Font name, size and default color
  string  fontName;
  int     fontHeight;
  int     fontWidth;
  double  fontSize;
  int     fontColor;

  // File name used to write the .png, .jpg, ...
  string  outputFileName;

  // File type: png, eps, jpg, ...
  string  fileType;


  RendererImage() : zoom(1.0), fileType("png") {};

};
///////////////////////////////////////////////////////////////////////////////
// PlotBase class
//
// Low level drawing wrapper base class.
//
class RendererBase {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // renderer data


  // Defines the renderer location. Defaults to '.'
  string    m_rendererLocation;
  // Defines the font files location. Defaults to '.'
  string    m_fontLocation;
  // defines if LD_LIBRARY_PATH is set to the same directory as to the same directory the executable is located in
  bool      m_setLibraryPath;

  //////////////////////////////////////////////////////////////////////////////
  // extra image properties

  // title data
  bool m_useTitle;
  RendererTitle m_title;


  //////////////////////////////////////////////////////////////////////////////
  // Areas data

  // areas in store to draw
  vector<RendererArea>  m_areas;

  // current area being rendered
  RendererArea          m_areaCurrent;

  // stores data for x broken axis
  //bool    m_breakYAxis;
  //double  m_yAxisStop, m_yAxisRestart;


  //////////////////////////////////////////////////////////////////////////////
  // auxiliary methods

  int m_debug;

  virtual void outputVector(vector<string> &) {};


 public:

  // Image properties
  RendererImage   m_rendererImage;


  //////////////////////////////////////////////////////////////////////////////
  // constructor and destructor

  // Constructors and destructor
  RendererBase(void);
  ~RendererBase(void);

  // Initialization method
  virtual void initialize(void);

  virtual void clear(void) {m_areas.clear();m_areaCurrent.reset();};

  virtual void setDebug(int d) {m_debug = d;};

  //
  virtual int execute(void)                                 {};

  // Internal method to calculate auxiliary values used in label position calculation
  virtual double calculateInterval(double distance, double factor, double zoom);

  // add a new drawing area
  virtual void addArea(void)                                {m_areas.push_back(m_areaCurrent); m_areaCurrent.reset();};

  // set the renderer executable location
  virtual void setRendererLocation(string &l)               {m_rendererLocation = l;};
  // set the font files location
  virtual void setFontLocation(string &f)                   {m_fontLocation = f;};
  // use the LD_LIBRARY_PATH
  virtual void setLibraryPath(void)                         {m_setLibraryPath = true;};
  // do not use the LD_LIBRARY_PATH
  virtual void clrLibraryPath(void)                         {m_setLibraryPath = false;};


  //sets the title
  virtual void setTitle(const RendererTitle &t)             {m_title = t;m_useTitle = true;};

  // Draw a single line
  virtual void drawLine(const RendererLine &line)           {m_areaCurrent.m_lines.push_back(line);};
  // Draw a single arrow
  virtual void drawArrow(const RendererLine &arrows)        {m_areaCurrent.m_arrows.push_back(arrows);};
  // Draw a curve
  virtual void drawCurve(void)                              {};
  // Draw graph axis
  virtual void drawAxis(void)                               {};
  // Draw a label
  virtual void drawLabel(void)                              {};

  // Draw data
  virtual void drawData(vector<RendererData> &data)         {m_areaCurrent.m_data.push_back(data);};

  // Define margins
  virtual void setMarginLeft(const double &margin)          {m_areaCurrent.m_marginLeft    = margin;};
  virtual void setMarginRight(const double &margin)         {m_areaCurrent.m_marginRight   = margin;};
  virtual void setMarginBottom(const double &margin)        {m_areaCurrent.m_marginBottom  = margin;};
  virtual void setMarginTop(const double &margin)           {m_areaCurrent.m_marginTop     = margin;};
  virtual void setMarginAll(const double &margin)           {m_areaCurrent.m_marginLeft    = m_areaCurrent.m_marginRight =
                                                             m_areaCurrent.m_marginBottom  = m_areaCurrent.m_marginTop   = margin;};
  // define origin and size
  virtual void setSize(double ox, double oy)                {m_areaCurrent.m_sizeX   = ox; m_areaCurrent.m_sizeY   = oy; m_areaCurrent.m_useSize   = true;};
  virtual void setOrigin(double sx, double sy)              {m_areaCurrent.m_originX = sx; m_areaCurrent.m_originY = sy; m_areaCurrent.m_useOrigin = true;};

  // Set Border
  virtual void setBorder(const int &border)                   {m_areaCurrent.m_border = border;};

  // define axis range - values covered by axis
  virtual void setXrange(const double &sx, const double &ex)  {m_areaCurrent.m_startX = sx; m_areaCurrent.m_endX = ex;};
  virtual void setXrange(const double &sx)                    {m_areaCurrent.m_startX = sx;};
  virtual void setYrange(const double &sy, const double &ey)  {m_areaCurrent.m_startY = sy; m_areaCurrent.m_endY = ey;};
  virtual void setYrange(const double &sy)                    {m_areaCurrent.m_startY = sy;};

  // define a broken axis
  virtual void breakAxisY(const double until, const double restart);

  // Set axis labels
  virtual void setXlabel(const string &l)                   {m_areaCurrent.m_xLabel = l;};
  virtual void setXlabel(const char *l)                     {m_areaCurrent.m_xLabel = l;};
  virtual void setYlabel(const string &l)                   {m_areaCurrent.m_yLabel = l;};
  virtual void setYlabel(const char *l)                     {m_areaCurrent.m_yLabel = l;};

  // Set an arbitrary label
  virtual void addLabel(const RendererLabel &l)             {m_areaCurrent.m_labels.push_back(l);};

  // set axis tics
  virtual void setXinterval(const double &v)                  {m_areaCurrent.m_ticsX = v; m_areaCurrent.m_useTicsX = true;};
  virtual void setYinterval(const double &v)                  {m_areaCurrent.m_ticsY = v; m_areaCurrent.m_useTicsY = true;};

  virtual void setXinterval(vector<pair<string,string> > &v)  {m_areaCurrent.m_ticsExtendedX = v; m_areaCurrent.m_useTicsX = true;};
  virtual void setYinterval(vector<pair<string,string> > &v)  {m_areaCurrent.m_ticsExtendedY = v; m_areaCurrent.m_useTicsY = true;};

  virtual void disableXinterval(void)                         {m_areaCurrent.m_useTicsX = false;};
  virtual void disableYinterval(void)                         {m_areaCurrent.m_useTicsY = false;};

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
