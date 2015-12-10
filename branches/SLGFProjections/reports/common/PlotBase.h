///////////////////////////////////////////////////////////////////////////////
#ifndef __PLOT_BASE_H__
#define __PLOT_BASE_H__
///////////////////////////////////////////////////////////////////////////////

#include "RendererBase.h"
#include "Defines.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
#define ANNOTATION_STYLE_SPECNETS   1
#define ANNOTATION_STYLE_INSPECT    2

#define DEFAULT_ANNOTATION_FILE_LOCATION "."




///////////////////////////////////////////////////////////////////////////////
// PlotBase class
//
// Base class for all plot object classes. Defines base methods for drawing
//
class PlotBase {

 protected:

  // debug state
  int m_debug;


  // Base class pointer for renderer object. Should be initialized as a parameter on object contruction.
  RendererBase   *m_rendererObject;

  // is this object owner of the renderer object?
  bool            m_ownRenderer;

  // Image properties
  RendererImage   m_rendererImage;

  // Image contents

  // axis sizes
  double m_axisSizeX;
  double m_axisSizeY;

  // Specifies if the title is printed
  int     m_titlePos;
  // Specifies the title
  string  m_title;

  // Top margin position
  double m_topMarginPosition;

  // Title offsets
  double m_titleOffsetX;
  double m_titleOffsetY;

  // Mass shift value
  double  m_massShift;

  // Peak mass tolerance
  double  m_peakMassTol;

  // peptide annotation style
  int m_annotationStyle;


  //////////////////////////////////////////////////////////////////////////////
  // input file location and name

  // Amino-acids file
  string  m_aminoAcidsFile;
  //  Amino-acids file location
  string  m_aminoAcidsFileDirectory;

  // Annotation model used
  string  m_annotationModel;
  // Annotation model file location
  string  m_annotationModelDirectory;
  // if true, the annotation model is acquired from the file
  bool m_readModelFromData;

  // Renderer directory
  string m_rendererLocation;
  // Font file directory
  string m_fontLocation;




  // Image file attributes

  // File name
  string m_fileName;
  // output directory
  string m_fileOutDir;
  // File prefix for file name formatting
  string m_filePrefix;
  // File suffix for file name formatting
  string m_fileNameSuffix;
  // File extension (defines file format)
  string m_fileFormat;
  // File encoding
  string m_encoding;
  // output target
  string m_target;

  // generated image stored internally
  string m_image;



  virtual void initializePixelSize(void);


  // Draws a generic label
  virtual void drawLabel(string theLabel,
                          double x1, RendererCoordType coordX1,
                          double y1, RendererCoordType coordY1,
                          double offSetX,
                          double offSetY,
                          RendererLabelOffsetType  location,
                          string color,
                          RendererDataOrder order = RENDERER_DATA_ORDER_NONE);

  // Draws a generic line
  virtual void drawLine(int lineWidth, int lineType,
                        double x1, RendererCoordType coordX1,
                        double y1, RendererCoordType coordY1,
                        double x2, RendererCoordType coordX2,
                        double y2, RendererCoordType coordY2,
                        string color,
                        RendererDataOrder order = RENDERER_DATA_ORDER_NONE);

  // draw a generic arrow
  virtual void drawArrow(int lineWidth, int lineType,
                        double x1, RendererCoordType coordX1,
                        double y1, RendererCoordType coordY1,
                        double x2, RendererCoordType coordX2,
                        double y2, RendererCoordType coordY2,
                        string color,
                        RendererDataOrder order = RENDERER_DATA_ORDER_NONE);

 public:


  // Default constructor
  PlotBase(void);
  // Class Destructor.
  ~PlotBase(void);

  // debug method
  virtual void setDebug(int d) {m_debug = d;};

  // Parameter setting methods
  virtual void setFileName(const string &fn)    {m_fileName   = fn;};
  virtual void setFnSuffix(const string &fns)   {m_fileNameSuffix = fns;};

  virtual void setFileOutDir(const string &od)  {m_fileOutDir = od;};
  virtual void setFilePrefix(const string &fp)  {m_filePrefix = fp;};
  virtual void setFileFormat(const string &ff)  {m_fileFormat = ff;};
  virtual void setEncoding(const string &e)     {m_encoding   = e; };
  virtual void setTarget(const string &t)       {m_target     = t; };

  virtual void setZoom(const double &z)         {m_rendererImage.zoom = z;};
  virtual void setImageDimensions(const int x, const int y)     {m_rendererImage.imageSizeX = x;  m_rendererImage.imageSizeY  = y; initializePixelSize();};
  virtual void setImageHeight(const int y)                      {m_rendererImage.imageSizeY = y; initializePixelSize();};
  virtual void setImageWidth(const int x)                       {m_rendererImage.imageSizeX = x; initializePixelSize();};
  virtual void setPixelSize(const double &pw, const double &ph) {m_rendererImage.pixelWidth = pw; m_rendererImage.pixelHeight = ph;};

  virtual void setImageStretchHeight(const int y)               {m_rendererImage.imageStretchLimitY = y;};
  virtual void setImageStretchWidth(const int x)                {m_rendererImage.imageStretchLimitX = x;};


  virtual void setFontName(const string &fn)    {m_rendererImage.fontName    = fn;};
  virtual void setFontHeight(const int fh)      {m_rendererImage.fontHeight  = fh;};
  virtual void setFontWidth(const int fw)       {m_rendererImage.fontWidth   = fw;};
  virtual void setFontSize(const double &fs)    {m_rendererImage.fontSize    = fs;};
  virtual void setFontColor(const int fc)       {m_rendererImage.fontColor   = fc;};

  virtual void setTitlePresence(bool tp)                {m_titlePos = tp;}; // True for no title; false to use the title
  virtual void setTitle(const string &t)                {m_title    = t;};
  virtual void setTopMarginPosition(const double &tmp)  {m_topMarginPosition = tmp;};

  // Mass
  virtual void setPeakMassTol(float m)                {m_peakMassTol              = m;    };
  virtual void setMassShift(float m)                  {m_massShift                = m;    };


  // Peptide
  virtual void setAnnotatinStyle(int s)               {m_annotationStyle          = s;    };


  // File & directory

  virtual void setAminoAcidsFile(string f)            {m_aminoAcidsFile           = f;    };
  virtual void setAminoAcidsFileDirectory(string l)   {m_aminoAcidsFileDirectory  = l;    };
  virtual void setAnnotationModel(string f)           {m_annotationModel          = f;    };
  virtual void setAnnotationModelDirectory(string l)  {m_annotationModelDirectory = l;    };
  virtual void setModelFromFile(void)                 {m_readModelFromData        = false;};
  virtual void setRendererLocation(string l)          {m_rendererLocation         = l;    };
  virtual void setFontLocation(string l)              {m_fontLocation             = l;    };

  // renderer object delegation


  // Parameter access methods
  virtual string &getFileFormat(void)           {return m_fileFormat;};
  virtual double  getZoom(void)                 {return m_rendererImage.zoom;};

  // get a reference to the internal image
  virtual const string &getImage(void) const {return m_image;};

  // Object initializer
  virtual void initialize(void);

  // clear the object
  virtual void clear(void) {if(m_rendererObject) m_rendererObject->clear();};

  // Initial image properties
  virtual void setDrawProperties(void);

  // Define renderer
  virtual void setDrawMethod(RendererType);

  // Default draw method
  virtual void draw(void);
  // Draw method, using a file as target.
  virtual void draw(char *filename);
  // draw execution
  virtual void drawExec(void) = 0;

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
