///////////////////////////////////////////////////////////////////////////////
#ifndef __PLOT_SPECTRUM_H__
#define __PLOT_SPECTRUM_H__
///////////////////////////////////////////////////////////////////////////////

#include "PlotBase.h"

#include "spectrum.h"
#include "aminoacid.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  AUxiliary data structures used
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Internal spectrum representation with several annotations
struct SpectrumAnnotationData {
  int     peptideIndex;
  int     charge;
  int     seriesIndex;
  float   prob;
  string  annotation;
};
////////////////////////////////////////////////////////////////////////////////
struct SpectrumItem {
  double                          peakMass;
  double                          peakIntensity;
  vector<SpectrumAnnotationData>  annotations;
};

// unsigned is the peak index, considering mass incresing ordering
typedef map<unsigned, SpectrumItem> SpectrumData;

////////////////////////////////////////////////////////////////////////////////
// Label data container for label drawing and placement over the graph
struct PeakLabelItem {
  // peak mass value
  double mass;
  // peak intensity value
  double intesity;
  // Annotation probability
  float  prob;
  // label coordinates after processing
  double x,y;
  // Is label placed?
  bool placed;
  // Does the label need a connection line?
  bool drawConnectionLine;
  // label to print
  string annotation;
  // superscript label part
  string superscript;
  // subscript label part
  string subscript;
  // color
  string color;
  // line color
  string lineColor;

  // used to order items - higher to lower intensity, in this case
  bool operator<(const PeakLabelItem &o) const {return intesity < o.intesity;};

};
////////////////////////////////////////////////////////////////////////////////
// Auxiliary data structure used to simplify data drawing and searching
struct DrawingData {
  // ion to search for
  string  ion;
  // specifies if ion is used as a prefix or an exact string to search for
  bool    isPrefix;
  // string to draw
  string  show;
  // color to draw
  int     colorIdx;
  // Color to draw conneting line
  int     lineColorIdx;
  // use B series color schema if true (use y if false)
  bool    useBColors;

  DrawingData(const char *i, bool p, const char *s, int c, int lc, bool bc) :
    ion(i), isPrefix(p), show(s), colorIdx(c), lineColorIdx(lc), useBColors(bc)
    {};

};
////////////////////////////////////////////////////////////////////////////////
// Annotation data container per mass item
// Used to speed up annotations searchs
struct AnnotationData {
  double  peakMass;
  double  peakIntensity;
  int     charge;
  int     seriesIndex;
  float   prob;
  string  annotation;
};
// Mass item data container. Used as an alternative data structure to speed up searches
struct PeptideAnnotation {
  string                  peptideItem;
  double                  mass;
  double                  previousMass;
  vector<AnnotationData>  annotations;
};
////////////////////////////////////////////////////////////////////////////////
// Data structure used to simplify and speed up parameter passing
struct ParamsData1 {
  vector<string>            *ions;
  vector<PeptideAnnotation> *peptideAnnotation;
  string                    *peptide;

  double                    shift;
  double                    yPosition;
  double                    peptideMass;
  string                    color;

  ParamsData1(vector<string> *i, vector<PeptideAnnotation> *a, string *p) : ions(i), peptideAnnotation(a), peptide(p) {};
};
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
// spectrumOutputElem contains information about a spectrum item to be printed
struct spectrumOutputElem {
  string color;
  float  mass;
  float  intensity;
  bool   annotated;

  bool operator<(const spectrumOutputElem &o) const {return intensity > o.intensity;};
};
////////////////////////////////////////////////////////////////////////////////
// ordering data
struct orderData {
  int index;
  double data;

  bool operator<(const orderData &o) const {return data > o.data;};
};
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class PlotSpectrum : public PlotBase {


 protected:


  //////////////////////////////////////////////////////////////////////////////
  // Data properties

  // Intensity upper limit
  double m_intensityLimit;

  // M/Z limits
  double m_mzUpperLimit;
  double m_mzLowerLimit;



  //////////////////////////////////////////////////////////////////////////////
  // Auxiliary data structures used for processing

  // Header positions (2 lines per record)
  // first is b peptide Y axis position
  // second is y peptide Y axis position
  map<int, pair<double,double> >  m_headerPosition;

  //prefix masses, i.e. b ions
  map<int, vector<float> >        m_prmMasses;
  //suffix masses, i.e. y ions
  map<int, vector<float> >        m_srmMasses;
  // peptide global masses
  map<int, float>                 m_peptideMass;

  // Holds the spectrum data with all the annotations per peak
  SpectrumData                    m_spectrumData;

  // processed annotations for module usage, to speed up searches
  map<int, vector<PeptideAnnotation> >  peptideAnnotations;
  map<int, vector<PeptideAnnotation> >  peptideAnnotationsReverse;

  // vector to hold peak information (for peak label generation)
  vector<PeakLabelItem>           peakLabelList;

  // Vector to hold peak information to be output
  vector<spectrumOutputElem>      peaksOutput;


  //////////////////////////////////////////////////////////////////////////////
  // Annotation model

  MS2ScoringModel m_model;


  //////////////////////////////////////////////////////////////////////////////
  // Annotation definition containers

  // Mass dots
  vector<string> ionsSearchDotsB, ionsSearchDotsY;
  // Mass segments
  vector<string> ionsSearchSegmentsB, ionsSearchSegmentsY;
  // Vertical lines
  vector<string> ionsSearchVerticalB, ionsSearchVerticalY;



  //////////////////////////////////////////////////////////////////////////////
  // Spectrum data

  // Data spectrum for output
  Spectrum *m_spectrum;

  // Spectrum to draw
  int     m_spectrumIndex;

  // Spectrum scan number
  string  m_spectrumScan;

  //
  string  m_spectrumInfo;


  //////////////////////////////////////////////////////////////////////////////
  // Peptides data

  // Peptide sequences
  vector<string>  m_peptide;

  //////////////////////////////////////////////////////////////////////////////
  // input file location and name

  // Annotation output data
  stringstream oAnnotFileContent;

  // Annotations output file name
  string m_annotationsOutputFile;
  // Annotations input file name
  string m_annotationsInputFile;

  // annotation order
  bool m_annotationByIntensity;

  //////////////////////////////////////////////////////////////////////////////
  // Graph & axis variables

  // Specified m/z axis range
  double  m_range_min;
  double  m_range_max;

  // Graph size, in pixels. Used in label placement
  int graphSizePixelsX, graphSizePixelsY;

  // main ion representations for N-Term
  string m_mainIon, m_mainIonN;

  // mass offsets
  double m_offset, m_offsetN;



  //////////////////////////////////////////////////////////////////////////////
  // The following methods contain low level calls to the renderer object

  // Used to draw the "y" and "b" labels
  void drawUpperGuideLabel(string label, double coord, string color);
  // Draws a upper line (dashed)
  void drawUpperGuideLine(double x1, double x2, double y1, string color);
  // Draw an upper amino-acids (or peptide items)
  void drawUpperGuideLabelSegments(string label, double x, double y, string color);
  // Draw a upper mass dots
  void drawUpperGuideDots(double x, double y, int size, string color);
  // Draw a upper mass dots for EPS format
  void drawUpperGuideDotsEps(double x, double y, int size, string color);
  // Draw an upper mass line segment
  void drawUpperGuideSegments(double x1, double x2, double y1, string color);

  // Draws the vertical line from the dot to the bottom of the graph (depracated)
  void drawUpperGuideVertical(double x1, double y1, string color);
  // Draws the vertical line from the dot to the top of the peak
  void drawUpperGuideVertical2(double x1, double y1, double x2, double y2, string color);


  //////////////////////////////////////////////////////////////////////////////
  // Methos used to build auxiliary (annotation) data structures

  // Main annotations draw procedure
  void drawAnnotationPeptides(void);

  // initialize spectrumData data structure from Spectrum object
  void initSpectrumData(void);
  // Add annotation data to spectrumData structure
  void addAnnotationData(int index);

  // Annotate spectra using annotation module
  void annotate(string &peptide, psmPtr);

  // Builds annotation vector used for drawing
  void buildAnnotationVector(int peptideIndex, string &peptide, float peptideMass, vector<float> &masses, vector<PeptideAnnotation> &peptideAnnotation, bool reverse);

  // Sets annotations vectors
  void defineAnnotations(void);
  // Sets annotations vector for peak labels
  void definePeakAnnotations(vector<DrawingData> &ionLabels);

    // Parses the peptide string, and returns the next item
  string getMassLabel(string &peptide, int &position, bool reverse);
  // Checks if a peak is annotated
  int  peakAnnotated(vector<DrawingData> &series, string &annotation);
  // Checks if a peak is annotated (2)
  int  peakAnnotated(vector<DrawingData> &series, int peakIndex);
  // Checks if a peak is annotated (3)
  int  findAnnotatedIndex(PeptideAnnotation &annotation, vector<string> &series, bool chargeOne = false);


  //////////////////////////////////////////////////////////////////////////////
  // Methods used to draw graph elements

  // Draws upper peptide
  void drawUpperLabelsAll(ParamsData1 &);
  // Draws mass dots
  void drawUpperMassDots(ParamsData1 &);
  // Draws upper mass segments
  void drawUpperSegments(ParamsData1 &);
  // Draws vertical lines on annotated peaks
  void drawVerticalLines(ParamsData1 &);
  // Draws annotations
  void drawAnnotations(int index);


  //////////////////////////////////////////////////////////////////////////////
  // Methods used to process and draw peak labels

  // Draw peak labels entry point
  void drawPeakLabelsMain(void);

  // process label for lower zoom levels -- removing of mass values and keeping AAs
  void processUpperLabel(string &label);
  // process label for upper zoom levels -- rouding of mass values to 1 decimal places
  void processUpperLabel2(string &label);


  // aquire main peak labels
  void defineMainAnnotations(void);
  // Build peak list to print labels
  void aquirePeakLabels(vector<DrawingData> &ions);
  // Place peak labels so that they don't overlap
  void placePeakLabels(void);
  // Draw the peaks labels
  void drawPeakLabels(void);

  // Checks if a specific position is already taken (peak label placement)
  bool checkPlacementPosition(ParamsData2 &params);
  // Checks all the possible positions for a peak label
  void checkPlacementPositions(ParamsData2 &params);


  //////////////////////////////////////////////////////////////////////////////
  // Methods used to process and draw peaks

  // Draw peaks entry point
  void drawSpectrumPeaks(void);

  // Builds a structure based on peakColors. This structure is used to produce output
  void preparePeaksToOutput(void);

  // Draws the spectrum lines sending "drawLine" commands to the draw engine. Slow.
  void drawSpectrumLines(void);
  // Draws the spectum lines using a binary array data method. Fast.
  void drawSpectrumLinesFast(void);


  //////////////////////////////////////////////////////////////////////////////
  // Auxiliary methods (misc)

  // Gets a color set used for rendering, given a peptide index
  void getColors(int index, char **cb[], char **cy[]);
  // Get a single color at a specified protein index for a specified position. Flag means b or y reference
  char *getColor(int index, int position, bool b = true);

  // Tests if a given mass point is within parameters
  bool testMass(double &value);
  // Tests if a given mass set is within parameters
  bool testMass(double &left, double &right);
  // Tests if a given mass point is within parameters, when peptide mass is not known.
  bool testMass2(double &value);


  //////////////////////////////////////////////////////////////////////////////
  // Annotation file output methods

  // Add annotation data to output stream
  void buildAnnotationOutput(void);
  // Write the annotation data file
  void writeAnnotationFile(void);
  // Read the annotation data file
  void readAnnotationFile(void);


  //////////////////////////////////////////////////////////////////////////////
  // Initialization methods

  // Calculates graph bounds
  void calcultateLimits(void);
  // Initialize several graph parameters
  void initializeGraph(void);
  // Sets the graph title
  void setGraphTitle(void);
  // Calculate b and y Y coordinate positions, and graph height
  void calcDrawPosition(void);


  //////////////////////////////////////////////////////////////////////////////
  // Debug

  void dumpPeptideAndMasses(void);


  //////////////////////////////////////////////////////////////////////////////
  // Main draw routine
  void plotSpectrum(void);


 public:


  //////////////////////////////////////////////////////////////////////////////
  // Constructors and destructor
  PlotSpectrum(void);
  ~PlotSpectrum(void);

  //////////////////////////////////////////////////////////////////////////////
  // Variable initializers

  // clear object
  virtual void clear(void)
  {
    PlotBase::clear();

    m_headerPosition.clear();
    m_prmMasses.clear();
    m_srmMasses.clear();
    m_peptideMass.clear();
    peptideAnnotations.clear();
    peptideAnnotationsReverse.clear();
    peakLabelList.clear();
    peaksOutput.clear();
    ionsSearchDotsB.clear();
    ionsSearchDotsY.clear();
    ionsSearchSegmentsB.clear();
    ionsSearchSegmentsY.clear();
    ionsSearchVerticalB.clear();
    ionsSearchVerticalY.clear();
    m_peptide.clear();
  };

  //////////////////////////
  // Spectrum
  virtual void setSpectrum(Spectrum *s)               {m_spectrum                 = s;    };
  virtual void setSpectrumIndex(int idx)              {m_spectrumIndex            = idx;  };
  virtual void setSpectrumScan(string &ss)            {m_spectrumScan             = ss;   };
  virtual void setSpectrumInfo(string ss)             {m_spectrumInfo             = ss;   };

  // Peptide
  virtual void addPeptide(string &p)                  {m_peptide.push_back(p);            };

  // File & directory
  virtual void setAnnotationOutputFile(string a)      {m_annotationsOutputFile    = a;    };
  virtual void setAnnotationInputFile(string a)       {m_annotationsInputFile     = a;    };

  // misc
  virtual void setAnnotationByIntensity(void)         {m_annotationByIntensity = true;};

  // Graph
  virtual void setRangeMax(float r)                   {m_range_max                = r;    };
  virtual void setRangeMin(float r)                   {m_range_min                = r;    };


  //////////////////////////////////////////////////////////////////////////////
  // Default draw object
  virtual void drawExec(void);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
