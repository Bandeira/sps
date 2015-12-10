////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_PDF_H__
#define __REPORT_RENDERER_PDF_H__
////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererBase.h"
#include "hpdf.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// ReportRendererPdf Helper classes
////////////////////////////////////////////////////////////////////////////////
#define PDF_DEFAULT_FONT_NAME "Helvetica"
#define PDF_PAGE_MARGIN 20
#define PDF_TABLE_CELL_SPACING 1
#define PDF_TABLE_CELL_PADDING 2


struct PdfPoint {

  float m_x ,m_y;


  PdfPoint() {m_x = m_y = 0.0;};
  //PdfPoint(const PdfPoint &p)  : m_x(p.m_x), m_y(p.m_y) {};
  //PdfPoint(float x, float y) : m_x(x),  m_y(y)  {};
  //~PdfPoint() {};

  PdfPoint  operator+(const PdfPoint &o);
  //PdfPoint &operator=(const PdfPoint &o);

};
////////////////////////////////////////////////////////////////////////////////
struct PdfLine {

  PdfPoint m_start, m_end;
  float m_width;
  float m_color_red, m_color_gree, m_color_blue;

  PdfLine() : m_width(1.0), m_color_red(0.0), m_color_gree(0.0), m_color_blue(0.0) {};

};
////////////////////////////////////////////////////////////////////////////////
struct PdfRect {
  float left, right, top, bottom;
};
////////////////////////////////////////////////////////////////////////////////
struct PdfRGB {
  float m_red, m_green, m_blue;

  PdfRGB() : m_red(0.0), m_green(0.0), m_blue(0.0) {};
};
////////////////////////////////////////////////////////////////////////////////
struct PdfTextData {
  PdfRect rect;
  PdfRect area;
  string  text;
  PdfRGB  color;
  string  font;
  int     fontSize;
  int     align;
};
////////////////////////////////////////////////////////////////////////////////
struct PdfImageData {

  PdfRect     rect;
  HPDF_Image  image;

};
////////////////////////////////////////////////////////////////////////////////
class PdfItem {

 protected:

  PdfPoint  m_size;
  PdfPoint  m_location;

 public:

  PdfItem() {};
  PdfItem(PdfPoint &l) : m_location(l) {};
  PdfItem(PdfPoint &l, PdfPoint &s) : m_location(l), m_size(s) {};
  ~PdfItem() {};

  virtual PdfPoint getSize(void) {return m_size + m_location;};

  virtual void calcSize(HPDF_Doc &pdf) {};

};
////////////////////////////////////////////////////////////////////////////////
class PdfText : public PdfItem {

  string    m_text;
  string    m_font;
  int       m_fontSize;
  int       m_font_descent;
  int       m_font_ascent;

 public:

  PdfText() {init();};
  PdfText(const PdfText &t) {init(); m_text = t.m_text;};
  PdfText(string &t) : m_text(t) {init();};
  PdfText(string &t, PdfPoint &l) : m_text(t), PdfItem(l) {init();};
  PdfText(string &t, PdfPoint &l, PdfPoint &s) : m_text(t), PdfItem(l, s) {init();};

  void init(void) {m_font = PDF_DEFAULT_FONT_NAME; m_fontSize = 10;};

  void setText(string &t)   {m_text = t;};
  void setFont(string &f)   {m_font = f;};
  void setFontSize(int fs)  {m_fontSize = fs;};

  int getAscent(void)       {return m_font_ascent;};
  int getDescent(void)      {return m_font_descent;};
  string &getText(void)     {return m_text;};
  string &getFont(void)     {return m_font;};
  int getFontSize(void)     {return m_fontSize;};

  virtual void calcSize(HPDF_Doc &pdf);

};
////////////////////////////////////////////////////////////////////////////////
class PdfImage : public PdfItem {

  HPDF_Image  m_image;
  float       m_factor;

 public:

  PdfImage() : m_factor(1.0) {};
  PdfImage(const PdfImage &i) {m_image = i.m_image; m_factor = i.m_factor;};
  PdfImage(HPDF_Image &i) : m_image(i), m_factor(1.0) {};
  PdfImage(HPDF_Image &i, PdfPoint &l) : m_image(i), PdfItem(l), m_factor(1.0) {};
  PdfImage(HPDF_Image &i, PdfPoint &l, PdfPoint &s) : m_image(i), PdfItem(l, s), m_factor(1.0) {};

  virtual void calcSize(HPDF_Doc &pdf);
  virtual HPDF_Image getImage(void) {return m_image;};

  virtual void setFactor(float f) {m_factor = f;};
  virtual float getFactor(void) {return m_factor;};
};
////////////////////////////////////////////////////////////////////////////////
class PdfTableElem {

 protected:

  PdfPoint          m_size;

 public:

  PdfTableElem() {};
  PdfTableElem(PdfPoint &s) : m_size(s) {};

  virtual void calcSize(HPDF_Doc &pdf) {};
  virtual PdfPoint &getSize(void) {return m_size;};

};
////////////////////////////////////////////////////////////////////////////////
class PdfTableCell : public PdfTableElem {

  vector<PdfText>   m_text;
  vector<PdfImage>  m_image;

 public:

  PdfTableCell() {};
  PdfTableCell(const PdfTableCell &c) {m_text = c.m_text ; m_image = c.m_image;};

  void addText(PdfText &t)   {m_text.push_back(t); };
  void addImage(PdfImage &i) {m_image.push_back(i);};

  int getTextLength(void) {return m_text.size();};
  int getImageLength(void) {return m_image.size();};

  PdfText  &getTextAt(int i) {return m_text[i];};
  PdfImage &getImageAt(int i) {return m_image[i];};

  virtual void calcSize(HPDF_Doc &pdf, float cellPadding);

};
////////////////////////////////////////////////////////////////////////////////
class PdfTableRow : public PdfTableElem {

  vector<PdfTableCell>  m_cell;

 public:

  PdfTableRow() {};
  PdfTableRow(const PdfTableRow &r) {m_cell = r.m_cell;};

  PdfTableCell &operator[](int i) {return m_cell[i];};

  void addCell(PdfTableCell &c) {m_cell.push_back(c);};

  virtual int length(void) {return m_cell.size();};
  virtual void calcSize(HPDF_Doc &pdf, float cellSpacing, float cellPadding);
  virtual PdfPoint &getSizeOf(int idx) {return m_cell[idx].getSize();};

};
////////////////////////////////////////////////////////////////////////////////
class PdfTable : public PdfTableElem {

  vector<PdfTableRow> m_row;

  vector<float>  m_Rows;
  vector<float>  m_Cols;

  float m_tableWidth;
  float m_tableHeight;
  float m_cellSpacing;
  float m_cellPadding;

 public:

  PdfTable() : m_cellSpacing(PDF_TABLE_CELL_SPACING), m_cellPadding(PDF_TABLE_CELL_PADDING) {};

  PdfTableRow &operator[](int i) {return m_row[i];};

  void addRow(PdfTableRow &r) {m_row.push_back(r);};

  virtual void calcSizes(HPDF_Doc &pdf);

  virtual vector<float> &getCols(void) {return m_Cols;};
  virtual float getRowAt(int idx) {return m_Rows[idx];};

  virtual float getTableWidth(void) {return m_tableWidth;};
  virtual float getCellSpacing(void) {return m_cellSpacing;};
  virtual float getCellPadding(void) {return m_cellPadding;};

  virtual int length(void) {return m_row.size();};

};
////////////////////////////////////////////////////////////////////////////////
// ReportRendererPdf
////////////////////////////////////////////////////////////////////////////////
class ReportRendererPdf : public ReportRendererBase {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // Data structures

  // pointer to the PDF document
  HPDF_Doc *m_pdf;

  // the tables preprocessed for PDF gwnwration.
  vector<PdfTable> pdfTable;

  // page header image
  // Logo image
  HPDF_Image imageHeaderRight;
  HPDF_Image imageHeaderLeft;
  HPDF_Image imageHeaderCenter;

  // page counter
  int m_pageIndex;

  // PDF current position (height)
  float positionH;

  // PDF page bottom limit - no writting below this line
  float PdfLimitBottom;

  // page pargin
  float m_pageMargin;


  //////////////////////////////////////////////////////////////////////////////
  // PDF interface added 0:methods
  int drawLine(HPDF_Page &page, PdfLine &line);

  int drawText(HPDF_Doc &pdf, HPDF_Page &page, PdfTextData &data);

  int drawTextRect(HPDF_Doc &pdf, HPDF_Page &page, PdfTextData &data);

  int drawImage(HPDF_Doc &pdf, HPDF_Page &page, PdfImageData &data);



  //////////////////////////////////////////////////////////////////////////////
  // Report high level mehtods

  virtual int createPdfPage(HPDF_Doc &pdf);

  // prolog
  virtual int renderReportProlog(HPDF_Doc &pdf);


  //////////////////////////////////////////////////////////////////////////////
  // Report page prolog and epilogue

  // prolog
  virtual int renderPageProlog(HPDF_Doc &pdf, HPDF_Page &page, HPDF_Image image, PdfRect &r);
  // epilog
  virtual int renderPageEpilog(HPDF_Doc &pdf, HPDF_Page &page);

  //////////////////////////////////////////////////////////////////////////////
  // Report table prolog and epilogue

  // prolog
  virtual int renderTableProlog(ReportBase *table, HPDF_Doc &pdf);
  // epilog
  virtual int renderTableEpilog(ReportBase *table, HPDF_Doc &pdf);
  // between tables
  virtual int renderInterTable(ReportBase *table, HPDF_Doc &pdf);


  //////////////////////////////////////////////////////////////////////////////
  // Table content renderers

  // render table common method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
  virtual int renderTable(ReportTableBase *table, PdfTable &pdfTable, HPDF_Doc &pdf, float pageMargin);
  // render a header row method. Renders a table header row.
  virtual int renderTableHeaderRow(PdfTableRow &row, HPDF_Doc &pdf, float rowHeight, vector<float> &sizes, float pageMargin, float tableWidth, float cellSpacing, float cellPadding);
  // render header cell method. Renders a header cell
  virtual int renderTableHeaderCell(PdfTableCell &cell, HPDF_Doc &pdf, HPDF_Page &page, float, float cellPadding, float height, float width);
  // build a row common method. Cycles thu all columnType cells and buildCell for each ColumnType item
  virtual int renderTableRow(PdfTableRow &row, HPDF_Doc &pdf, float rowHeight, vector<float> &sizes, float pageMargin, float tableWidth, float cellSpacing, float cellPadding);
  // build cell common method. builds a specific cell based on ColumnType specifications and a row of data
  virtual int renderTableCell(PdfTableCell &cell, HPDF_Doc &pdf, HPDF_Page &page, float, float cellPadding);


  //////////////////////////////////////////////////////////////////////////////
  // Table header builders

  // render a header row common method. Cycles thu all columnType cells and renders each one.
  virtual int buildTableHeaderRow(ReportTableBase *table, PdfTableRow &pdfTableRow);
  // render header cell common method. Renders all ColumnTypes
  virtual int buildTableHeaderCell(ReportColumnTypeBase *base, PdfTableCell &pdfTableCell);


  //////////////////////////////////////////////////////////////////////////////
  // Table content builders

  // build table common method. Cycles thru all the table rows and invokes renderTableRow() method to render a row
  virtual int buildTable(ReportTableBase *table, PdfTable &pdfTable);
  // build a row common method. Cycles thu all columnType cells and renderCell for each ColumnType item
  virtual int buildTableRow(ReportTableBase *table, vector<string> *row, PdfTableRow  &pdfTableRow);
  // build cell common method. Renders a specific cell based on ColumnType specifications and a row of data
  virtual int buildTableCell(ReportColumnTypeBase *base, vector<string> *row, PdfTableCell &pdfTableCell);

  // Builders for the diferent cell types
  virtual int buildCellImageOnDemand(ReportColumnTypeImageOnDemand *, vector<string> *row, PdfTableCell &pdfTableCell);
  virtual int buildCellString(ReportColumnTypeString *, vector<string> *row, PdfTableCell &pdfTableCell);
  virtual int buildCellBox(ReportColumnTypeBox *, vector<string> *row, PdfTableCell &pdfTableCell);
  // Method to process 'multiple' table cells, to have several lines of data in the same table cell
  virtual int buildCellStringMultiple(ReportColumnTypeStringMultiple *ct, vector<string> *row, PdfTableCell &pdfTableCell);


  //////////////////////////////////////////////////////////////////////////////
  // Image renderers

  // method to generate images
  virtual int renderImage(const string &object, vector<string> & params);


  //////////////////////////////////////////////////////////////////////////////
  // Table exception renderers

  // exception for rendering main page
  virtual int renderTableExceptionMainPage(ReportTableBase *table, HPDF_Doc &pdf);
  // renderer for protein table header
  virtual int renderTableExceptionProteinHeader(ReportTableBase *table, HPDF_Doc &pdf);
  // helper for renderer for protein table header
  virtual void colorProteinString(vector<string> &in, string &out);
  // helper for renderer for protein table header
  virtual void breakProteinIntoChunks(vector<string> &in, int &count);


 public:


  // Constructors and destructor
  ReportRendererPdf()  {};
  ~ReportRendererPdf() {};

  virtual int initReport(HPDF_Doc &pdf);

  virtual int renderReport(ReportBase *, HPDF_Doc &pdf, bool renderDecoration = true, int start = 0, int count = -1);

  virtual int generateReport(ReportGeneratorData &data);

  // pagination
  virtual int paginate(ReportGeneratorData &data, ReportRendererBase &rr, ReportBase &rep);

};
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
