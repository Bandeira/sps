////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_HTML_SINGLE_H__
#define __REPORT_RENDERER_HTML_SINGLE_H__
////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererHtml.h"



////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// ReportRendererHtmlSingle
//
//
////////////////////////////////////////////////////////////////////////////////
class ReportRendererHtmlSingle : public ReportRendererHtml {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // Reports pages rendering methods (specific pages)

  virtual int renderTable(ReportTableBase *table, ostream &outstream);

  virtual int renderTableRow(vector<string> &row, ostream &outstream);

  //////////////////////////////////////////////////////////////////////////////
  // Report page prolog and epilogue

  // prolog
  virtual int renderProlog(ostream &outstream);
  virtual int renderProlog2(ostream &outstream);
  // epilog
  virtual int renderEpilog(ostream &outstream);


  //////////////////////////////////////////////////////////////////////////////
  // Table header builders

  //////////////////////////////////////////////////////////////////////////////
  // Table content builders

  // main page
  virtual int renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream);

  // table cell renderer
  virtual int buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &ss);

  // Builders for the diferent cell types
  virtual int buildCellImageOnDemand(ReportColumnTypeImageOnDemand *, vector<string> *row, stringstream &ss);

  // write the table
  virtual int writeTable(ReportTableBase *table, ostream &outstream, int idx);


 public:


  // Constructors and destructor
  ReportRendererHtmlSingle()  {};
  ~ReportRendererHtmlSingle() {};

  virtual int  generateReport(ReportGeneratorData &data);

};
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
