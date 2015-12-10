///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_PDF_H__
#define __REPORT_RENDERER_PDF_H__
///////////////////////////////////////////////////////////////////////////////
#include "ReportRendererBase.h"

///////////////////////////////////////////////////////////////////////////////
class ReportRendererPdf : public ReportRendererBase {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // Report page prolog and epilogue

  // prolog
  virtual int renderProlog(ReportBase *table, ostream &outstream);
  // epilog
  virtual int renderEpilog(ReportBase *table, ostream &outstream);
  // between tables
  virtual int renderInterTable(ReportBase *table, ostream &outstream);


  //////////////////////////////////////////////////////////////////////////////
  // Table header renderers

  // render a header row comon method. Cycles thu all columnType cells and renders each one.
  virtual int renderTableHeaderRow(ReportTableBase *table, ostream &outstream);
  // render header cell comon method. Renders all ColumnTypes
  virtual int renderTableHeaderCell(ReportColumnTypeBase *base, ostream &outstream);



  //////////////////////////////////////////////////////////////////////////////
  // Table content renderers

  // render table comon method. Cycles thru all the table rows and invokes renderTableRow() method to render a row
  virtual int renderTable(ReportTableBase *table, ostream &outstream);
  // render a row comon method. Cycles thu all columnType cells and renderCell for each ColumnType item
  virtual int renderTableRow(ReportTableBase *table, vector<string> *row, ostream &outstream);
  // render cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
  virtual int renderTableCell(ReportColumnTypeBase *base, vector<string> *row, ostream &outstream);

  // Renderers for the diferent cell types
  virtual int renderCellImageOnDemand(ReportColumnTypeImageOnDemand *)  {};
  virtual int renderCellString(ReportColumnTypeString *)                {};
  virtual int renderCellSequencesBox(ReportColumnTypeSequencesBox *)    {};

 public:

  // Constructors and destructor
  ReportRendererPdf()  {};
  ~ReportRendererPdf() {};

};
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
