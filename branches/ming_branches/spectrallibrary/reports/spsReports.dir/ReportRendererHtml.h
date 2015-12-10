////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_HTML_H__
#define __REPORT_RENDERER_HTML_H__
////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererBase.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// ReportRendererHtml
//
// According to each columnType, a specific HTML sequence is output to output stream.
//
// While building column headers
//
//
// <td>
//   <columnLabel>
// </td>
//
//
// -----------------------------------------------------------------------------
// While building table
//
// ReportColumnTypeString:
// <td class="<cssClass>">
//   <text> or
//   <input type="button" value="<text>" onClick="<onClick>" /> or
//   <input ID="<ID>" type="text" style="text-transform: uppercase; width:100%" />
// </td>
//
//
// ReportColumnTypeImageOnDemand:
// <td class="<cssClass>">
//   <a href="<>" onclick="<onClick>">
//       <label>  or <img src='<icon>'>
//   </a>
// </td>
//
//
// ReportColumnTypeSequencesBox:
// <td class="<cssClass>">
//   <table>
//                --- begining of repeat block
//     <tr>
//       call to ReportColumnType
//     </tr>
//                --- repeat per ReportColumnType entry on vector
//   </table>
// </td>
//
////////////////////////////////////////////////////////////////////////////////
class ReportRendererHtml : public ReportRendererBase {

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
  // Table content renderers

  // render table comon method. Cycles thru all the table rows and invokes buildTableRow() method to render a row
  virtual int renderTable(ReportTableBase *table, ostream &outstream);
  // render a header row method. Renders a table header row.
  virtual int renderTableHeaderRow(vector<string> &row, ostream &outstream);
  // render header cell method. Renders a header cell
  //virtual int renderTableHeaderCell(string &cell, ostream &outstream);
  // build a row comon method. Cycles thu all columnType cells and buildCell for each ColumnType item
  virtual int renderTableRow(vector<string> &row, ostream &outstream);
  // build cell comon method. builds a specific cell based on ColumnType specifications and a row of data
  //virtual int renderTableCell(string &cell, ostream &outstream);


  //////////////////////////////////////////////////////////////////////////////
  // Table header builders

  // render a header row comon method. Cycles thu all columnType cells and renders each one.
  //virtual int buildTableHeaderRow(ReportTableBase *table, vector<string> &renderedRow);
  // render header cell comon method. Renders all ColumnTypes
  virtual int buildTableHeaderCell(ReportColumnTypeBase *base, stringstream &ss);


  //////////////////////////////////////////////////////////////////////////////
  // Table content builders

  // build table comon method. Cycles thru all the table rows and invokes renderTableRow() method to render a row
  //virtual int buildTable(ReportTableBase *table);
  // build a row comon method. Cycles thu all columnType cells and renderCell for each ColumnType item
  //virtual int buildTableRow(ReportTableBase *table, vector<string> *row, vector<string> &renderedRow);
  // build cell comon method. Renders a specific cell based on ColumnType specifications and a row of data
  virtual int buildTableCell(ReportColumnTypeBase *base, vector<string> *row, stringstream &ss);

  // Builders for the diferent cell types
  virtual int buildCellImageOnDemand(ReportColumnTypeImageOnDemand *, vector<string> *row, stringstream &ss);
  virtual int buildCellString(ReportColumnTypeString *, vector<string> *row, stringstream &ss);
  virtual int buildCellBox(ReportColumnTypeBox *, vector<string> *row, stringstream &ss);
  // Method to process 'multiple' table cells, to have several lines of data in the same table cell
  virtual int buildCellStringMultiple(ReportColumnTypeStringMultiple *ct, vector<string> *row, stringstream &ss);


  //////////////////////////////////////////////////////////////////////////////
  // Image renderers

  // method to generate images
  virtual int renderImage(const string &object, vector<string> & params);


  //////////////////////////////////////////////////////////////////////////////
  // Table exception renderers

  // exception for rendering main page
  virtual int renderTableExceptionMainPage(ReportTableBase *table, ostream &outstream);
  // renderer for protein table header
  virtual int renderTableExceptionProteinHeader(ReportTableBase *table, ostream &outstream);
  // helper for renderer for protein table header
  virtual void colorProteinString(vector<string> &in, string &out);
  // helper for renderer for protein table header
  virtual void breakProteinIntoChunks(vector<string> &in, int &count);


 public:


  // Constructors and destructor
  ReportRendererHtml()  {};
  ~ReportRendererHtml() {};

  virtual void buildNavigationBar(string &navBar, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames);

};
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
