////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_RENDERER_HTML_H__
#define __REPORT_RENDERER_HTML_H__
////////////////////////////////////////////////////////////////////////////////
#include "ReportRendererBase.h"


#define DEFAULT_ROWS_PER_TABLE  20

#define TABLE_SEP_L1  "|"
#define TABLE_SEP_L2  "@"
#define TABLE_SEP_L3  "&"
#define TABLE_SEP_L4  "!"

#define TABLE_COVERAGE_FIELD_ID               0
#define TABLE_COVERAGE_FIELD_NAME             1
#define TABLE_COVERAGE_FIELD_SEQ_REFERENCE    2
#define TABLE_COVERAGE_FIELD_PROT_SEQUENCE    3
#define TABLE_COVERAGE_CSPS_DATA              4
#define TABLE_COVERAGE_SPS_DATA               5


#define CELLS_PER_LINE 20


#define IMAGE_ICON_ID_PREFIX        "im_"
#define IMAGE_ICON_CTRL_ID_PREFIX   "imc_"
#define IMAGE_LARGE_ID_PREFIX       "io_"
#define IMAGE_LARGE_CTRL_ID_PREFIX  "ioc_"


////////////////////////////////////////////////////////////////////////////////
// Holds a contig/protein mass index pair unit
struct contigMatchElem {
  int     start;
  int     colSpan;
  string  data;
};

////////////////////////////////////////////////////////////////////////////////
//
struct ContigData {

  int                     id;
  string                  name;
  int                     start;
  int                     end;
  vector<contigMatchElem> contigMatch;

};
////////////////////////////////////////////////////////////////////////////////
// Ordering vector element. Used to estabelish contig order when rendering protein coverage report
struct ContigOrdering {
  int contigIndex;
  int startIndex;
  int endIndex;

  // Used by sort method. First order by beggining index, then by ending index.
  bool operator<(const ContigOrdering &o) const
  {return (startIndex == o.startIndex ? endIndex < o.endIndex : startIndex < o.startIndex);};

};

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

  // variable to keep track of the number of lauched processes
  int m_nprocess;

  // generated file number, used to support multi processing image generation
  int m_fileNumber;

  //////////////////////////////////////////////////////////////////////////////
  // Reports pages rendering methods (specific pages)

  // Header page generation
  virtual int  generateReportPageHeader(ReportGeneratorData &data);
  // Render protein list page
  virtual int  generateReportPageProteins(ReportGeneratorData &data);
  // Render single protein page
  virtual int  generateReportPageProtein(ReportGeneratorData &data);
  // Render the protein coverage pages
  virtual int  generateReportPageProteinCoverage(ReportGeneratorData &data);
  // Render the protein coverage CSV text files
  virtual int  generateReportPageProteinCoverageCSV(ReportGeneratorData &data);
  // Render contig list page
  virtual int  generateReportPageContigs(ReportGeneratorData &data);
  // Render single contig page
  virtual int  generateReportPageContig(ReportGeneratorData &data);
  // Render single cluster page
  virtual int  generateReportPageCluster(ReportGeneratorData &data);
  // render spectra page
  virtual int  generateReportPageSpectra(ReportGeneratorData &data);

  //////////////////////////////////////////////////////////////////////////////
  // Report rendering general methods

  // paginnation method used to split long pages in smaller ones
  virtual int paginate(ReportGeneratorData &data, ReportBase &rep, int ipp, string fnamePrefix, string pageName, int barType);
  // page building routine in multi-process
  //virtual int buildPage(ReportBase &rep, string &fn, int start, int len);
  virtual int buildPage(ReportBase &rep, string &fn, int start, int len, string &navBarTop);

  // build the navigation bar
  virtual int buildNavigationBar(ReportBase &rep, stringstream &out, int barType);
  // add an element to the navigation bar
  virtual int addBarLink(stringstream &out, string &img, string &text, string &link);


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

  // build table common method. Cycles thru all the table rows and invokes renderTableRow() method to render a row
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
  virtual int renderImage(const string &object, vector<string> &params, vector<ReportParamsFiles> &files);


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

  virtual void genRandomString(string &s, const int len);


  //////////////////////////////////////////////////////////////////////////////
  // Protein coverage renderer methods

  //
  int  renderTableExceptionProteinCoverage(ReportTableBase *table, ostream &outstream);
  //
  int  renderTableExceptionProteinCoverageCSV(ReportTableBase *table, ostream &outstream);
  //
  void generateProteinSequence(ostream &outstream, int i, vector<string> &proteinData, int proteinLength);
  //
  void generateProteinSequenceCSV(ostream &outstream, int i, vector<string> &proteinData, int proteinLength);
  //
  void generateOutputContig(int i, int proteinSize, vector<int> &vectorID, ostream &outstream, int cellPerLine, vector<ContigData> &contig, bool link);
  //
  void generateOutputContigCSV(int i, int proteinSize, vector<int> &vectorID, ostream &outstream, int cellPerLine, vector<ContigData> &contig, bool link);
  //
  void getOrder(vector<ContigData> &contig, int i, int size, vector<int> &order);
  //
  string getIntFromSeqName(const string &seq);
  //
  void buildContigDataStructure(vector<ContigData> &contigDataArray, string &contig);


  void parseFilesVector(vector<string> *row, vector<ReportParamsFiles> &files);



 public:


  // Constructors and destructor
  ReportRendererHtml() : m_nprocess(0), m_fileNumber(0) {};
  ~ReportRendererHtml() {};

  virtual int  generateReport(ReportGeneratorData &data);

  virtual void buildNavigationBar(string &navBar, vector<string> &IDs, vector<string> &IDsEnd, vector<string> &fNames);

};
///////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
