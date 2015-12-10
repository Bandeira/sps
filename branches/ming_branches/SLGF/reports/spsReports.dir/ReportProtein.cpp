///////////////////////////////////////////////////////////////////////////////
#include "ReportProtein.h"
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#include "ReportTableProtein.h"
#include "ReportTableContig.h"

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportProtein::ReportProtein(const string &projectPath, const string &proteinTableFilename)
{
  ReportTableBase *is = new ReportTableProtein(projectPath, proteinTableFilename);
  addTable(is);
}
////////////////////////////////////////////////////////////////////////////////
ReportProtein::ReportProtein(const string &projectPath, const string &proteinTableFilename, const string &contigTableFilename)
{
  // Add contig table for top contig object, filtered by contig
  ReportTableBase *c = new ReportTableProtein(projectPath, proteinTableFilename, TABLE_PROTEIN_FILTER_COL_PROTEIN);
  // Define view for top contig
  c->defineView2();
  // no borders and no header
  c->noBorders(); c->noHeaders();
  // specify exclusive renderer method
  c->setRenderingException(RENDERING_EXCEPTION_PROTEIN_HEADER);
  // Add table to table list
  addTable(c);



  // add cluster table for contig list, filtered by contig
  ReportTableBase *cc = new ReportTableContig(projectPath, contigTableFilename, TABLE_CONTIG_FILTER_COL_PROTEIN);
  // Add the table to table list
  addTable(cc);
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////