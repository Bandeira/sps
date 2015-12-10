///////////////////////////////////////////////////////////////////////////////
#include "ReportContig.h"
#include "ReportTableContig.h"
#include "ReportTableClusterConsensus.h"
///////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportContig::ReportContig(const string &projectPath, const string &contigTableFilename)
{
  // Define contig list view, not filtered
  ReportTableBase *c = new ReportTableContig(projectPath, contigTableFilename);
  // Add the table to table list
  addTable(c);
}
////////////////////////////////////////////////////////////////////////////////
ReportContig::ReportContig(const string &projectPath, const string &contigTableFilename,const string &clusterTableFilename)
{
  // Add contig table for top contig object, filtered by contig
  ReportTableBase *c = new ReportTableContig(projectPath, contigTableFilename, TABLE_CONTIG_FILTER_COL_CONTIG);
  // Define view for top contig
  c->defineView2();
  // no borders and no header
  c->noBorders(); c->noHeaders();
  // Add table to table list
  addTable(c);

  // add cluster table for contig list, filtered by contig
  ReportTableBase *cc = new ReportTableClusterConsensus(projectPath, clusterTableFilename);
  // Add the table to table list
  addTable(cc);
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
///////////////////////////////////////////////////////////////////////////////
