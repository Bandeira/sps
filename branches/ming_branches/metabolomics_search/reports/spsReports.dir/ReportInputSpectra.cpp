///////////////////////////////////////////////////////////////////////////////
#include "ReportInputSpectra.h"
#include "ReportTableInputSpectra.h"

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportInputSpectra::ReportInputSpectra(const string &projectPath, const string &inputSpectraTableFilename)
{
  // Define the table, set the filter column to filter by file index
  ReportTableBase *is = new ReportTableInputSpectra(projectPath, inputSpectraTableFilename, TABLE_SPECTRA_FILTER_COL_FILE);
  // add the table to the report
  addTable(is);
}
///////////////////////////////////////////////////////////////////////////////
}; //namespace
///////////////////////////////////////////////////////////////////////////////
