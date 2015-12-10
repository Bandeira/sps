///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_SERVER_INTERFACE_H__
#define __REPORT_SERVER_INTERFACE_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "ParameterList.h"
#include "ReportTableBase.h"

#include "ReportData.h"

///////////////////////////////////////////////////////////////////////////////
using namespace std;
///////////////////////////////////////////////////////////////////////////////
namespace spsReports {
///////////////////////////////////////////////////////////////////////////////
// Defines
#define STR(s) #s
#define XSTR(s) STR(s)

#define DEFAULT_ROWS_PER_TABLE  20
///////////////////////////////////////////////////////////////////////////////
class ReportServerInterface {

  // tables names
  string m_tableNameHeader;
  string m_tableNameProtein;
  string m_tableNameProteinCoverage;
  string m_tableNameContig;
  string m_tableNameCluster;
  string m_tableNameSpectra;

  // build directory path by concatenation with given path
  string composeFileName(const string &projectDir, const string &fileName);

  // get a table object
  ReportTableBase *getTableObject(ReportData &data);


 public:


  // Constructors and destructor
  ReportServerInterface();
  ~ReportServerInterface();

  // Option parsing
  int parseOptions(int argc, char **argv);
  // Option processing
  int processOptions(specnets::ParameterList &commandLineParams);

  // Output help
  int help(ostream &);

  // Output version information
  int version(ostream &);

  // Output error messages
  int error(const string &);

  // Spectrum drawing method
  int updateTable(ReportData &data);

  // Spectrum drawing method
  int getTableData(ReportData &data);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
