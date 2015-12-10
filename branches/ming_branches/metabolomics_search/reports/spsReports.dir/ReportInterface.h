///////////////////////////////////////////////////////////////////////////////
#ifndef __SPS_PLOT_H__
#define __SPS_PLOT_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>

#include "ReportTableGenerator.h"
#include "ParameterList.h"
#include "ReportBase.h"
#include "ReportRendererBase.h"


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
class ReportInterface {

  // tables names
  string m_tableNameHeader;
  string m_tableNameProtein;
  string m_tableNameProteinCoverage;
  string m_tableNameContig;
  string m_tableNameCluster;
  string m_tableNameSpectra;

  // verbose flag
  bool m_verbose;


  // build directory path by concatenation with given path
  string composeFileName(const string &projectDir, const string &fileName);

  // build directory path by concatenation with path
  int buildDirectoryPath(const string &dir);

  // pagination
  int paginate(ReportGeneratorData &data, ReportRendererBase &rr, ReportBase &rep, string fnamePrefix, string pageName);

  // verbose output
  void verboseOutput(ostream &os, const double &v, bool nl = true);
  void verboseOutput(ostream &os, const char *str, bool nl = true);
  void verboseOutput(ostream &os, const char *str, const char *val, bool nl = true);
  void verboseOutput(ostream &os, const char *str, const char *val, const char *term, bool nl = true);

  void dump_abruijn(ReportGeneratorData &data);
  void dump_binArray(ReportGeneratorData &data);


 public:


  // Constructors and destructor
  ReportInterface();
  ~ReportInterface();

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
  int buildTables(ReportGeneratorData &data);

  // Generate HTML report pages
  int generateReportHtml(ReportGeneratorData &data);

  // Generate dynamic HTML report entry page
  int generateReportHtmlDynamic(ReportGeneratorData &data);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
