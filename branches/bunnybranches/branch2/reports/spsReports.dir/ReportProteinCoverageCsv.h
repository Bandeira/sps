///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_PROTEIN_COVERAGE_CSV_H__
#define __REPORT_PROTEIN_COVERAGE_CSV_H__
///////////////////////////////////////////////////////////////////////////////
#include "ReportBase.h"
#include "ReportTableProteinCoverage.h"

///////////////////////////////////////////////////////////////////////////////
// ReportProteinCoverage
//
// Tables are dynamically created, one for each sequence stretch
//
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
class ReportProteinCoverageCsv : public ReportBase {

 protected:


 public:

  // Constructors and destructor
  ReportProteinCoverageCsv(const string &projectPath, const string &proteinTableFilename)
  {
    ReportTableBase *is = new ReportTableProteinCoverage(projectPath, proteinTableFilename);
    // specify exclusive renderer method
    is->setRenderingException(RENDERING_EXCEPTION_PROTEIN_COVERAGE_CSV_PAGE);
    addTable(is);
  };

  virtual ~ReportProteinCoverageCsv() {};

};
///////////////////////////////////////////////////////////////////////////////
};
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
