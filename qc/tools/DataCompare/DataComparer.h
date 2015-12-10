///////////////////////////////////////////////////////////////////////////////
#ifndef __DATA_COMPARE_H__
#define __DATA_COMPARE_H__
///////////////////////////////////////////////////////////////////////////////
#include "SpecSet.h"
#include "abruijn.h"
#include "ReportTableBase.h"
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace spsReports;
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
class DataCompare {

  // spectrum to draw
  specnets::SpecSet m_specset[2];

  // ABinfo data set
  specnets::abinfo_t m_abinfo[2];

  // 2x int bin array
  vector<vector<int> > m_binArray[2];

  // report tables
  ReportTableBase m_reportTable[2];


  // the files format
  int m_type;

  // the results detail level
  int m_level;

  // compare two specsets
  int compareSpecsets(void);

  // compare two abruijn graphs
  int compareAbruijn(void);

  // compare two bin arrays
  int compareBinArray(void);

  // compare two report tables
  int compareReportTables(void);


 public:

  // constructors and destructor
  DataCompare()  {};
  ~DataCompare() {};

  // set the objects type
  void setType(int t) {m_type = t;};

  // set the objects type
  void setType(string &t);

  // set the detail level
  void setLevel(int l) {m_level = l;};

  // add an object
  int addObject(string filename, int index);

  // execute the comparison
  int exec(void);

};
///////////////////////////////////////////////////////////////////////////////
}; // namespace
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
