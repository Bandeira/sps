////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_CLUSTER_CONSENSUS_H__
#define __REPORT_TABLE_CLUSTER_CONSENSUS_H__
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>

#include "spectrum.h"
#include "ReportTableBase.h"
#include "ClusterData.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// Specific defines
#define TABLE_CLUSTER_FILTER_COL_CLUSTER  0
#define TABLE_CLUSTER_FILTER_COL_CONTIG   1
////////////////////////////////////////////////////////////////////////////////
// Cluster consensus views
//
////////////////////////////////////////////////////////////////////////////////
// View for cluster consensus list
////////////////////////////////////////
//
// Table colTypes has the following structure:
//
// colTypes[0] -> (CTstring) Spectrum index
//   --> text = index in specset = "<0>"
//   --> columnLabel = "Scan"
// colTypes[1] -> (CTseqsBox)
//   --> columnLabel = "Sequences"
//     sequences[0]->(CTIOD) : Reference
//       --> icon  --> ""
//       --> label --> Reference sequence = "<3>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrum <0> --pklbin <12> --peptide <3>
//     sequences[1]->(CTIOD) : Homolog
//       --> icon  --> ""
//       --> label --> Homolog sequence = "<4>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrum <0> --pklbin <12> --peptide <4>
//     sequences[2]->(CTIOD) : deNovo
//       --> icon  --> ""
//       --> label --> deNovo sequence = "<5>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrum <0> --pklbin <12> --peptide <5>
//     sequences[3]->(CTIOD) : User
//       --> icon  --> ""
//       --> label --> User sequence = "<6>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrum <0> --pklbin <12> --peptide <6>
//       --> onClick = 'javascript:DoOnCick('user_<row>_<col>, <0>);';
//       --> id = 'user_<row>_<col>' --> <row> is row #; <col> is column #
//     sequences[4]-> isInput  = true
//       --> text = ""
//       --> id = 'input_<row>_<col>'
//     sequences[5]-> isButton = true
//       --> text = "Update"
//       --> onClick='javascript:DoOnCick('input_<row>_<col>, <0>, user_<row>_<col>);';
//
// colTypes[2] -> (CTstring) mass
//   --> text = value in specset[spectrum index] = "<7>"
//   --> columnLabel = "Mass"
// colTypes[3] -> (CTstring) charge
//   --> text = value in specset[spectrum index] = "<8>"
//   --> columnLabel = "Charge"
// colTypes[4] -> (CTstring) B%
//   --> text = "<9>"
//   --> columnLabel = "B %"
// colTypes[5] -> (CTstring) Y%
//   --> text = "<10>"
//   --> columnLabel = "Y %"
// colTypes[6] -> (CTstring) BY intensity %
//   --> text = "<11>"
//   --> columnLabel = "BY int %"
//
///////////////////////////////////////////////////////////////////////////////
// View for cluster consensus image and input box with update button
//
// Table colTypes has the following structure:
// colTypes[0] -> (CTseqsBox)
//   --> columnLabel = ""
//     sequences[0]->(CTIOD) : Reference
//       --> icon  --> "/cgi/specplot --projectdir <projectdir> --p <paramsFile> --spectrum <0> --peptide <1> --zoom 1 --output-format uu64 --output cout"
//       --> columnLabel = ""
//       --> url = ""
//       --> label = ""
//     sequences[4]-> isInput  = true
//       --> text = ""
//       --> id = 'input_<row>_<col>'
//     sequences[5]-> isButton = true
//       --> text = "Update"
//       --> onClick='javascript:DoOnCick('input_<row>_<col>', <0>, user_<row>_<col>);';
//       --> url --> --update --spectra/consensus <0> --peptide <1>
//////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text   --> index in specset (consensus spectra index)
// cells[row][1] -> text   --> contig ID
// cells[row][2] -> text   --> Reference sequence --> generateSequence()
// cells[row][3] -> text   --> Homolog sequence --> generateSequence()
// cells[row][4] -> text   --> DeNovo sequence --> generateSequence()
// cells[row][5] -> text   --> User sequence
// cells[row][6] -> text   --> mass value, from specs[i].parentMass
// cells[row][7] -> text   --> charge value, from specs[i].parentCharge
// cells[row][8] -> text   --> B%
// cells[row][9] -> text   --> Y&
// cells[row][10] -> text  --> BY Int %
// cells[row][11] -> text  --> protein name
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class ReportTableClusterConsensus : public ReportTableBase {

 protected:



 public:

  //////////////////////////////////////////////////////////////////////////////
  // Constructors and destructor

  ReportTableClusterConsensus(const string &projectPath, const string &tableFilename, int columnFilter = TABLE_CLUSTER_FILTER_COL_CONTIG);

  virtual ~ReportTableClusterConsensus() {};


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build views

  // default view
  virtual void defineView(void);
  // use a different view
  virtual void defineView2(void);

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
