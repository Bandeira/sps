///////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_CONTIG_H__
#define __REPORT_TABLE_CONTIG_H__
///////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>

#include "ReportTableBase.h"

///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// Specific defines
#define TABLE_CONTIG_FILTER_COL_CONTIG    0
#define TABLE_CONTIG_FILTER_COL_PROTEIN   1

///////////////////////////////////////////////////////////////////////////////
// Contig views
//
///////////////////////////////////////////////////////////////////////////////
// View for contig list
////////////////////////////////////////
//
// Table colTypes has the following structure:
//
// colTypes[0] -> (CTstring) Contig index
//   --> text = index in specset = "<0>"
//   --> columnLabel = "Index"
// colTypes[1] -> (CTstring) Number of spectra
//   --> text = number of spectra for this contig = "<1>"
//   --> columnLabel = "Spectra"
// colTypes[2] -> (CTIOD) Contig image
//   --> icon  --> "/cgi/contplot --projectdir <projectdir> --p <paramsfile> --contig <0> --reference <2> --homolog <3> --consensus <5> --user <5> --zoom .4 --target cout --spectra<8>"
//   --> columnLabel = "Contig"
//   --> url = "/cgi/spsplot --projectdir <projectdir> --p <paramsfile> --table contig --contig <0> --output cout"
//   --> label = ""
// colTypes[2] -> (CTseqsBox)
//   --> columnLabel = "Sequences"
//     sequences[0]->(CTIOD) : Reference
//       --> icon  --> ""
//       --> label --> Reference sequence = "<2>"
//       --> url   --> ""
//     sequences[1]->(CTIOD) : Homolog
//       --> icon  --> ""
//       --> label --> Homolog sequence = "<3>"
//       --> url   --> ""
//     sequences[2]->(CTIOD) : consensus
//       --> icon  --> ""
//       --> label --> consensus sequence = "<4>"
//       --> url   --> ""
//     sequences[3]->(CTIOD) : User
//       --> icon  --> ""
//       --> label --> User sequence = "<5>"
//       --> url   --> ""
//       --> id    --> 'user_<row>_<col>'
//     sequences[4]-> isInput  = true
//       --> text = ""
//       --> id = 'user_<row>_<col>'
//     sequences[5]-> isButton = true
//       --> text = "Update"
//       --> onClick="javascript:DoOnCick('input_<row>_<col>', <0>, 'user_<row>_<col>')";
//      #--> url --> --update <0> --pklbin <1> --peptide <6>
// colTypes[3] -> (CTstring) protein
//       --> text = protein name & description = "<6><nl><7>"  # <nl> means 'newline'
//       --> columnLabel = "Protein"
//
///////////////////////////////////////////////////////////////////////////////
// View for contig image with edit box
////////////////////////////////////////
//
// Table colTypes has the following structure:
// colTypes[0] -> (CTseqsBox)
//   --> columnLabel = ""
//     sequences[0]->(CTIOD) : Reference
//       --> icon  --> "/cgi/contplot --projectdir <projectdir> --contig <0> --reference <1> --homolog <2> --consensus <3> --user <4> --zoom 1 --output-format uu64 --output cout"
//       --> columnLabel = ""
//       --> url = ""
//       --> label = ""
//     sequences[1]-> isInput  = true
//       --> text = ""
//       --> id    --> 'user_<row>_<col>'
//     sequences[2]-> isButton = true
//       --> text = "Update"
//       --> onClick='javascript:DoOnCick('input_<row>_<col>', <0>, user_<row>_<col>);';
//       --> url --> --update --contig <0> --peptide <4>
//
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
//
// cells[row][0] -> text   --> contig index
// cells[row][1] -> text   --> number of spectra
// cells[row][2] -> text   --> Reference sequence --> generateSequence(i, return)
// cells[row][3] -> text   --> Homolog sequence   --> generateSequence(i, return)
// cells[row][4] -> text   --> Consensus sequence --> generateSequence(i, return)
// cells[row][5] -> text   --> User sequence
// cells[row][6] -> text   --> Protein name
// cells[row][7] -> text   --> Protein description
// cells[row][8] -> text   --> spectra (coma separated indexes)
// cells[row][9] -> text   --> protein index
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class ReportTableContig : public ReportTableBase {

 protected:

 public:

  //////////////////////////////////////////////////////////////////////////////
  // Constructors and destructor

  ReportTableContig(const string &projectPath, const string &tableFilename, int columnFilter = TABLE_FILTER_COL_NONE);

  virtual ~ReportTableContig() {};

  // element find methods
  virtual int find(string &data);
  virtual int find(int data);

  // IDs[0] = proteinID, IDs[1] = contigID
  virtual void getId(vector<string> &IDs) {getFilteredId(IDs, 1);ReportTableBase::getId(IDs);};

  //////////////////////////////////////////////////////////////////////////////
  // Methods to build views

  //default view
  virtual void defineView(void);
  // use a different view
  virtual void defineView2(void);
  // view to generate images
  virtual void defineViewImages(void);

};
////////////////////////////////////////////////////////////////////////////////
}; // namespace
///////////////////////////////////////////////////////////////////////////////
#endif
///////////////////////////////////////////////////////////////////////////////
