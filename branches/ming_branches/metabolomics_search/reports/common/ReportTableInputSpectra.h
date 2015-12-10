////////////////////////////////////////////////////////////////////////////////
#ifndef __REPORT_TABLE_INPUT_SPECTRA_H__
#define __REPORT_TABLE_INPUT_SPECTRA_H__
////////////////////////////////////////////////////////////////////////////////
// includes
#include <string>
#include <vector>

#include "spectrum.h"
#include "ReportTableBase.h"

////////////////////////////////////////////////////////////////////////////////
// defines
#define TABLE_INPUT_SPECTRA_FILENAME "tableInputSpectra.txt"

#define TABLE_SPECTRA_FILTER_COL_CLUSTER  3
#define TABLE_SPECTRA_FILTER_COL_FILE     5

namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
// View for input spectra list
////////////////////////////////////////
//
// Table colTypes has the following structure:
//
// colTypes[0] -> (CTstring) Spectrum index / spectrum scan
//   --> text = index in specset = "<1>"
//   --> columnLabel = "Scan"
// colTypes[1] -> (CTstring) Spectrum file name
//   --> text = file name where the spectrum is contained = "<6>"
//   --> columnLabel = "File name"
// colTypes[2] -> (CTseqsBox)
//   --> columnLabel = "Sequences"
//     sequences[0]->(CTIOD) : Reference
//       --> icon  --> ""
//       --> label --> Homolog sequence = "<7>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrumscan <2> --<16> <1> --peptide <7> --target cout
//     sequences[1]->(CTIOD) : Homolog
//       --> icon  --> ""
//       --> label --> Homolog sequence = "<8>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrumscan <2> --<16> <1> --peptide <8> --target cout
//     sequences[2]->(CTIOD) : deNovo
//       --> icon  --> ""
//       --> label --> deNovo sequence = "<9>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrumscan <2> --<16> <1> --peptide <9> --target cout
//     sequences[3]->(CTIOD) : User
//       --> icon  --> ""
//       --> label --> User sequence = "<10>"
//      #--> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrumscan <0> --<16> <1> --peptide <10> --target cout
//       --> onClick = 'javascript:DoOnCick('user_<row>_<col>, <2>);';
//       --> id = 'user_<row>_<col>' --> <row> is row #; <col> is column #
//     sequences[4]-> isInput  = true
//       --> text = ""
//       --> id = 'input_<row>_<col>'
//     sequences[5]-> isButton = true ;
//       --> text = "Update"
//       --> onClick = 'javascript:DoOnCick('input_<row>_<col>, <0>, user_<row>_<col>);';
//
// possible replacement
// -p <paramsfile> =
//
//
// colTypes[3] -> (CTstring) mass
//   --> text = value in specset[spectrum index] = "<11>"
//   --> columnLabel = "Mass"
// colTypes[4] -> (CTstring) charge
//   --> text = vaue in specset[spectrum index] = "<12>"
//   --> columnLabel = "Charge"
// colTypes[5] -> (CTstring) B%
//   --> text = "<13>"
//   --> columnLabel = "B %"
// colTypes[6] -> (CTstring) Y%
//   --> text = "<14>"
//   --> columnLabel = "Y %"
// colTypes[7] -> (CTstring) BY intensity %
//   --> text = "<15>"
//   --> columnLabel = "BY int %"
//
////////////////////////////////////////////////////////////////////////////////
// View for input spectra list (static files)
////////////////////////////////////////
//
// colTypes[0] -> (CTstring) Spectrum index
//   --> text = index in specset = "<1>"
//   --> columnLabel = "Scan"
// colTypes[1] -> (CTstring) Spectrum file name
//   --> text = file name where the spectrum is contained = "<6>"
//   --> columnLabel = "File name"
// colTypes[2] -> (CTseqsBox)
//   --> columnLabel = "Sequences"
//     sequences[0]->(CTIOD) : Reference
//       --> icon  --> ""
//       --> label --> Homolog sequence = "<7>"
//       #--> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrumscan <2> --<16> <1> --peptide <7> --target cout
//       --> params --> --projectdir <projectdir> --p <paramsfile> --spectrumscan <2> --<16> <1> --peptide <7> --target cout
//       --> renderer --> specplot
//     sequences[1]->(CTIOD) : Homolog
//       --> icon  --> ""
//       --> label --> Homolog sequence = "<8>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrumscan <2> --<16> <1> --peptide <8> --target cout
//     sequences[2]->(CTIOD) : deNovo
//       --> icon  --> ""
//       --> label --> deNovo sequence = "<9>"
//       --> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrumscan <2> --<16> <1> --peptide <9> --target cout
//     sequences[3]->(CTIOD) : User
//       --> icon  --> ""
//       --> label --> User sequence = "<10>"
//      #--> url   --> /cgi-bin/specplot --projectdir <projectdir> --p <paramsfile> --spectrumscan <0> --<16> <1> --peptide <10> --target cout
//       --> onClick = 'javascript:DoOnCick('user_<row>_<col>, <2>);';
//       --> id = 'user_<row>_<col>' --> <row> is row #; <col> is column #
//     sequences[4]-> isInput  = true
//       --> text = ""
//       --> id = 'input_<row>_<col>'
//     sequences[5]-> isButton = true ;
//       --> text = "Update"
//       --> onClick = 'javascript:DoOnCick('input_<row>_<col>, <0>, user_<row>_<col>);';
//
// possible replacement
// -p <paramsfile> =
//
//
// colTypes[3] -> (CTstring) mass
//   --> text = value in specset[spectrum index] = "<11>"
//   --> columnLabel = "Mass"
// colTypes[4] -> (CTstring) charge
//   --> text = vaue in specset[spectrum index] = "<12>"
//   --> columnLabel = "Charge"
// colTypes[5] -> (CTstring) B%
//   --> text = "<13>"
//   --> columnLabel = "B %"
// colTypes[6] -> (CTstring) Y%
//   --> text = "<14>"
//   --> columnLabel = "Y %"
// colTypes[7] -> (CTstring) BY intensity %
//   --> text = "<15>"
//   --> columnLabel = "BY int %"
//
////////////////////////////////////////////////////////////////////////////////
// Generated table, per row:
////////////////////////////////////////
//
// cells[row][0] -> text   --> index in specset
// cells[row][1] -> text   --> index in specset
// cells[row][2] -> text   --> scan #
// cells[row][3] -> text   --> cluster index (used when filtered by cluster consensus)
// cells[row][4] -> text   --> Protein index
// cells[row][5] -> text   --> Protein name
// cells[row][6] -> text   --> spectrum file index
// cells[row][7] -> text   --> spectrum file name
// cells[row][8] -> text   --> Reference sequence --> generateSequence()
// cells[row][9] -> text   --> Homolog sequence --> generateSequence()
// cells[row][10] -> text  --> DeNovo sequence --> generateSequence()
// cells[row][11] -> text  --> User sequence
// cells[row][12] -> text  --> mass value, from specs[i].parentMass
// cells[row][13] -> text  --> charge value, from specs[i].parentCharge
// cells[row][14] -> text  --> B%
// cells[row][15] -> text  --> Y%
// cells[row][16] -> text  --> BY Int %
// cells[row][17] -> text  --> input file extension: mgf or mzxml
//
////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////
class ReportTableInputSpectra : public ReportTableBase {

 protected:

  //////////////////////////////////////////////////////////////////////////////
  // Sequence generation methods

  // Generate a DeNovo sequence given a spectrum index
  // int generateSequenceDeNovo(int spectrumIdx, string &sequence);
  // Generate a Homolog sequence given a spectrum index
  // int generateSequenceHomolog(int spectrumIdx, string &sequence);
  // Generate a Reference sequence given a spectrum index
  // int generateSequenceReference(int spectrumIdx, string &sequence);


  //////////////////////////////////////////////////////////////////////////////
  // Helper table generation methods

  // Calculate B percentage
  //int calculateBper();
  // Calculate Y percentage
  //int calculateYper();
  // Calculate BY intensity percentage
  //int calculateBYintPer();


 public:


  //////////////////////////////////////////////////////////////////////////////
  // Constructors and destructor

  ReportTableInputSpectra(const string &projectPath, const string &tableFilename, int columnFilter = TABLE_SPECTRA_FILTER_COL_CLUSTER);

  virtual ~ReportTableInputSpectra() {};


  //////////////////////////////////////////////////////////////////////////////
  // Methods to build views

  // default view (input spectra page)
  virtual void defineView(void);
  // cluster page view
  virtual void defineView2(void);
  // populate IDs column
  void populateIDs(void);
};
////////////////////////////////////////////////////////////////////////////////
}; //namespace
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
