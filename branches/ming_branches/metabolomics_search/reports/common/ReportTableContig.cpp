////////////////////////////////////////////////////////////////////////////////
#include "ReportTableContig.h"

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////
ReportTableContig::ReportTableContig(const string &projectPath, const string &tableFilename, int columnFilter)
: ReportTableBase(projectPath, tableFilename)
{
  // set the filter column
  setFilterColumn(columnFilter);
  // define view (cluster list)
  defineView();
  // set sort column
  setIdColumn(0);
}
////////////////////////////////////////////////////////////////////////////////
void ReportTableContig::defineView(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;

  ///////////////////////////////////////////////////////////////////////////////
  // View for contig list
  ////////////////////////////////////////
  // Table colTypes has the following structure:

  // colTypes[0] -> (CTstring) Contig index
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Index";
  auxS->text         = "<0>"; //"�0�";
  m_colTypes.push_back(auxS);

  // colTypes[1] -> (CTstring) Number of spectra
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Spectra";
  auxS->text         = "<2>";
  m_colTypes.push_back(auxS);

  // colTypes[2] -> (CTIOD) Contig image
  auxI = new ReportColumnTypeImageOnDemand();
  auxI->link              = "contig.<0>.0.html";
  auxI->columnLabel       = "Contig";
  auxI->iconRenderer      = "contplot";

  auxI->iconParams.push_back( ReportParamsOption("--star",            "<16>",    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--abruijn",         "<15>",    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--seqs",            "<17>",    ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--star",            "<projectdir>/<16>",    ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--abruijn",         "<projectdir>/<15>",    ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--seqs",            "<projectdir>/<17>",    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--contig",          "<0>",                                      "<0>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--reference",       "<3>",                                      "<3>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--homolog",         "<4>",                                      "<4>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--denovo",          "<5>",                                      "<5>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--user",            "<6>",                                      "<6>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--mass-reference",  "<10>",                                     "<10>") );
  auxI->iconParams.push_back( ReportParamsOption("--mass-homolog",    "<11>",                                     "<11>") );
  auxI->iconParams.push_back( ReportParamsOption("--offset-reference","<12>",                                     "<12>") );
  auxI->iconParams.push_back( ReportParamsOption("--offset-homolog",  "<13>",                                     "<13>") );
  auxI->iconParams.push_back( ReportParamsOption("--reverse",         "",                                         "<14>") );
  auxI->iconParams.push_back( ReportParamsOption("--target",          "internal",                                 ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",        "uu64",                                     ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--notitle",         "",                                         ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--zoom",            "0.4",                                      ""    ) );
  //auxI->iconParams        = "--star spectra/stars.pklbin --abruijn assembly/component_info.bin --seqs assembly/sps_seqs.pklbin --contig <0> --reference \"<3>\" --homolog \"<4>\" --denovo \"<5>\" --user \"<6>\" --mass-reference \"<10>\" --mass-homolog \"<11>\"  --target internal --encoding uu64 --zoom 0.35";
  //auxI->iconParams         = "contigs.<0>.png";
  m_colTypes.push_back(auxI);


  // colTypes[2] -> (CTseqsBox)
  auxB = new ReportColumnTypeBox();
  auxB->columnLabel  = "Contig sequence";

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel  = "Reference";
  auxI->label        = "<3>";
  auxI->validator    = "<3>";
  //auxI->renderer     = "contplot";
  //auxI->params       = "--projectdir <projectdir> --p <paramsfile> --spectrumscan <0> --<12> <1> --peptide <4> --target cout";
  auxB->sequences.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel  = "Homolog";
  auxI->label        = "<4>";
  auxI->validator    = "<4>";
  //auxI->renderer     = "contplot";
  //auxI->params       = "--projectdir <projectdir> --p <paramsfile> --spectrumscan <0> --<12> <1> --peptide <5> --target cout";
  auxB->sequences.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel  = "de Novo";
  auxI->label        = "<5>";
  auxI->validator    = "<5>";
  //auxI->renderer     = "contplot";
  //auxI->params       = "--projectdir <projectdir> --p <paramsfile> --spectrumscan <0> --<12> <1> --peptide <6> --target cout";
  auxB->sequences.push_back(auxI);

  auxI = new ReportColumnTypeImageOnDemand();
  auxI->dynamic      = true;
  auxI->columnLabel  = "User";
  auxI->label        = "<6>";
  //auxI->renderer     = "contplot";
  //auxI->params       = "--projectdir <projectdir> --p <paramsfile> --spectrumscan <0> --<12> <1> --peptide <6> --target cout";
  auxB->sequences.push_back(auxI);

  auxS = new ReportColumnTypeString();
  auxS->isInput      = true;
  auxS->dynamic      = true;
  auxS->id           = "input_contig_<row>_<col>";
  auxB->sequences.push_back(auxS);

  auxS = new ReportColumnTypeString();
  auxS->isButton     = true;
  auxS->dynamic      = true;
  auxS->text         = "Update";
  auxS->link         = "--update <0> --column 6 --data <7>";
  auxS->onClick      = "javascript:DoOnCick(\"input_contig_<row>_<col>\", n, 6);";
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);


  // colTypes[3] -> (CTstring) protein
  auxS = new ReportColumnTypeString();
  auxS->columnLabel  = "Protein";
  auxS->text         = "<7><br><br><8>";
  auxS->link         = "protein.<1>.0.html";
  m_colTypes.push_back(auxS);
}
////////////////////////////////////////////////////////////////////////////////
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
void ReportTableContig::defineView2(void)
{
  // clear view
  clearView();

  // auxiliary object for column type holding
  ReportColumnTypeString        *auxS;
  ReportColumnTypeImageOnDemand *auxI;
  ReportColumnTypeBox           *auxB;


  auxB = new ReportColumnTypeBox();
  auxB->columnLabel  = "";

  // the image
  auxI = new ReportColumnTypeImageOnDemand();
  auxI->columnLabel       = "";
  //auxI->label             = "<br>Reference<br><3><br>Homolog<br><4><br>de Novo<br><5>";
  auxI->iconRenderer      = "contplot";
  //auxI->iconParams        = "--projectdir <projectdir> --contig --peptide <3> --target internal --encoding uu64";
  //auxI->iconParams         = "contig.<0>.png";

  auxI->iconParams.push_back( ReportParamsOption("--star",            "<16>",    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--abruijn",         "<15>",    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--seqs",            "<17>",    ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--star",            "<projectdir>/<16>",    ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--abruijn",         "<projectdir>/<15>",    ""    ) );
//  auxI->iconParams.push_back( ReportParamsOption("--seqs",            "<projectdir>/<17>",    ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--title",           "Contig <0>",                               ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--contig",          "<0>",                                      "<0>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--reference",       "<3>",                                      "<3>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--homolog",         "<4>",                                      "<4>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--denovo",          "<5>",                                      "<5>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--user",            "<6>",                                      "<6>" ) );
  auxI->iconParams.push_back( ReportParamsOption("--mass-reference",  "<10>",                                     "<10>") );
  auxI->iconParams.push_back( ReportParamsOption("--mass-homolog",    "<11>",                                     "<11>") );
  auxI->iconParams.push_back( ReportParamsOption("--offset-reference","<12>",                                     "<12>") );
  auxI->iconParams.push_back( ReportParamsOption("--offset-homolog",  "<13>",                                     "<13>") );
  auxI->iconParams.push_back( ReportParamsOption("--reverse",         "",                                         "<14>") );
  auxI->iconParams.push_back( ReportParamsOption("--target",          "internal",                                 ""    ) );
  auxI->iconParams.push_back( ReportParamsOption("--encoding",        "uu64",                                     ""    ) );

  auxB->sequences.push_back(auxI);


  // the input sequence
  auxS = new ReportColumnTypeString();
  auxS->isInput      = true;
  auxS->dynamic      = true;
  auxS->id           = "input_contig_<row>_<col>";
  auxB->sequences.push_back(auxS);

  // the "update" button
  auxS = new ReportColumnTypeString();
  auxS->isButton     = true;
  auxS->dynamic      = true;
  auxS->text         = "Update";
  auxS->link         = "--update <0> --pklbin <1> --peptide <5>";
  auxS->onClick      = "javascript:DoOnCick(\"input_contig_<row>_<col>\", n, 5);";
  auxB->sequences.push_back(auxS);

  m_colTypes.push_back(auxB);
}
////////////////////////////////////////////////////////////////////////////////
}; // namespace
////////////////////////////////////////////////////////////////////////////////
