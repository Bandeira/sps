////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "ExecReportSpsplot.h"
#include "Logger.h"

#include "aminoacid.h"
#include "utils.h"

#include "ReportInterface.h"


using namespace specnets;
using namespace std;
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
  ExecReportSpsplot::ExecReportSpsplot(void)
  {
    m_name = "ExecReportSpsplot";
    m_type = "ExecReportSpsplot";
  }

  // -------------------------------------------------------------------------
  ExecReportSpsplot::ExecReportSpsplot(const ParameterList & inputParams) :
    ExecBase(inputParams)
  {
    m_name = "ExecReportSpsplot";
    m_type = "ExecReportSpsplot";
  }
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportSpsplot::validateParams(std::string & error)
{
  m_isValid = false;

  //VALIDATE_PARAM_EXIST("OUTDIR");
  VALIDATE_PARAM_EXIST("FILE_STARS");
  VALIDATE_PARAM_EXIST("FILE_COMP");
  VALIDATE_PARAM_EXIST("FILE_SEQS");
  VALIDATE_PARAM_EXIST("FILE_MS");
  VALIDATE_PARAM_EXIST("FILE_MIDX");
  VALIDATE_PARAM_EXIST("FILE_FASTA");
  VALIDATE_PARAM_EXIST("FILE_REFINDEX");
  VALIDATE_PARAM_EXIST("FILE_REFMP");

  m_isValid = true;
  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportSpsplot::merge(void)
{
  return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
vector<ExecBase *> const & ExecReportSpsplot::split(int numSplit)
{
	m_subModules.resize(0);
	return m_subModules;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
ExecBase * ExecReportSpsplot::clone(const ParameterList & input_params) const
{
  return new ExecReportSpsplot(input_params);
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportSpsplot::loadOutputData(void)
{
  return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportSpsplot::saveInputData(std::vector<std::string> & filenames)
{
  return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportSpsplot::saveOutputData(void)
{
  return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
bool ExecReportSpsplot::loadInputData(void)
{
	return false;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////
#define SET_PARAM(var,name) if(m_params.exists(name)){var=name;var+='=';var+=m_params.getValue(name);var+="\n";}


bool ExecReportSpsplot::invoke(void)
{
  //////////////////////////////////////////////////////////////////////////////
  // Execute spsReports

  ParameterList spsReportsParams;
  spsReportsParams.addIfExists(m_params,       "PROJECT_DIR");
  spsReportsParams.addIfExists(m_params,       "EXE_DIR");

  if(m_params.exists("REPORT_JOB"))
    spsReportsParams.setValue("JOB", m_params.getValue("REPORT_JOB"));
    
  if(m_params.exists("REPORT_USER"))
    spsReportsParams.setValue("USER", m_params.getValue("REPORT_USER"));

  if(m_params.exists("REPORT_SERVER"))
    spsReportsParams.setValue("SERVER", m_params.getValue("REPORT_SERVER"));

  if(m_params.exists("REPORT_DIR_SERVER"))
    spsReportsParams.setValue("PROJECT_DIR_SERVER", m_params.getValue("REPORT_DIR_SERVER"));

  if(m_params.exists("REPORT_CELLS_PER_LINE"))
    spsReportsParams.setValue("CELLS_PER_LINE", m_params.getValue("REPORT_CELLS_PER_LINE"));

  spsReportsParams.addIfExists(m_params,       "HTML_DEFS");
  spsReportsParams.addIfExists(m_params,       "FILE_REFMIDX");
  spsReportsParams.addIfExists(m_params,       "FILE_REFMP");
  spsReportsParams.addIfExists(m_params,       "FILE_REFINDEX");
  spsReportsParams.addIfExists(m_params,       "FILE_FASTA");
  spsReportsParams.addIfExists(m_params,       "FILE_MIDX");
  spsReportsParams.addIfExists(m_params,       "FILE_MP");
  spsReportsParams.addIfExists(m_params,       "FILE_INDEX");
  spsReportsParams.addIfExists(m_params,       "FILE_MS");
  spsReportsParams.addIfExists(m_params,       "FILE_SEQS");
  spsReportsParams.addIfExists(m_params,       "FILE_COMP");
  spsReportsParams.addIfExists(m_params,       "FILE_STARS");
  spsReportsParams.addIfExists(m_params,       "FILE_CLUSTERMS");
  spsReportsParams.addIfExists(m_params,       "FILE_CLUSTERSCAN");
  spsReportsParams.addIfExists(m_params,       "MATCHED_CONTIGS");
  spsReportsParams.addIfExists(m_params,       "MATCHED_CONTIGS_IDX");

  spsReportsParams.addIfExists(m_params,       "FILE_MP2");
  spsReportsParams.addIfExists(m_params,       "FILE_MIDX2");

  spsReportsParams.addIfExists(m_params,       "FONT_PATH");

  spsReportsParams.setValue("buildTables",     ""       );

  if(m_params.exists("REPORT_MSMS_IMAGES")) {
    int aux = m_params.getValueInt("REPORT_MSMS_IMAGES");
    if(aux == 0)
      spsReportsParams.setValue("NO_MSMS_IMAGES", "");
  }

  if(m_params.exists("REPORT_DYNAMIC")) {
    //cout << "Building dynamic HTML reports" << endl;
    spsReportsParams.setValue("REPORT_DYNAMIC",""       );
  } else {
    //cout << "Building static HTML reports" << endl;
    spsReportsParams.setValue("buildhtml",       ""       );
  }

  if(m_params.exists("OUTDIR")) {

    DEBUG_MSG("Invoking spsReports");

    string aux = m_params.getValue("OUTDIR");
    spsReportsParams.setValue("REPORT_DIR", aux);
    //spsReportsParams.setValue("TABLES_DIR", aux);

    spsReports::ReportInterface  rep;
    // run
    int ret = rep.processOptions(spsReportsParams);
  }

  // Exit program
  return true;
}
////////////////////////////////////////////////////////////////////////////////
/*----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////////////

