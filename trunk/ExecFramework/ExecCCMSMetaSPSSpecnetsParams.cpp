//Module Includes
#include "SpectralLibrary.h"
#include "utils.h"
#include "projectionutils.h"

// Header Include
#include "ExecCCMSMetaSPSSpecnetsParams.h"

// System Include
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

using namespace specnets;
using namespace std;

namespace specnets
{

  ExecCCMSMetaSPSSpecnetsParams::ExecCCMSMetaSPSSpecnetsParams(void)
  {
    m_name = "ExecCCMSMetaSPSSpecnetsParams";
    m_type = "ExecCCMSMetaSPSSpecnetsParams";
  }

  ExecCCMSMetaSPSSpecnetsParams::ExecCCMSMetaSPSSpecnetsParams(const ParameterList & inputParams)
  {
    m_name = "ExecCCMSMetaSPSSpecnetsParams";
    m_type = "ExecCCMSMetaSPSSpecnetsParams";
    //inputParams.print(cout);
    m_params = inputParams;

    // pk tol
    // pm tol Da
    // ppm tol Da
    // instrument type
    // merge same prec
    // fasta database
  }

  ExecCCMSMetaSPSSpecnetsParams::~ExecCCMSMetaSPSSpecnetsParams(void)
  {
  }

  ExecBase * ExecCCMSMetaSPSSpecnetsParams::clone(const ParameterList & inputParams) const
  {
    return new ExecCCMSMetaSPSSpecnetsParams(inputParams);
  }

  bool ExecCCMSMetaSPSSpecnetsParams::invoke(void)
  {
    if (!m_params.exists("Default.instrument"))
    {
      ERROR_MSG("Parameter \'Default.instrument\' not specified!!");
      return false;
    }
    string IT_ID = "ESI-ION-TRAP";
    string FT_ID = "FT-HYBRID";
    string Inst_ID = m_params.getValue("Default.instrument");
    if (Inst_ID == IT_ID)
    {
      m_params.setValue("DECONV_MS2", "0");
      m_params.setValue("INSTRUMENT_TYPE", "IT");
    }
    else if (Inst_ID == FT_ID)
    {
      m_params.setValue("DECONV_MS2", "1");
      m_params.setValue("INSTRUMENT_TYPE", "FT");
    }
    else
    {
      ERROR_MSG("Unrecognized value \'" << Inst_ID << " for parameter \'Default.instrument\'!!");
      return false;
    }

    if (!m_params.exists("Default.fragmentation"))
    {
      ERROR_MSG("Parameter \'Default.fragmentation\' not specified!!");
      return false;
    }

    m_params.setValue("CLUSTER_MIN_SIZE", "1");

    string Mult_ID = "Multiple";
    string Merge1_ID = "CID/ETD";
    string Merge2_ID = "HCD/ETD";
    string Merge3_ID = "CID/HCD/ETD";
    string CID_ID = "CID";
    string HCD_ID = "HCD";
    string ETD_ID = "ETD";
    string Frag_ID = m_params.getValue("Default.fragmentation");
    if (Frag_ID == Mult_ID)
    {
      m_params.setValue("CLUSTER_TOOL", "PrmClust");
      m_params.setValue("NUM_CONSECUTIVE", "1");
      m_params.setValue("MERGE_SAME_PREC", "0");
      m_params.setValue("MAX_PVALUE", "0.045");
      m_params.setValue("MIN_METACONTIG_SCORE", "3.0");
    }
    else if (Frag_ID == CID_ID)
    {
      m_params.setValue("ACTIVATION", "CID");
      m_params.setValue("CLUSTER_TOOL", "MSCluster");
      m_params.setValue("NUM_CONSECUTIVE", "1");
      m_params.setValue("MERGE_SAME_PREC", "0");
      m_params.setValue("MAX_PVALUE", "0.045");
      m_params.setValue("MIN_METACONTIG_SCORE", "3.0");
    }
    else if (Frag_ID == HCD_ID)
    {
      m_params.setValue("ACTIVATION", "HCD");
      m_params.setValue("CLUSTER_TOOL", "PrmClust");
      m_params.setValue("NUM_CONSECUTIVE", "1");
      m_params.setValue("MERGE_SAME_PREC", "0");
      m_params.setValue("MAX_PVALUE", "0.045");
      m_params.setValue("MIN_METACONTIG_SCORE", "3.0");
    }
    else if (Frag_ID == ETD_ID)
    {
      m_params.setValue("ACTIVATION", "ETD");
      m_params.setValue("CLUSTER_TOOL", "PrmClust");
      m_params.setValue("NUM_CONSECUTIVE", "1");
      m_params.setValue("MERGE_SAME_PREC", "0");
      m_params.setValue("MAX_PVALUE", "0.045");
      m_params.setValue("MIN_METACONTIG_SCORE", "3.0");
    }
    else if (Frag_ID == Merge1_ID || Frag_ID == Merge2_ID)
    {
      m_params.setValue("CLUSTER_TOOL", "PrmClust");
      m_params.setValue("NUM_CONSECUTIVE", "2");
      m_params.setValue("MERGE_SAME_PREC", "1");
      m_params.setValue("MAX_PVALUE", "0.05");
      m_params.setValue("MIN_METACONTIG_SCORE", "3.3");
    }
    else if (Frag_ID == Merge3_ID)
    {
      m_params.setValue("CLUSTER_TOOL", "PrmClust");
      m_params.setValue("NUM_CONSECUTIVE", "3");
      m_params.setValue("MERGE_SAME_PREC", "1");
      m_params.setValue("MAX_PVALUE", "0.05");
      m_params.setValue("MIN_METACONTIG_SCORE", "3.3");
    }
    else
    {
      ERROR_MSG("Unrecognized value \'" << Frag_ID << " for parameter \'Default.fragmentation\'!!");
      return false;
    }

    if (!m_params.exists("tolerance.PM_tolerance"))
    {
      ERROR_MSG("Parameter \'tolerance.PM_tolerance\' not specified!!");
      return false;
    }
    if (!m_params.exists("tolerance.Ion_tolerance"))
    {
      ERROR_MSG("Parameter \'tolerance.Ion_tolerance\' not specified!!");
      return false;
    }
    if (!m_params.exists("tolerance.ppm_PM_tolerance"))
    {
      ERROR_MSG("Parameter \'tolerance.ppm_PM_tolerance\' not specified!!");
      return false;
    }
    if (!m_params.exists("tolerance.ppm_Ion_tolerance"))
    {
      ERROR_MSG("Parameter \'tolerance.ppm_Ion_tolerance\' not specified!!");
      return false;
    }

    m_params.setValue("TOLERANCE_PM_PPM",
                      m_params.getValue("tolerance.ppm_PM_tolerance"));
    m_params.setValue("TOLERANCE_PEAK_PPM",
                      m_params.getValue("tolerance.ppm_Ion_tolerance"));
    m_params.setValue("TOLERANCE_PM",
                      m_params.getValue("tolerance.PM_tolerance"));
    m_params.setValue("TOLERANCE_PEAK",
                      m_params.getValue("tolerance.Ion_tolerance"));

    if (!m_params.exists("Default.complexity"))
    {
      ERROR_MSG("Parameter \'Default.complexity\' not specified!!");
      return false;
    }
    string Simple_ID = "purified simple mixture";
    string Complex_ID = "complex mixture";
    string complexity_ID = m_params.getValue("Default.complexity");
    if (complexity_ID == Simple_ID)
    {
      m_params.setValue("SPS_MIN_EDGES_TO_COMPONENT", "1");
    }
    else if (complexity_ID == Complex_ID)
    {
      m_params.setValue("SPS_MIN_EDGES_TO_COMPONENT", "3");
    }
    else
    {
      ERROR_MSG("Unrecognized value \'" << complexity_ID << " for parameter \'Default.complexity\'!!");
      return false;
    }

    m_params.setValue("PEPNOVO_PTMS", "D+1:N+1:Q-17");
    m_params.setValue("REPORT_SERVER", "http://rodney.ucsd.edu/cgi-bin/");
    m_params.setValue("REPORT_USER", "aguthals");
    m_params.setValue("REPORT_DYNAMIC", "1");
    m_params.setValue("REPORT_DIR", "report");
    m_params.setValue("GRID_NUMNODES", "-1");
    m_params.setValue("GRID_NUMCPUS", "1");
    m_params.setValue("GRID_SGE_EXE_DIR", "/opt/ge2011.11/bin/linux-x64");
    m_params.setValue("MIN_OVERLAP_AREA", "0.450000");
    m_params.setValue("FILTER_TRIGS", "no");
    m_params.setValue("MAX_MOD_MASS", "100");
    m_params.setValue("MIN_MOD_MASS", "-100");
    m_params.setValue("MIN_RATIO", "0.35");
    m_params.setValue("MIN_SPECTRUM_QUALITY", "0");
    m_params.setValue("FIX_CHARGE_ZEROS", "1");
    m_params.setValue("CORRECT_PM", "no");
    m_params.setValue("GUESS_CHARGE", "no");
    m_params.setValue("MIN_MATCHED_PEAKS", "6");
    m_params.setValue("RESOLUTION", "0.01");
    m_params.setValue("TAG_LEN", "6");
    m_params.setValue("MAX_NUM_TAGS", "400");
    m_params.setValue("MAX_NUM_MODS", "2");
    m_params.setValue("MIN_MATCHED_PEAKS_DB", "7");
    m_params.setValue("MAX_PARSIMONY", "1");
    m_params.setValue("CLUSTALW_MINSCORE", "500");
    m_params.setValue("MIN_METACONTIG_SIZE", "1");
    m_params.setValue("SPSPATH_MIN_NUM_PEAKS", "5");
    m_params.setValue("SPSPATH_MIN_NUM_SPECS", "2");
    m_params.setValue("PENALTY_PTM", "-2000");
    m_params.setValue("PENALTY_SAME_VERTEX=", "-1000000");
    m_params.setValue("SPEC_TYPE_MSMS=", "0");

    //Reading File Names
    int mapping_count = 0;
    map<string, string> file_name_mapping;
    std::vector<std::string> search_spectra_names;
    while (1)
    {
      char buf[100];
      sprintf(buf, "upload_file_mapping%i", mapping_count);

      mapping_count++;

      if (!m_params.exists(buf))
        break;

      std::string mapping = m_params.getValue(buf);
      std::string mangled_name = mapping.substr(0, mapping.find("|"));
      std::string original_name = mapping.substr(mapping.find("|") + 1);

      if (mangled_name.find("spec-") == string::npos)
        continue;

      file_name_mapping[original_name] = mangled_name;

      std::string path_to_spectra = m_params.getValue("SPECTRA_DIR", "");
      search_spectra_names.push_back(path_to_spectra + "/" + mangled_name);
    }

    //Constructing input spectra
    string all_input_spectra_specnetsformat = "";
    for (int i = 0; i < search_spectra_names.size(); i++)
    {
      all_input_spectra_specnetsformat += search_spectra_names[i];
      all_input_spectra_specnetsformat += ";";
    }

    m_params.setValue("INPUT_SPECS_MS",
                      all_input_spectra_specnetsformat.c_str());

    char buf[1024];
    getcwd(buf, 1024);
    cout << "CWD" << "\t" << buf << endl;
    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::loadInputData(void)
  {

    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::saveOutputData(void)
  {
    string param_outfile = m_params.getValue("RESULTS_DIR");

    DEBUG_MSG("EXE_DIR"<<"\t"<<m_params.getValue("EXE_DIR"));

    DEBUG_MSG("SAVING\t"<<param_outfile);
    m_params.writeToFile(param_outfile);
    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::saveInputData(std::vector<std::string> & filenames)
  {

    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::loadOutputData(void)
  {

    return true;
  }

  std::vector<ExecBase *> const & ExecCCMSMetaSPSSpecnetsParams::split(int numSplit)
  {

    //m_subModules.push_back
    return m_subModules;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::merge(void)
  {

    return true;
  }

  bool ExecCCMSMetaSPSSpecnetsParams::validateParams(std::string & error)
  {
    return true;
  }

}
