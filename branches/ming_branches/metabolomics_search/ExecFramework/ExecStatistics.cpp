// Header Include
#include "ExecStatistics.h"

// Module Includes
#include "Logger.h"
#include "FileUtils.h"
#include "utils.h"
#include "ExecMergeConvert.h"

// System Includes
#include <stdio.h>

using namespace std;

namespace specnets
{
  // -------------------------------------------------------------------------
  ExecStatistics::ExecStatistics(void) :
    m_spectra(0x0), m_model(0x0), m_statsParams(0x0), m_peptideResults(0x0),
        m_spectraStats(0x0), m_spectraStatsHeader(0x0), ownInput(true),
        ownOutput(true)
  {
    m_name = "ExecStatistics";
    m_type = "ExecStatistics";
  }

  // -------------------------------------------------------------------------
  ExecStatistics::ExecStatistics(const ParameterList & inputParams) :
    ExecBase(inputParams), m_spectra(0x0), m_model(0x0), m_statsParams(0x0),
        m_peptideResults(0x0), m_spectraStats(0x0), m_spectraStatsHeader(0x0),
        ownInput(true), ownOutput(true)
  {
    m_name = "ExecStatistics";
    m_type = "ExecStatistics";
  }

  // -------------------------------------------------------------------------
  ExecStatistics::ExecStatistics(const ParameterList & inputParams,
                                 SpecSet * spectra,
                                 MS2ScoringModel * model,
                                 SpectrumAnnotParameterList * statsParams,
                                 PeptideSpectrumMatchSet * peptideResults,
                                 vector<vector<float> > * spectraStats,
                                 vector<string> * spectraHeader) :
    ExecBase(inputParams), m_spectra(spectra), m_model(model),
        m_statsParams(statsParams), m_peptideResults(peptideResults),
        m_spectraStats(spectraStats), m_spectraStatsHeader(spectraHeader),
        ownInput(false), ownOutput(false)
  {
    m_name = "ExecStatistics";
    m_type = "ExecStatistics";
  }

  // -------------------------------------------------------------------------
  ExecStatistics::~ExecStatistics(void)
  {
    if (ownInput)
    {
      if (m_spectra)
      {
        delete m_spectra;
        m_spectra = 0x0;
      }
      if (m_model)
      {
        delete m_model;
        m_model = 0x0;
      }
      if (m_statsParams)
      {
        delete m_statsParams;
        m_statsParams = 0x0;
      }
      if (m_peptideResults)
      {
        delete m_peptideResults;
        m_peptideResults = 0x0;
      }
    }

    if (ownOutput)
    {
      if (m_spectraStats)
      {
        delete m_spectraStats;
        m_spectraStats = 0x0;
      }
      if (m_spectraStatsHeader)
      {
        delete m_spectraStatsHeader;
        m_spectraStatsHeader = 0x0;
      }
    }
  }

  // -------------------------------------------------------------------------
  ExecBase * ExecStatistics::clone(const ParameterList & inputParams) const
  {
    return new ExecStatistics(inputParams);
  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::invoke(void)
  {
    if (m_spectraStats == 0x0)
    {
      ownOutput = true;
      m_spectraStats = new vector<vector<float> > ();
      m_spectraStatsHeader = new vector<string> ();
    }

    //Make sure required parameters are set:
    if (!m_spectra or m_spectra->size() == 0)
    {
      DEBUG_MSG("Missing input spectra!");
      return false;
    }

    if (!m_peptideResults or !m_statsParams)
    {
      DEBUG_MSG("Missing stats configuration or annotations!")
      return false;
    }

    if (m_statsParams->m_params.size() == 0)
    {
      DEBUG_MSG("No defined stats parameters!");
      return false;
    }


    float prmOffset = m_params.getValueFloat("PRM_OFFSET");
    float srmOffset = m_params.getValueFloat("SRM_OFFSET");
    float peakTol = 0.05;
    if (m_params.exists("TOLERANCE_PEAK")) {
		peakTol = m_params.getValueFloat("TOLERANCE_PEAK");
	}

    vector<psmPtr>::iterator resultIterator;

    SpectrumAnnotStatistics stats;

    AAJumps jumps(1);
    string ionTypes = "all";

    for (resultIterator = m_peptideResults->m_psmSet.begin(); resultIterator
			!= m_peptideResults->m_psmSet.end(); resultIterator++) {
		PeptideSpectrumMatch* psm = resultIterator->get();
		//annotate
		psm->annotate(psm->m_annotation, ionTypes, *m_model, prmOffset,
				srmOffset, jumps);
	}

    bool addedHeaders = false;
    if (m_params.getValueInt("REPORT_CUMMULATIVE", 0) > 0) {
		vector<float> currStatistics(m_statsParams->m_params.size() + 1);
		currStatistics[0] = -1.0;
		addedHeaders = true;
		m_spectraStatsHeader->push_back("Scan#");
		for (int i = 0; i < m_statsParams->m_params.size(); i++) {
			SpectrumAnnotParameter * currParam = &(m_statsParams->m_params[i]);

			//add to header
			m_spectraStatsHeader->push_back(currParam->statisticName);

			if (currParam->statistic.compare("%explained intensity") == 0) {
				if (currParam->ionNames.compare("na") == 0) { // we're ignoring this parameter
					//do nothing
				} else {
					float explainedIntensity = stats.percentExplainedIntensity(
							*m_peptideResults, currParam->ionNames);
					currStatistics[i + 1] = explainedIntensity;
				}
			} else if (currParam->statistic.compare("%explained peaks") == 0) {
				float explainedPeaks = stats.percentExplainedPeaks(*m_peptideResults,
						currParam->ionNames);
				currStatistics[i + 1] = explainedPeaks;
			} else if (currParam->statistic.compare("%observed ions") == 0) {
				if (currParam->ionNames.compare("na") == 0) { // check if we're ignoring this parameter
					// do nothing
				} else {
					float observedIons = stats.percentObservedIons(*m_peptideResults,
							currParam->ionNames);
					currStatistics[i + 1] = observedIons;
				}
			} else if (currParam->statistic.compare("%observed breaks") == 0) {
				float observedBreaks = stats.observedBreaks(*m_peptideResults,
						currParam->ionNames);
				currStatistics[i + 1] = observedBreaks;
			} else if (currParam->statistic.compare("") != 0) {
				currStatistics[i + 1] = -1.0;
			}
		}
		m_spectraStats->push_back(currStatistics);
	}
    //annotate spectra and output stats
    for (resultIterator = m_peptideResults->m_psmSet.begin(); resultIterator
        != m_peptideResults->m_psmSet.end(); resultIterator++)
    {
      PeptideSpectrumMatch& psm = *(resultIterator->get());

      int scanNum = psm.m_scanNum;

      Spectrum * currSpec = psm.m_spectrum;

      if (currSpec == NULL)
      {
        WARN_MSG("Spectrum for scan " << psm.m_scanNum << " not defined!");
        continue;
      }

      vector<float> currStatistics(m_statsParams->m_params.size() + 1); //vector to hold stats for current spectra

      currStatistics[0] = (float)scanNum;

      //add to header
      if (!addedHeaders && resultIterator == m_peptideResults->m_psmSet.begin())
      {
        m_spectraStatsHeader->push_back("Scan#");
      }

      for (int i = 0; i < m_statsParams->m_params.size(); i++)
      {
        SpectrumAnnotParameter * currParam = &(m_statsParams->m_params[i]);

        //add to header
        if (!addedHeaders && resultIterator == m_peptideResults->m_psmSet.begin())
        {
          m_spectraStatsHeader->push_back(currParam->statisticName);
        }

        if (currParam->statistic.compare("%explained intensity") == 0)
        {
          if (currParam->ionNames.compare("na") == 0)
          { // we're ignoring this parameter
            //do nothing
          }
          else
          {
            float explainedIntensity =
                stats.percentExplainedIntensity(psm, currParam->ionNames);
            currStatistics[i + 1] = explainedIntensity;
          }
        }
        else if (currParam->statistic.compare("%explained peaks") == 0)
        {
          float explainedPeaks =
              stats.percentExplainedPeaks(psm, currParam->ionNames);
          currStatistics[i + 1] = explainedPeaks;
        }
        else if (currParam->statistic.compare("%observed ions") == 0)
        {
          if (currParam->ionNames.compare("na") == 0)
          { // check if we're ignoring this parameter
            // do nothing
          }
          else
          {
            float observedIons = stats.percentObservedIons(psm,
                                                           currParam->ionNames);
            currStatistics[i + 1] = observedIons;
          }
        }
        else if (currParam->statistic.compare("total peaks") == 0)
        {
          int totalPeaks = stats.totalPeaks(psm);
          currStatistics[i + 1] = (float)totalPeaks;
        }
        else if (currParam->statistic.compare("parent mass error ppm") == 0)
        {
          float parentMassErrorPpm =
              stats.parentMassErrorPPM(psm, psm.m_charge);
          parentMassErrorPpm = abs(parentMassErrorPpm / (1 - peakTol)); //output error adjusted by tolerance
          currStatistics[i + 1] = parentMassErrorPpm;
        }
        else if (currParam->statistic.compare("parent mass error da") == 0)
        {
          float parentMassErrorDa =
              stats.parentMassErrorDa(psm, psm.m_charge);
          parentMassErrorDa = abs(parentMassErrorDa / (1 - peakTol)); //output error adjusted by tolerance
          currStatistics[i + 1] = parentMassErrorDa;
        }
        else if (currParam->statistic.compare("%observed breaks") == 0)
        {
          float observedBreaks = stats.observedBreaks(psm,
                                                      currParam->ionNames);
          currStatistics[i + 1] = observedBreaks;
        }
        else if (currParam->statistic.compare("%observed difference") == 0)
        {
          vector<string> ionList;
          stringSplit(currParam->ionNames, ionList, ",");
          float firstValue = stats.percentObservedIons(psm, ionList[0]);
          float secondValue = stats.percentObservedIons(psm, ionList[1]);
          float observedDifference = abs(firstValue - secondValue);
          currStatistics[i + 1] = observedDifference;
        }
        else if (currParam->statistic.compare("%explained difference") == 0)
        {
          vector<string> ionList;
          stringSplit(currParam->ionNames, ionList, ",");
          float firstValue = stats.percentExplainedIntensity(psm,
                                                             ionList[0]);
          float secondValue = stats.percentExplainedIntensity(psm,
                                                              ionList[1]);
          float explainedDifference = abs(firstValue - secondValue);
          currStatistics[i + 1] = explainedDifference;
        }
        else if (currParam->statistic.compare("") != 0) {
        	DEBUG_MSG("Unknown parameter type: <" << currParam->statistic << ">");
        }
      }
      addedHeaders = true;
      m_spectraStats->push_back(currStatistics);
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::loadInputData(void)
  {
    //Load in spectra if they haven't been passed in by
    //another function
    if (m_spectra == 0x0)
    {
      ownInput = true;
      m_spectra = new SpecSet;
      m_statsParams = new SpectrumAnnotParameterList;
      m_peptideResults = new PeptideSpectrumMatchSet;
      m_model = new MS2ScoringModel;
    }

    if (m_params.exists("INPUT_SPECS"))
    {
      if (!m_spectra->LoadSpecSet_pkl(m_params.getValue("INPUT_SPECS").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("INPUT_SPECS"));
        return false;
      }
    }
    else if (m_params.exists("INPUT_SPECS_PKLBIN"))
    {
      if (m_spectra->LoadSpecSet_pklbin(m_params.getValue("INPUT_SPECS_PKLBIN").c_str())
          <= 0)
      {
        DEBUG_MSG("Could not load " << m_params.getValue("INPUT_SPECS_PKLBIN"));
        return false;
      }
    } else if (m_params.exists("INPUT_SPECTRA")) {
		if (!ExecMergeConvert::loadSpecset(m_params.getValue("INPUT_SPECTRA"),
				m_spectra)) {
			return false;
		}
	}

    if (m_spectra->size() == 0)
    {
      DEBUG_MSG("Input spectra size is 0!");
      return false;
    }

    DEBUG_VAR(m_spectra->size());

    if (m_params.getValueInt("SET_SCAN_NUMS", 0) > 0) {
		DEBUG_MSG("Setting the scan number of each spectrum to its one-based index");
		for (int i = 0; i < m_spectra->size(); i++) {
			(*m_spectra)[i].scan = (unsigned int) (i + 1);
		}
	}

    if (m_params.exists("TOLERANCE_PEAK")) {
		float peakTol = m_params.getValueFloat("TOLERANCE_PEAK");
		m_spectra->setPeakTolerance(peakTol, false);
	} else if (m_params.exists("TOLERANCE_PEAK_PPM")) {
		float ppmTol = m_params.getValueFloat("TOLERANCE_PEAK_PPM");
		m_spectra->setPeakTolerance(ppmTol, true);
	}

    //Load in statistics generation parameters
    if (m_params.exists("STATISTICS_CONFIG"))
    {
      if (!m_statsParams->loadSpectrumAnnotFile(m_params.getValue("STATISTICS_CONFIG").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("STATISTICS_CONFIG"));
        return false;
      }
    }

    if (m_statsParams->size() == 0)
    {
      DEBUG_MSG("No statistics asked for!");
      return false;
    }

    if (m_params.exists("MS2_SCORING_MODEL"))
    {
      if (!m_model->LoadModel(m_params.getValue("MS2_SCORING_MODEL").c_str()))
      {
        DEBUG_MSG("Could not load " << m_params.getValue("MS2_SCORING_MODEL"));
        return false;
      }
    }

    if (m_model->probs.size() == 0)
    {
      DEBUG_MSG("No model parameters!");
      return false;
    }

    DEBUG_VAR(m_statsParams->size());

    //Load in inspect / specnets results
    if (m_params.exists("INPUT_RESULTS"))
    {
      DEBUG_MSG("Input results: " << m_params.getValue("INPUT_RESULTS"));
      string resultsType(m_params.getValue("INPUT_RESULTS_TYPE"));
      std::transform(resultsType.begin(), resultsType.end(), resultsType.begin(),
      	  			::tolower);
      if (resultsType.compare("inspect") == 0)
      {
        if (!m_peptideResults->loadInspectResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                      m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      } else if (resultsType.compare("msgfdb") == 0) {
		if (!m_peptideResults->loadMSGFDBResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
				                                     true,
				                                     true,
				                                     m_params.getValueBool("SCAN_ZERO_INDEX"))) {
		  DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
				return false;
		 }
	  } else {
        if (!m_peptideResults->loadSpecnetsResultsFile(m_params.getValue("INPUT_RESULTS").c_str(),
                                                       m_params.getValueBool("SCAN_ZERO_INDEX")))
        {
          DEBUG_MSG("Could not load" << m_params.getValue("INPUT_RESULTS"));
          return false;
        }
      }
    }

    if (m_peptideResults->size() == 0)
    {
      DEBUG_MSG("No peptide annotation results!");
      return false;
    }

    DEBUG_VAR(m_peptideResults->size());

    //associate results with spectra
    m_peptideResults->addSpectra(m_spectra);

    return true;
  }
  // -------------------------------------------------------------------------
  bool ExecStatistics::saveOutputData(void)
  {
    if (m_params.exists("OUTPUT_SPECTRA_STATS_BIN"))
    {
      DEBUG_MSG("Outputting statistics...");
      if (!Save_binArray(m_params.getValue("OUTPUT_SPECTRA_STATS_BIN").c_str(),
                         *m_spectraStats))
      {
        ERROR_MSG("Unable to write bin file!");
        return false;
      }
      DEBUG_MSG("done...");
    }

    if (m_params.exists("OUTPUT_SPECTRA_STATS"))
    {
      ofstream outputFile;

      DEBUG_MSG("Opening output stats file...");
      outputFile.open(m_params.getValue("OUTPUT_SPECTRA_STATS").c_str(),
                      ios::out | ios::trunc);
      if (outputFile.fail())
      {
        ERROR_MSG("Unable to open stats file! " << m_params.getValue("OUTPUT_SPECTRA_STATS"));
        return false;
      }

      if (outputFile.is_open() && outputFile.good())
      {
        //output header

        outputFile << (*m_spectraStatsHeader)[0];
        for (int i = 1; i < m_spectraStatsHeader->size(); i++)
        {
          outputFile << "," << (*m_spectraStatsHeader)[i];
        }
        outputFile << endl;

        for (int i = 0; i < m_spectraStats->size(); i++)
        {
          outputFile << (*m_spectraStats)[i][0]; //scan number

          for (int j = 1; j < (*m_spectraStats)[i].size(); j++)
          {
            outputFile << "," << (*m_spectraStats)[i][j];
          }
          outputFile << endl;
        }
      }
      else
      {
        ERROR_MSG("Unable to open file!");
        return false;
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::saveInputData(std::vector<std::string> & filenames)
  {

  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::loadOutputData(void)
  {
  }

  // -------------------------------------------------------------------------
  vector<ExecBase*> const & ExecStatistics::split(int numSplit)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::merge(void)
  {
    //Mr. Stubbs
  }

  // -------------------------------------------------------------------------
  bool ExecStatistics::validateParams(std::string & error)
  {
    m_isValid = false;

    //VALIDATE_PARAM_EXIST("TOLERANCE_PEAK");

    m_isValid = true;
    return true;
  }
}
