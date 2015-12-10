/*
 * spectrum_statistics.cpp
 *
 *  Created on: Sep 23, 2010
 *      Author: jsnedecor
 */

#include "spectrum.h"
#include "spectrum_annot_statistics.h"
#include "spectrum_annot_statistics_parse.h"
#include "inspect_parse.h"

using namespace std;

void print_help(const char *message)
{
  cout
      << "***************************************************************************\n\n%s\n"
      << message << endl;

  cout << "Required arguments:" << endl;
  cout << "-------------------" << endl;
  cout << endl;

  cout << "--inspect  - inspect results. Note! Annotation and Scan# must be defined!" << endl;
  cout << "--snets    - snets results. Note! Annotation and Scan# must be defined!" << endl;
  cout << endl;
  cout << "--pklbin <path to input file>  - pklbin format file" << endl;
  cout << "--ms2_model <path to input file> - ms2 model format file" << endl;
  cout << "--spectra_stats <path to input file> - spectra_stats file" << endl;
  cout << "--srm_offset default 0, should be -1.007 for stars/prm" << endl;
  cout << "--prm_offset default 0, should be -1.007 for stars/prm" << endl;
  cout << "--peak_tol default 0.5" << endl;
  cout << endl << endl;
}

int main(int argc, char **argv)
{
  int i;
  char inspect_file[256];
  char snets_file[256];
  char pklbin_file[256];
  char model_file[256];
  char spectra_stats_file[256];
  float srm_offset = 0.0;
  float prm_offset = 0.0;
  float peak_tol = 0.5;

  bool got_inspect_file = false;
  bool got_snets_file = false;
  bool got_pklbin_file = false;
  bool got_model = false;
  bool got_spectra_stats = false;

  i = 1;

  while (argc > i)
  {
    if (strcmp(argv[i], "--inspect") == 0)
    {
      if (++i == argc)
      {
        print_help("Missing inspect result file!");
        return 1;
      }

      strcpy(inspect_file, argv[i]);

      got_inspect_file = true;

    }
    else if (strcmp(argv[i], "--snets") == 0)
    {
      if (++i == argc)
      {
        print_help("Missing snets result file!");
        return 1;
      }

      strcpy(snets_file, argv[i]);

      got_snets_file = true;
    }
    else if (strcmp(argv[i], "--pklbin") == 0)
    {
      if (++i == argc)
      {
        print_help("Missing pklbin file!");
        return 1;
      }

      strcpy(pklbin_file, argv[i]);

      got_pklbin_file = true;
    }
    else if (strcmp(argv[i], "--ms2_model") == 0)
    {
      if (++i == argc)
      {
        print_help("Missing ms2model file!");
        return 1;
      }

      strcpy(model_file, argv[i]);
      got_model = true;
    }
    else if (strcmp(argv[i], "--spectra_stats") == 0)
    {
      if (++i == argc)
      {
        print_help("Missing spectra stats file!");
        return 1;
      }

      strcpy(spectra_stats_file, argv[i]);
      got_spectra_stats = true;
    }
    else if (strcmp(argv[i], "--prm_offset") == 0)
    {
      if (++i == argc)
      {
        print_help("Missing prm_offset variable!");
        return 1;
      }
      prm_offset = (float)atof(argv[i]);
    }
    else if (strcmp(argv[i], "--srm_offset") == 0)
    {
      if (++i == argc)
      {
        print_help("Missing srm_offset variable!");
        return 1;
      }
      srm_offset = (float)atof(argv[i]);
    }
    else if (strcmp(argv[i], "--peak_tol") == 0)
    {
      if (++i == argc)
      {
        print_help("Missing peak_tol variable!");
      }
      peak_tol = (float)atof(argv[i]);
    }
    else
    {
      i++;
    }
  }

  if ((got_inspect_file || got_snets_file) && got_pklbin_file && got_model
      && got_spectra_stats)
  {
    if (got_inspect_file)
    {
    cerr << "inspect_file " << inspect_file << " pklbin file " << pklbin_file
        << " ms2_model " << model_file << " spectra_stats_file"
        << spectra_stats_file << endl;
    }
    else
    {
      cerr << "inspect_file " << snets_file << " pklbin file " << pklbin_file
          << " ms2_model " << model_file << " spectra_stats_file"
          << spectra_stats_file << endl;
    }
  }
  else
  {
    print_help("Missing --inspect or --snets --ms2_model, --pklbin or --spectra_stats parameter!");
    return 1;
  }

  if (got_inspect_file && got_snets_file)
  {
    cerr << "Only one results file allowed, either snets or inspect" << endl;
    return 1;
  }

  SpecSet specs;

  if (specs.LoadSpecSet_pklbin(pklbin_file) <= 0)
  {
    cerr << "ERROR loading " << pklbin_file << "!" << endl;
    return 1;
  }

  SpecSet results;

  InspectResultsSet peptide_results;

  if (got_inspect_file)
  {
    //load in Inspect file
    if (!peptide_results.loadInspectResultsFile(inspect_file))
    {
      cerr << "ERROR loading " << inspect_file << "!" << endl;
      return 1;
    }
  }

  if (got_snets_file)
  {
    if (!peptide_results.loadInspectResultsFile(snets_file))
    {
      cerr << "ERROR loading " << snets_file << "!" << endl;
      return 1;
    }
  }

  MS2ScoringModel model;

  if (!model.LoadModel(model_file))
  {
    cerr << "ERROR loading " << model_file << "!" << endl;
    return 1;
  }

  SpectrumAnnotParams spectra_stats;

  if (!spectra_stats.loadSpectrumAnnotFile(spectra_stats_file))
  {
    cerr << "ERROR loading" << spectra_stats_file << "!" << endl;
    return 1;
  }

  cerr << " size of inspect data " << peptide_results.results.size() << endl;

  SpectrumAnnotStatistics stats;

  map<int, InspectResultsLine>::iterator result_iterator;

  //annotate spectra and output stats
  for (result_iterator = peptide_results.results.begin(); result_iterator
      != peptide_results.results.end(); result_iterator++)
  {
    int scan_num = result_iterator->first;

    //set to 1 index
    scan_num = scan_num + 1;

    InspectResultsLine curr_line = result_iterator->second;

    string ion_types = "all";
    Spectrum* curr_spec = specs.getScan(scan_num);
    string specnets;

    if (got_inspect_file) 
    {
      peptide_results.inspectToSpecNets(curr_line.Annotation, specnets);
    }
    else
    {
      specnets = curr_line.Annotation;
    }
      
    //annotate
    curr_spec->annotate(specnets,
                        ion_types,
                        model,
                        prm_offset,
                        srm_offset,
                        peak_tol);

    //output stats
    cout << scan_num;

    for (int i = 0; i < spectra_stats.params.size(); i++)
    {
      SpectrumAnnotParam *curr_param = &(spectra_stats.params[i]);

      cout << ",";

      if (curr_param->statistic.compare("%explained intensity") == 0)
      {
        if (curr_param->ionNames.compare("na") == 0)
        { // we're ignoring this parameter
          //do nothing
        }
        else
        {
          float explained_intensity =
              stats.percentExplainedIntensity(*curr_spec, curr_param->ionNames);
          cout << explained_intensity;
        }
      }
      else if (curr_param->statistic.compare("%explained peaks") == 0)
      {
        float explained_peaks =
            stats.percentExplainedPeaks(*curr_spec, curr_param->ionNames);
        cout << explained_peaks;
      }
      else if (curr_param->statistic.compare("%observed ions") == 0)
      {
        if (curr_param->ionNames.compare("na") == 0)
        { // check if we're ignoring this parameter
          // do nothing
        }
        else
        {
          float observed_ions = stats.percentObservedIons(*curr_spec,
                                                          curr_param->ionNames);
          cout << observed_ions;
        }
      }
      else if (curr_param->statistic.compare("total peaks") == 0)
      {
        float total_peaks = stats.totalPeaks(*curr_spec);
        cout << total_peaks;
      }
      else if (curr_param->statistic.compare("parent mass error ppm") == 0)
      {
        float parent_mass_error_ppm =
            stats.parentMassErrorPPM(*curr_spec, curr_line.Charge);
        cout << parent_mass_error_ppm;
      }
      else if (curr_param->statistic.compare("parent mass error da") == 0)
      {
        float parent_mass_error_da = stats.parentMassErrorDa(*curr_spec,
                                                             curr_line.Charge);
        cout << parent_mass_error_da;
      }
      else if (curr_param->statistic.compare("%observed breaks") == 0)
      {
        float observed_breaks = stats.observedBreaks(*curr_spec,
                                                     curr_param->ionNames);
        cout << observed_breaks;
      }
      else
      {
        cerr << "Unknown parameter type! " << endl;
      }
    }
    cout << endl;
  }
}
