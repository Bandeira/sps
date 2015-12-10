/**
 * Tests for spectrum_annot_statistics class
 *
 */
#include "spectrum_annot_statistics.h"
#include "spectrum_annot_statistics_test.h"

#include <iostream>
#include <fstream>
#include <string>

#define BOOST_TEST_MODULE spectrum_annot spectrum_annot_test
#ifdef NO_DLL
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#else
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#endif

using namespace std;

/*pp QC format:
 * columns:
 * Spectrum index (one-based)
 * % explained intensity
 * % explained peaks
 * % observed b-ions (charge 1 only)
 * % observed y-ions (charge 1 only)
 * % observed neutral-loss ions (-H2O, -NH3, a-ions, b-iso, y-iso, charge 1 only)
 * % observed doubly-charged ions: b, b-iso, b-H2O, b-NH3, y, y-iso, y-H2O, y-NH3
 * % intensity explained by b-ions (charge 1 only)
 * % intensity explained by y-ions (charge 1 only)
 * % intensity explained by neutral-loss ions (-H2O, -NH3, a-ions, b-iso, y-iso, charge 1 only)
 * % intensity explained by doubly-charged ions: b, b-iso, b-H2O, b-NH3, y, y-iso, y-H2O, y-NH3
 * Number of peaks
 * PPM precursor mass error
 * Parent mass error in Da
 */

/*
 * Inspect format:
 * file
 * scan#
 * Annotation
 * a bunch of stuff we don't care about.
 */

bool loadQCFile(const char* file,
                vector<vector<string> >& input_vec,
                const char* delim)
{
  ifstream file_handle(file);

  if (!file_handle.is_open() || !file_handle.good())
  {
    BOOST_FAIL("Unable to open file:" << file);
  }

  //read fields for file
  string line_buff;
  getline(file_handle, line_buff);

  if ('\r' == line_buff[line_buff.size() - 1])
  { //strip off carriage return for dos
    line_buff.resize(line_buff.size() - 1);
  }

  vector<string> curr_line;
  splitText(line_buff.c_str(), curr_line, delim);
  if (curr_line.size() > 1)
    input_vec.push_back(curr_line);

  while (!file_handle.eof())
  {
    getline(file_handle, line_buff);
    vector < string > curr_line;
    splitText(line_buff.c_str(), curr_line, delim);
    if (curr_line.size() > 1)
      input_vec.push_back(curr_line);
  }

  return true;
}

//fixture for pp and inspect information
struct qc_fixture
{
  qc_fixture()
  {
    BOOST_TEST_MESSAGE("Setup fixture.");
    BOOST_TEST_MESSAGE("Load pp ms file:");
    BOOST_REQUIRE(loadQCFile(pp_stats_ms_file, pp_stats_ms, (const char*)","));
    BOOST_TEST_MESSAGE("Load pp prm file:");
    BOOST_REQUIRE(loadQCFile(pp_stats_prm_file, pp_stats_prm, (const char*)","));
    BOOST_TEST_MESSAGE("Load pp stars file:");
    BOOST_REQUIRE(loadQCFile(pp_stats_stars_file,
                             pp_stats_stars,
                             (const char*)","));
    BOOST_TEST_MESSAGE("Load inspect file:");
    BOOST_REQUIRE(loadQCFile(inspect_results_file,
                             inspect_results,
                             (const char*)"\t"));
    BOOST_TEST_MESSAGE("Loading ms:");
    BOOST_REQUIRE(ms_specs.LoadSpecSet_pklbin(specs_ms_file));
    BOOST_TEST_MESSAGE("Loading prm:");
    BOOST_REQUIRE(prm_specs.LoadSpecSet_pklbin(specs_prm_file));
    BOOST_TEST_MESSAGE("Loading stars:");
    BOOST_REQUIRE(stars_specs.LoadSpecSet_pklbin(specs_stars_file));
    BOOST_TEST_MESSAGE("Loading model:");
    BOOST_REQUIRE(model.LoadModel(model_file));

    BOOST_TEST_MESSAGE("Annotating spectra: ");

    //annotate all spectra
    for (int i = 1; i < inspect_results.size(); i++)
    {
      int scan_num = atoi(inspect_results[i][1].c_str()) + 1;
      Spectrum* ms_scan = ms_specs.getScan(scan_num);
      string ion_types = "all";
      ms_scan->annotate(inspect_results[i][2], ion_types, model, 0.0, 0.0, .45);
      Spectrum* prm_scan = prm_specs.getScan(scan_num);
      prm_scan->annotate(inspect_results[i][2],
                         ion_types,
                         model,
                         -1.007,
                         -1.007,
                         .45);
      Spectrum* stars_scan = stars_specs.getScan(scan_num);
      stars_scan->annotate(inspect_results[i][2],
                           ion_types,
                           model,
                           -1.007,
                           -1.007,
                           .45);
    }
  }
  ~qc_fixture()
  {
    BOOST_TEST_MESSAGE("Teardown fixture.");
  }

  vector<vector<string> > pp_stats_ms;
  vector<vector<string> > pp_stats_prm;
  vector<vector<string> > pp_stats_stars;
  vector<vector<string> > inspect_results;
  SpecSet ms_specs;
  SpecSet prm_specs;
  SpecSet stars_specs;
  MS2ScoringModel model;
};

//check to make sure we're correctly catching unannotated spectra
BOOST_AUTO_TEST_CASE(unannotated_spectra)
{
  Spectrum test_spec;
  SpectrumAnnotStatistics stats;
  string ion_types = "all";

  BOOST_CHECK(stats.percentExplainedIntensity(test_spec, ion_types) < 0);
  BOOST_CHECK(stats.percentExplainedPeaks(test_spec, ion_types) < 0);
  BOOST_CHECK(stats.percentObservedIons(test_spec, ion_types) < 0);
  BOOST_CHECK(stats.parentMassErrorPPM(test_spec) < 0);
  BOOST_CHECK(stats.parentMassErrorDa(test_spec) < 0);
  BOOST_CHECK(stats.observedBreaks(test_spec, ion_types) < 0);
}

//next few tests require fixture

BOOST_FIXTURE_TEST_CASE(percentExplainedIntensity, qc_fixture)
{
  SpectrumAnnotStatistics stats;
  string ion_types;
  AAJumps jumps(1);

  BOOST_TEST_MESSAGE("Checking explained intensity ms:");
  for(int i=1; i<pp_stats_ms.size(); i++)
  {
    int scan_num = atoi(pp_stats_ms[i][0].c_str());
    Spectrum* curr_spec = ms_specs.getScan(scan_num);
    ion_types = "all";
    float explained_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_intensity = atof(pp_stats_ms[i][1].c_str());
    BOOST_CHECK_CLOSE(std_explained_intensity, explained_intensity,10);
    ion_types = "b";
    float explained_b_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_b_intensity = atof(pp_stats_ms[i][7].c_str());
    BOOST_CHECK_CLOSE(std_explained_b_intensity, explained_b_intensity,10);
    ion_types = "y";
    float explained_y_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_y_intensity = atof(pp_stats_ms[i][8].c_str());
    BOOST_CHECK_CLOSE(std_explained_y_intensity, explained_y_intensity,10);
    ion_types = "b-iso,b-NH3,b-H2O,b-H2O-H2O,b-H2O-NH3,a,a-H2O,a-NH3,y-iso,y-NH3,y-H2O,y-H2O-NH3,y-H2O-H2O";
    float explained_neutralloss_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_neutralloss_intensity = atof(pp_stats_ms[i][9].c_str());
    BOOST_CHECK_CLOSE(std_explained_neutralloss_intensity, explained_neutralloss_intensity,10);
    ion_types = "b++,b++-NH3,b++-H2O,b++-H2O-H20,b++-H2O-NH3,y++,y++-NH3,y++-H2O,y++-H2O-NH3,y++-H2O-H2O,P++-H2O,P++-NH3,P++";
    float explained_double_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_double_intensity = atof(pp_stats_ms[i][10].c_str());
    BOOST_CHECK_CLOSE(std_explained_double_intensity, explained_double_intensity,10);
  }

  BOOST_TEST_MESSAGE("Checking explained intensity prm: ");
  for(int i=1; i< pp_stats_prm.size(); i++)
  {
    int scan_num = atoi(pp_stats_prm[i][0].c_str());
    Spectrum* curr_spec = prm_specs.getScan(scan_num);
    ion_types = "b,y";
    float explained_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_intensity = atof(pp_stats_prm[i][1].c_str());
    BOOST_CHECK_CLOSE(std_explained_intensity, explained_intensity,10);
    ion_types = "b";
    float explained_b_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_b_intensity = atof(pp_stats_prm[i][7].c_str());
    BOOST_CHECK_CLOSE(std_explained_b_intensity, explained_b_intensity,10);
    ion_types = "y";
    float explained_y_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_y_intensity = atof(pp_stats_prm[i][8].c_str());
    BOOST_CHECK_CLOSE(std_explained_y_intensity, explained_y_intensity,10);
  }

  BOOST_TEST_MESSAGE("Checking explained intensity stars: ");
  for(int i=1; i< pp_stats_stars.size(); i++)
  {
    ion_types = "b,y";
    int scan_num = atoi(pp_stats_stars[i][0].c_str());
    Spectrum* curr_spec = stars_specs.getScan(scan_num);
    float explained_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_intensity = atof(pp_stats_stars[i][1].c_str());
    BOOST_CHECK_CLOSE(std_explained_intensity, explained_intensity,10);
    ion_types = "b";
    float explained_b_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_b_intensity = atof(pp_stats_stars[i][7].c_str());
    BOOST_CHECK_CLOSE(std_explained_b_intensity, explained_b_intensity,10);
    ion_types = "y";
    float explained_y_intensity = stats.percentExplainedIntensity(*curr_spec,ion_types);
    float std_explained_y_intensity = atof(pp_stats_stars[i][8].c_str());
    BOOST_CHECK_CLOSE(std_explained_y_intensity, explained_y_intensity,10);
  }
}

BOOST_FIXTURE_TEST_CASE(percentExplainedPeaks, qc_fixture)
{
  SpectrumAnnotStatistics stats;
  string ion_types;

  BOOST_TEST_MESSAGE("Checking explained peaks ms:");
  for(int i=1; i<pp_stats_ms.size(); i++)
  {
    int scan_num = atoi(pp_stats_ms[i][0].c_str());
    Spectrum* curr_spec = ms_specs.getScan(scan_num);
    ion_types = "all";
    float explained_peaks = stats.percentExplainedPeaks(*curr_spec,ion_types);
    float std_explained_peaks = atof(pp_stats_ms[i][2].c_str());
    BOOST_CHECK_CLOSE(std_explained_peaks, explained_peaks,10);
  }

  BOOST_TEST_MESSAGE("Checking explained peaks prm: ");
  for(int i=1; i< pp_stats_prm.size(); i++)
  {
    int scan_num = atoi(pp_stats_prm[i][0].c_str());
    Spectrum* curr_spec = prm_specs.getScan(scan_num);
    ion_types = "b,y";
    float explained_peaks = stats.percentExplainedPeaks(*curr_spec,ion_types);
    float std_explained_peaks = atof(pp_stats_prm[i][2].c_str());
    BOOST_CHECK_CLOSE(std_explained_peaks, explained_peaks,10);
  }

  BOOST_TEST_MESSAGE("Checking explained peaks stars: ");
  for(int i=1; i< pp_stats_stars.size(); i++)
  {
    ion_types = "b,y";
    int scan_num = atoi(pp_stats_stars[i][0].c_str());
    Spectrum* curr_spec = stars_specs.getScan(scan_num);
    float explained_peaks = stats.percentExplainedPeaks(*curr_spec,ion_types);
    float std_explained_peaks = atof(pp_stats_stars[i][2].c_str());
    BOOST_CHECK_CLOSE(std_explained_peaks, explained_peaks,10);
  }
}

BOOST_FIXTURE_TEST_CASE(percentExplainedIons, qc_fixture)
{
  SpectrumAnnotStatistics stats;
  string ion_types;

  BOOST_TEST_MESSAGE("Checking explained ions ms:");
  for(int i=1; i<pp_stats_ms.size(); i++)
  {
    int scan_num = atoi(pp_stats_ms[i][0].c_str());
    Spectrum* curr_spec = ms_specs.getScan(scan_num);
    ion_types = "b";
    float explained_b_ions = stats.percentObservedIons(*curr_spec,ion_types);
    float std_explained_b_ions = atof(pp_stats_ms[i][3].c_str());
    BOOST_CHECK_CLOSE(std_explained_b_ions, explained_b_ions,10);
    ion_types = "y";
    float explained_y_ions = stats.percentObservedIons(*curr_spec,ion_types);
    float std_explained_y_ions = atof(pp_stats_ms[i][4].c_str());
    BOOST_CHECK_CLOSE(std_explained_y_ions, explained_y_ions,10);
    ion_types = "b-iso,b-NH3,b-H2O,b-H2O-H2O,b-H2O-NH3,a,a-H2O,a-NH3,y-iso,y-NH3,y-H2O,y-H2O-NH3,y-H2O-H2O";
    float explained_neutralloss_ions = stats.percentObservedIons(*curr_spec,ion_types);
    float std_explained_neutralloss_ions = atof(pp_stats_ms[i][5].c_str());
    BOOST_CHECK_CLOSE(std_explained_neutralloss_ions, explained_neutralloss_ions,10);
    ion_types = "b++,b++-NH3,b++-H2O,b++-H2O-H20,b++-H2O-NH3,y++,y++-NH3,y++-H2O,y++-H2O-NH3,y++-H2O-H2O,P++-H2O,P++-NH3,P++";
    float explained_double_ions = stats.percentObservedIons(*curr_spec,ion_types);
    float std_explained_double_ions = atof(pp_stats_ms[i][6].c_str());
    BOOST_CHECK_CLOSE(std_explained_double_ions, explained_double_ions,10);
  }

  BOOST_TEST_MESSAGE("Checking explained ions prm: ");
  for(int i=1; i< pp_stats_prm.size(); i++)
  {
    int scan_num = atoi(pp_stats_prm[i][0].c_str());
    Spectrum* curr_spec = prm_specs.getScan(scan_num);
    ion_types = "b";
    float explained_b_ions = stats.percentObservedIons(*curr_spec,ion_types);
    float std_explained_b_ions = atof(pp_stats_prm[i][3].c_str());
    BOOST_CHECK_CLOSE(std_explained_b_ions, explained_b_ions,10);
    ion_types = "y";
    float explained_y_ions = stats.percentObservedIons(*curr_spec,ion_types);
    float std_explained_y_ions = atof(pp_stats_prm[i][4].c_str());
    BOOST_CHECK_CLOSE(std_explained_y_ions, explained_y_ions,10);
  }

  BOOST_TEST_MESSAGE("Checking explained ions stars: ");
  for(int i=1; i< pp_stats_stars.size(); i++)
  {
    int scan_num = atoi(pp_stats_stars[i][0].c_str());
    Spectrum* curr_spec = stars_specs.getScan(scan_num);
    ion_types = "b";
    float explained_b_ions = stats.percentObservedIons(*curr_spec,ion_types);
    float std_explained_b_ions = atof(pp_stats_stars[i][3].c_str());
    BOOST_CHECK_CLOSE(std_explained_b_ions, explained_b_ions,10);
    ion_types = "y";
    float explained_y_ions = stats.percentObservedIons(*curr_spec,ion_types);
    float std_explained_y_ions = atof(pp_stats_stars[i][4].c_str());
    BOOST_CHECK_CLOSE(std_explained_y_ions, explained_y_ions,10);
  }
}

BOOST_FIXTURE_TEST_CASE(parentMassError, qc_fixture)
{
  SpectrumAnnotStatistics stats;

  BOOST_TEST_MESSAGE("Checking ppm/da error ms:");
  for(int i=1; i<pp_stats_ms.size(); i++)
  {
    int scan_num = atoi(pp_stats_ms[i][0].c_str());
    Spectrum* curr_spec = ms_specs.getScan(scan_num);

    int charge = atoi(inspect_results[i][4].c_str());
    float ppm_error = stats.parentMassErrorPPM(*curr_spec,charge);

    float std_ppm_error = atof(pp_stats_ms[i][12].c_str());
    BOOST_CHECK_MESSAGE(abs(std_ppm_error) + 10 >= abs(ppm_error) && abs(std_ppm_error) - 10 <= abs(ppm_error), "std_ppm_error " << std_ppm_error << " ppm_error " << ppm_error << " scan " << scan_num);

    charge = atoi(inspect_results[i][4].c_str());
    float da_error = stats.parentMassErrorDa(*curr_spec,charge);

    float std_da_error = atof(pp_stats_ms[i][13].c_str());
    BOOST_CHECK_MESSAGE(std_da_error + 1 >= da_error && std_da_error - 1 <= da_error, "std_da_error " << std_da_error << " da_error " << da_error << " scan " << scan_num);
  }

  BOOST_TEST_MESSAGE("Checking ppm/da error prm:");
  for(int i=1; i<pp_stats_prm.size(); i++)
  {
    int scan_num = atoi(pp_stats_prm[i][0].c_str());
    Spectrum* curr_spec = prm_specs.getScan(scan_num);

    int charge = atoi(inspect_results[i][4].c_str());
    float ppm_error = stats.parentMassErrorPPM(*curr_spec,charge);

    float std_ppm_error = atof(pp_stats_prm[i][12].c_str());
    BOOST_CHECK_MESSAGE(abs(std_ppm_error) + 10 >= abs(ppm_error) && abs(std_ppm_error) - 10 <= abs(ppm_error), "std_ppm_error " << std_ppm_error << " ppm_error " << ppm_error << " scan " << scan_num);

    charge = atoi(inspect_results[i][4].c_str());
    float da_error = stats.parentMassErrorDa(*curr_spec,charge);

    float std_da_error = atof(pp_stats_prm[i][13].c_str());
    BOOST_CHECK_MESSAGE(std_da_error + 1 >= da_error && std_da_error - 1 <= da_error, "std_da_error " << std_da_error << " da_error " << da_error << " scan " << scan_num);
  }

  BOOST_TEST_MESSAGE("Checking ppm/da error stars:");
  for(int i=1; i<pp_stats_stars.size(); i++)
  {
    int scan_num = atoi(pp_stats_stars[i][0].c_str());
    Spectrum* curr_spec = stars_specs.getScan(scan_num);

    int charge = atoi(inspect_results[i][4].c_str());
    float ppm_error = stats.parentMassErrorPPM(*curr_spec,charge);

    float std_ppm_error = atof(pp_stats_stars[i][12].c_str());
    BOOST_CHECK_MESSAGE(abs(std_ppm_error) + 10 >= abs(ppm_error) && abs(std_ppm_error) - 10 <= abs(ppm_error), "std_ppm_error " << std_ppm_error << " ppm_error " << ppm_error << " scan " << scan_num);

    charge = atoi(inspect_results[i][4].c_str());
    float da_error = stats.parentMassErrorDa(*curr_spec,charge);

    float std_da_error = atof(pp_stats_stars[i][13].c_str());
    BOOST_CHECK_MESSAGE(std_da_error + 1 >= da_error && std_da_error - 1 <= da_error, "std_da_error " << std_da_error << " da_error " << da_error << " scan " << scan_num);
  }
}

BOOST_FIXTURE_TEST_CASE(totalPeaks, qc_fixture)
{
  SpectrumAnnotStatistics stats;

  BOOST_TEST_MESSAGE("Checking total peaks ms:");
  for(int i=1; i<pp_stats_ms.size(); i++)
  {
    int scan_num = atoi(pp_stats_ms[i][0].c_str());
    Spectrum* curr_spec = ms_specs.getScan(scan_num);
    int total_peaks = stats.totalPeaks(*curr_spec);
    int std_total_peaks = atoi(pp_stats_ms[i][11].c_str());
    BOOST_CHECK_EQUAL(std_total_peaks,total_peaks);
  }

  BOOST_TEST_MESSAGE("Checking total peaks prm:");
  for(int i=1; i<pp_stats_prm.size(); i++)
  {
    int scan_num = atoi(pp_stats_prm[i][0].c_str());
    Spectrum* curr_spec = prm_specs.getScan(scan_num);
    int total_peaks = stats.totalPeaks(*curr_spec);
    int std_total_peaks = atoi(pp_stats_prm[i][11].c_str());
    BOOST_CHECK_EQUAL(std_total_peaks,total_peaks);
  }

  BOOST_TEST_MESSAGE("Checking total peaks stars:");
  for(int i=1; i<pp_stats_stars.size(); i++)
  {
    int scan_num = atoi(pp_stats_stars[i][0].c_str());
    Spectrum* curr_spec = stars_specs.getScan(scan_num);
    int total_peaks = stats.totalPeaks(*curr_spec);
    int std_total_peaks = atoi(pp_stats_stars[i][11].c_str());
    BOOST_CHECK_EQUAL(std_total_peaks,total_peaks);
  }
}

BOOST_FIXTURE_TEST_CASE(observedBreaks, qc_fixture) {
  SpectrumAnnotStatistics stats;
  BOOST_TEST_MESSAGE("Check observed breaks:");
  string ion_types = "b,y,b++,y++";

  BOOST_TEST_MESSAGE("Check observed breaks ms:");
  for(int i=1; i<pp_stats_ms.size();i++)
  {
    int scan_num = atoi(pp_stats_ms[i][0].c_str());
    Spectrum* curr_spec = ms_specs.getScan(scan_num);
    float observed_breaks = stats.observedBreaks(*curr_spec,ion_types);
    float std_observed_breaks = atof(pp_stats_ms[i][14].c_str());
    BOOST_CHECK_CLOSE(std_observed_breaks,observed_breaks,10);
  }

  ion_types = "b,y";
  BOOST_TEST_MESSAGE("Check observed breaks prm:");
  for(int i=1; i<pp_stats_prm.size();i++)
  {
    int scan_num = atoi(pp_stats_prm[i][0].c_str());
    Spectrum* curr_spec = prm_specs.getScan(scan_num);
    float observed_breaks = stats.observedBreaks(*curr_spec,ion_types);
    float std_observed_breaks = atof(pp_stats_prm[i][14].c_str());
    BOOST_CHECK_CLOSE(std_observed_breaks,observed_breaks,10);
  }

  ion_types = "b,y";
  BOOST_TEST_MESSAGE("Check observed breaks stars:");
  for(int i=1; i<pp_stats_stars.size();i++)
  {
    int scan_num = atoi(pp_stats_stars[i][0].c_str());
    Spectrum* curr_spec = stars_specs.getScan(scan_num);
    float observed_breaks = stats.observedBreaks(*curr_spec,ion_types);
    float std_observed_breaks = atof(pp_stats_stars[i][14].c_str());
    BOOST_CHECK_CLOSE(std_observed_breaks,observed_breaks,10);
  }
}
