/**
 * Unit testing for spectrum class.
 */

#include "spectrum.h"
#include "spectrum_test.h"
#define BOOST_TEST_MODULE spectrum spectrum_test
#ifdef NO_DLL
#include <boost/test/included/unit_test.hpp>
#else
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#endif

using namespace std;

void test_peaks(float curr_first_mass, float curr_first_intensity, float curr_last_mass,float curr_last_intensity) {
	BOOST_CHECK_CLOSE(curr_first_mass,first_mass,.3);
	BOOST_CHECK_CLOSE(curr_first_intensity,first_intensity,.3);
	BOOST_CHECK_CLOSE(curr_last_mass,last_mass,.3);
	BOOST_CHECK_CLOSE(curr_last_intensity,last_intensity,.3);
}

void test_spectrum(const Spectrum &test_spec) {
	int curr_num_peaks = test_spec.peakList.size();
	//check to make sure that we've got the right number of peaks and scan #
	BOOST_CHECK_EQUAL(curr_num_peaks,num_peaks);
	BOOST_CHECK_EQUAL(test_spec.scan, scan_num);
	//check peaks for first scan
	test_peaks(test_spec.peakList[0][0], test_spec.peakList[0][1],test_spec.peakList[curr_num_peaks-1][0],test_spec.peakList[curr_num_peaks-1][1]);
}

void test_copy(const Spectrum &spec1, const Spectrum &spec2) {
	BOOST_CHECK_EQUAL(spec1.scan,spec2.scan);
	BOOST_CHECK_EQUAL(spec1.parentMass,spec2.parentMass);
	BOOST_CHECK_EQUAL(spec1.parentCharge, spec2.parentCharge);
    BOOST_CHECK_EQUAL(spec1.resolution,spec2.resolution);
    BOOST_CHECK_EQUAL(spec1.idDist,spec2.idDist);
    BOOST_CHECK_EQUAL(spec1.annotation_peptide,spec2.annotation_peptide);
    for(unsigned int i=0; i< spec1.peakList.size(); i++ ) {
    	BOOST_CHECK_EQUAL(spec1.peakList[i][0],spec2.peakList[i][0]);
    	BOOST_CHECK_EQUAL(spec2.peakList[i][1],spec2.peakList[i][1]);
    }
    for(unsigned int i=0; i< spec1.ionTypes.size(); i++)  {
    	BOOST_CHECK_EQUAL(spec1.ionTypes[i].name,spec2.ionTypes[i].name);
    	BOOST_CHECK_EQUAL(spec1.ionTypes[i].isNTerm,spec2.ionTypes[i].isNTerm);
    	BOOST_CHECK_EQUAL(spec1.ionTypes[i].prob,spec2.ionTypes[i].prob);
    	BOOST_CHECK_EQUAL(spec1.ionTypes[i].charge,spec2.ionTypes[i].charge);
    	BOOST_CHECK_EQUAL(spec1.ionTypes[i].massOffset,spec2.ionTypes[i].massOffset);
    }
    for(unsigned int i=0; i< spec1.annotation.size(); i++) {
    	const ftIonFragment* curr_frag1 = spec1.annotation[i].first;
    	const ftIonFragment* curr_frag2 = spec2.annotation[i].first;
    	//make sure both pointers aren't going to the same object!
    	if (curr_frag1 == NULL && curr_frag2 == NULL) { //null pointers
    		// do nothing
    	}
    	else {
			BOOST_CHECK_NE(curr_frag1,curr_frag2);
			BOOST_CHECK_EQUAL(curr_frag1->name, curr_frag2->name);
			//check to see that breaks are the same
			BOOST_CHECK_EQUAL(spec1.annotation[i].second,spec2.annotation[i].second);
		}
    }
}

BOOST_AUTO_TEST_CASE( readers )  {
	int curr_num_scans;
	SpecSet specs;

	BOOST_TEST_MESSAGE("Testing file readers LoadSpecSet_mgf,LoadSpecSet_pkl,LoadSpecSet_pklbin ");
	BOOST_TEST_MESSAGE("LoadSpecSet_pkl:");
	BOOST_REQUIRE(specs.LoadSpecSet_pkl(pkl_spec));
	curr_num_scans = specs.size();
	BOOST_CHECK_EQUAL(curr_num_scans,num_scans);
	test_spectrum(specs[0]);

	BOOST_TEST_MESSAGE("LoadSpecSet_mgf:");
	BOOST_REQUIRE(specs.LoadSpecSet_mgf(mgf_spec));
	curr_num_scans = specs.size();
	BOOST_CHECK_EQUAL(curr_num_scans,num_scans);
	test_spectrum(specs[0]);

	BOOST_TEST_MESSAGE("LoadSpecSet_pklbin:");
	BOOST_REQUIRE(specs.LoadSpecSet_pklbin(pklbin_spec));
	curr_num_scans = specs.size();
	BOOST_CHECK_EQUAL(curr_num_scans,num_scans);
	test_spectrum(specs[0]);
}

BOOST_AUTO_TEST_CASE( annotation_test ) {
	int curr_num_peaks;
	int i;
	SpecSet specs;
	MS2ScoringModel model;

	BOOST_TEST_MESSAGE("Loading annotation test pkl:");
	BOOST_REQUIRE(specs.LoadSpecSet_pkl(pkl_spec));
	Spectrum test_spec = specs[0];
	curr_num_peaks = test_spec.size();
	BOOST_CHECK_EQUAL(curr_num_peaks,num_peaks);

	BOOST_TEST_MESSAGE("Loading annotation model file:");
	BOOST_REQUIRE(model.LoadModel(model_file));
	BOOST_REQUIRE_EQUAL(model.probs.size(), num_ion_types);
	string include_ions = "all";
	test_spec.annotate(annotation, include_ions, model, 0, 0, .45);
	BOOST_REQUIRE_EQUAL(test_spec.ionTypes.size(),num_ion_types);
	BOOST_REQUIRE_EQUAL(test_spec.annotation.size(),num_peaks);

	BOOST_TEST_MESSAGE("Testing ion annotations against expected:");

	for(i=0;i<curr_num_peaks;i++){
		const ftIonFragment* curr_frag = test_spec.annotation[i].first;
		const short curr_ion = test_spec.annotation[i].second;

		if (ion_index[0][i] != 0) { //ion should  be defined
			if (curr_frag != (ftIonFragment*) NULL) {
				BOOST_CHECK_MESSAGE( ion_annot[0][i].compare(curr_frag->name) == 0 && ion_index[0][i] == curr_ion,
						"Ions do not match! mass: " << test_spec.peakList[i][0] << " correct: " << ion_annot[0][i] << " " << ion_index[0][i] << " incorrect: " << curr_frag->name <<  " " << curr_ion);
			}
			else {
				BOOST_CHECK_MESSAGE(false,"Ion not defined for mass: " << test_spec.peakList[i][0] << " correct: " << ion_annot[0][i] << " " << ion_index[0][i]);
			}
		}
		else { //ion should not be defined
			BOOST_CHECK_MESSAGE(curr_frag == (ftIonFragment*) NULL, "Ion defined incorrectly for mass: " << test_spec.peakList[i][0] << " incorrect: " << curr_frag->name << " " << curr_ion);
		}
	}

	BOOST_TEST_MESSAGE("Testing annotation of spectrum");
	BOOST_CHECK(test_spec.annotation_peptide.compare(annotation) == 0);
}

BOOST_AUTO_TEST_CASE( annotation_start_mod_test ) {
  int curr_num_peaks;
  int i;
  SpecSet specs;
  MS2ScoringModel model;

  BOOST_TEST_MESSAGE("Loading annotation test pkl:");
  BOOST_REQUIRE(specs.LoadSpecSet_pkl(pkl_spec));
  Spectrum test_spec = specs[0];
  curr_num_peaks = test_spec.size();
  BOOST_CHECK_EQUAL(curr_num_peaks,num_peaks);

  BOOST_TEST_MESSAGE("Loading annotation model file:");
  BOOST_REQUIRE(model.LoadModel(model_file));
  BOOST_REQUIRE_EQUAL(model.probs.size(), num_ion_types);
  string include_ions = "all";
  string test_peptide = "*.[-0.3]AKTTPPSV.*";
  BOOST_TEST_MESSAGE("Annotating peptide:");
  test_spec.annotate(test_peptide, include_ions, model, 0, 0, .45);
}

BOOST_AUTO_TEST_CASE(get_scans_test) {
	SpecSet specs;
	Spectrum* curr_spec;

	BOOST_TEST_MESSAGE("Loading annotation test pkl:");
	BOOST_REQUIRE(specs.LoadSpecSet_pkl(pkl_spec));
	BOOST_REQUIRE_EQUAL(specs.size(),2);
	SpecSet specs2;
	specs2.resize(2);

	specs2[0] = specs[1];
	specs2[1] = specs[0];

	BOOST_TEST_MESSAGE("Checking scans: ");
	curr_spec = specs.getScan(1);
	BOOST_CHECK_EQUAL(curr_spec->scan,1);
	BOOST_CHECK_CLOSE(curr_spec->parentMass,parent_mass,.1);

	curr_spec = specs2.getScan(1);
	BOOST_CHECK_EQUAL(curr_spec->scan,1);
	BOOST_CHECK_CLOSE(curr_spec->parentMass,parent_mass,.1);
}

BOOST_AUTO_TEST_CASE(copy_spectrum_test) {
	SpecSet specs;
	MS2ScoringModel model;
	Spectrum spec_copy;


	BOOST_TEST_MESSAGE("Loading annotation test pkl:");
	BOOST_REQUIRE(specs.LoadSpecSet_pkl(pkl_spec));
	BOOST_REQUIRE_EQUAL(specs.size(),2);
	spec_copy = specs[0];

	BOOST_TEST_MESSAGE("Checking copy: ");
	test_copy(spec_copy,specs[0]);

	BOOST_TEST_MESSAGE("Loading annotation model file:");
	BOOST_REQUIRE(model.LoadModel(model_file));
	string include_ions = "all";
	specs[0].annotate(annotation, include_ions, model, 0, 0, .45);

	spec_copy = specs[0];
	BOOST_TEST_MESSAGE("Checking copy with annotate: ");
	test_copy(spec_copy,specs[0]);
}
