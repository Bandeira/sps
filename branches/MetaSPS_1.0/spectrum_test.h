/*
 * spectrum_test.h
 *
 *  Created on: Sep 1, 2010
 *      Author: jsnedecor
 */

#ifndef SPECTRUM_TEST_H_
#define SPECTRUM_TEST_H_

#include "spectrum.h"

const char* pkl_spec = "../qc/spectrum_test.pkl";
const char* pklbin_spec = "../qc/spectrum_test.pklbin";
const char* mgf_spec = "../qc/spectrum_test.mgf";
char* model_file = "../qc/dancik_model.txt";

//first peak info
const int num_peaks = 67;
const int num_scans = 2;
const float parent_mass = 400.729;
const float first_mass = 118.1;
const float first_intensity = 45.4;
const float last_mass = 766.6;
const float last_intensity = 8.6;
int scan_num = 1;

//model info
const int num_ion_types = 28;

//annotation constants
std::string annotation = "*.AKTTPPSV.*";
//this should probably be less terrible...
const string ion_annot[2][67] = {
		{"y","","","y-H2O","","b"  ,"y","","","b-H2O","y-H2O","","b"  ,"y"    ,"y-iso","","","","b++","","","b-H2O-NH3","","y-H2O","y-NH3"  ,"","b-H2O","b-NH3","P++-H2O"   ,"","","","a-NH3","","b-H2O-H2O","b-H2O-NH3","","b-H2O","b","y"    ,"y-iso","y-H2O-H2O","b-H2O","","y-H2O","y-NH3","b","b-iso","y","y-iso","","","b-H2O","b-NH3","","","b","b-iso","","","","","","","","",""},
		{"" ,"","",""     ,"","y++","" ,"","",""     ,"b-NH3","","y++","b-iso",""     ,"","","",""   ,"","",""         ,"",""     ,""       ,"",""     ,""     ,""          ,"","","",""     ,"",""         ,"y-H2O-H2O","",""     ,"" ,"b-iso",""     ,""         ,""     ,"",""     ,""     ,"" ,""     ,"" ,""     ,"","",""     ,""     ,"","","" ,""     ,"","","","","","","","",""}
};

const int ion_index[2][67] = {
		{1,0,0,2,0,2,2,0,0,3,3,0,3,3,3,0,0,0,7,0,0,4,0,4,4,0,4,4,8,0,0,0,5,0,5,5,0,5,5,5,5,6,6,0,6,6,6,6,6,6,0,0,7,7,0,0,7,7,0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,4,0,0,0,0,3,0,6,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
};

#endif /* SPECTRUM_TEST_H_ */
