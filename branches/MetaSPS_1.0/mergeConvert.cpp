
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>
#include <list>
#include <vector>

#include "spectrum.h"
#include "inputParams.h"
#include "utils.h"
#include "abruijn.h"

using namespace std;

const char* in_pklbin = "INPUT_PKLBIN (*:*)";
const char* in_mgf = "INPUT_MGF (*:*)";
const char* in_pkl = "INPUT_PKL (*:*)";
const char* in_pkl_dir = "INPUT_PKL_DIR";
const char* in_mic_dir = "INPUT_MIC_DIR";
const char* in_ms2 = "INPUT_MS2 (*:*)";
const char* in_prms = "INPUT_PRMS (*:*)";
const char* in_prmsv3 = "INPUT_PRMSV3 (*:*)";
const char* in_abinfo = "INPUT_ABINFO (*:*)";
const char* in_overlaps = "INPUT_OVERLAPS (*_midx.pklbin) (*:*)";
const char* in_prot_match = "INPUT_MATCHED_PROTS (*_mp.bin) (*:*)";

const char* out_pklbin = "OUTPUT_PKLBIN";
const char* activ = "ACTIVATION";
const char* out_mgf = "OUTPUT_MGF";
const char* out_pkl = "OUPUT_PKL";
const char* out_ms2 = "OUTPUT_MS2";
const char* out_abinfo = "OUTPUT_ABINFO";
const char* out_overlaps = "OUTPUT_OVERLAPS (*_midx.pklbin)";
const char* out_prot_match = "OUTPUT_MATCHED_PROTS (*_mp.bin)";

int main(int argc, char *argv[], char **envp ) {
	InputParams params;
	bool paramsOk;
	const char* inputParamFile;
		if(argc<=1) {
		paramsOk=params.readParams("mergeConvert.params");
		inputParamFile = "mergeConvert.params";
	} else {
		paramsOk=params.readParams(argv[1]);
		inputParamFile = argv[1];
	}
	if(!paramsOk) {
		cerr << "Error opening parameters file " << inputParamFile << "\n";
		return -1;
	}
	
	SpecSet allSpecs;
	abinfo_t all_abinfo;
	vector<unsigned int> in_spec_sizes;
	SpecSet overlaps;
	vector<vector<int> > prot_match;
	
	ostringstream outs;
	outs << inputParamFile << ".info";
	
	FILE* output = fopen( outs.str().c_str(), "w");
	
	fprintf(output, "Spectra Merge Order:\n<filename>\t<# spectra>\n");
	// if (! Load_binArray(argv[5], prot_match1)) cerr << "ERROR: Failed to load " << argv[5] << "\n";

	int idxUse = 0;
	int idxOv = 0;
	int idxPr = 0;
	vector<string> filenames;
	if (params.paramPresent(in_pklbin)) {
		cout << "Merging pklbin spectra:\n";
		splitText(params.getValue(in_pklbin), filenames, ":");
		unsigned int prev = 0;
		for (int i = 0; i < filenames.size(); i++) {
			SpecSet in_specs;
			if (!in_specs.LoadSpecSet_pklbin(filenames[i].c_str())) {cerr<<"ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			in_spec_sizes.push_back(prev);
			prev += in_specs.size();
			cout << in_specs.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
			allSpecs.specs.resize(allSpecs.size() + in_specs.size());
			for (int j = 0; j < in_specs.size(); j++) {
				allSpecs[idxUse] = in_specs[j];
				idxUse ++;
			}
		}
		cout << "\n\n";
	}

	filenames.resize(0);
	if (params.paramPresent(in_mgf)) {
		cout << "Merging mgf spectra:\n";
		splitText(params.getValue(in_mgf), filenames, ":");
		for (int i = 0; i < filenames.size(); i++) {
			SpecSet in_specs;
			if (!in_specs.LoadSpecSet_mgf(filenames[i].c_str())) {cerr<<"ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			cout << in_specs.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
			allSpecs.specs.resize(allSpecs.size() + in_specs.size());
			for (int j = 0; j < in_specs.size(); j++) {
				allSpecs[idxUse] = in_specs[j];
				idxUse ++;
			}
		}
		cout << "\n\n";
	}
	
	filenames.resize(0);
	if (params.paramPresent(in_pkl)) {
		cout << "Merging pkl spectra:\n";
		splitText(params.getValue(in_pkl), filenames, ":");
		for (int i = 0; i < filenames.size(); i++) {
			SpecSet in_specs;
			if (!in_specs.LoadSpecSet_pkl(filenames[i].c_str())) {cerr<<"ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			cout << in_specs.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
			allSpecs.specs.resize(allSpecs.size() + in_specs.size());
			for (int j = 0; j < in_specs.size(); j++) {
				allSpecs[idxUse] = in_specs[j];
				idxUse ++;
			}
		}
		cout << "\n\n";
	}
	
	if (params.paramPresent(in_pkl_dir)) {
		list<string> dir_files;
		cout << "Merging pkl spectra from directory " << params.getValue(in_pkl_dir) << ":\n";
		if (! getdir(params.getValue(in_pkl_dir), dir_files)) {cerr << "ERROR: Failed to load <" << params.getValue(in_pkl_dir) << ">\n"; return -1;}
		int specCount = 0;
		for (list<string>::iterator dirIt = dir_files.begin(); dirIt != dir_files.end(); dirIt ++) {
			SpecSet in_specs;
			string filename = *dirIt;
			if (filename.find(".pkl") == string::npos) continue;
			string root = params.getValue(in_pkl_dir);
			root.append("/");
			root.append(filename);
			if (!in_specs.LoadSpecSet_pkl_mic(root.c_str())) {cerr<<"ERROR: Failed to load <" << root << ">\n"; return -1;}
			fprintf(output, "%s\t%d\n", root.c_str(), in_specs.size());
			allSpecs.specs.resize(allSpecs.size() + in_specs.size());
			specCount += in_specs.size();
			for (int j = 0; j < in_specs.size(); j++) {
				allSpecs[idxUse] = in_specs[j];
				idxUse ++;
			}
		}
		cout << specCount << " from .pkl files in " << params.getValue(in_pkl_dir) << "\n\n";
	}
	
	if (params.paramPresent(in_mic_dir)) {
		list<string> dir_files;
		cout << "Merging mic spectra from directory " << params.getValue(in_mic_dir) << ":\n";
		if (! getdir(params.getValue(in_mic_dir), dir_files)) {cerr << "ERROR: Failed to load <" << params.getValue(in_mic_dir) << ">\n"; return -1;}
		int specCount = 0;
		for (list<string>::iterator dirIt = dir_files.begin(); dirIt != dir_files.end(); dirIt ++) {
			SpecSet in_specs;
			string filename = *dirIt;
			if (filename.find(".mic") == string::npos) continue;
			string root = params.getValue(in_mic_dir);
			root.append("/");
			root.append(filename);
			if (!in_specs.LoadSpecSet_pkl_mic(root.c_str())) {cerr<<"ERROR: Failed to load <" << root << ">\n"; return -1;}
			fprintf(output, "%s\t%d\n", root.c_str(), in_specs.size());
			allSpecs.specs.resize(allSpecs.size() + in_specs.size());
			specCount += in_specs.size();
			for (int j = 0; j < in_specs.size(); j++) {
				allSpecs[idxUse] = in_specs[j];
				idxUse ++;
			}
		}
		cout << specCount << " from .mic files in " << params.getValue(in_mic_dir) << "\n\n";
	}

	filenames.resize(0);
	if (params.paramPresent(in_ms2)) {
		cout << "Merging ms2 spectra:\n";
		splitText(params.getValue(in_ms2), filenames, ":");
		for (int i = 0; i < filenames.size(); i++) {
			SpecSet in_specs;
			if (!in_specs.LoadSpecSet_ms2(filenames[i].c_str())) {cerr<<"ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			cout << in_specs.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
			allSpecs.specs.resize(allSpecs.size() + in_specs.size());
			for (int j = 0; j < in_specs.size(); j++) {
				allSpecs[idxUse] = in_specs[j];
				idxUse ++;
			}
		}
		cout << "\n\n";
	}

	filenames.resize(0);
	if (params.paramPresent(in_prms)) {
		cout << "Merging prms spectra:\n";
		splitText(params.getValue(in_prms), filenames, ":");
		for (int i = 0; i < filenames.size(); i++) {
			SpecSet in_specs;
			if (!in_specs.LoadSpecSet_prms(filenames[i].c_str())) {cerr<<"ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			cout << in_specs.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
			allSpecs.specs.resize(allSpecs.size() + in_specs.size());
			for (int j = 0; j < in_specs.size(); j++) {
				allSpecs[idxUse] = in_specs[j];
				idxUse ++;
			}
		}
		cout << "\n\n";
	}

	filenames.resize(0);
	if (params.paramPresent(in_prmsv3)) {
		cout << "Merging prmsv3 spectra:\n";
		splitText(params.getValue(in_prmsv3), filenames, ":");
		for (int i = 0; i < filenames.size(); i++) {
			SpecSet in_specs;
			if (!in_specs.LoadSpecSet_prmsv3(filenames[i].c_str())) {cerr<<"ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			cout << in_specs.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
			allSpecs.specs.resize(allSpecs.size() + in_specs.size());
			for (int j = 0; j < in_specs.size(); j++) {
				allSpecs[idxUse] = in_specs[j];
				idxUse ++;
			}
		}
		cout << "\n\n";
	}
	
	filenames.resize(0);
	if (params.paramPresent(in_abinfo)) {
		fprintf(output, "\nAbinfo Merge Order:\n<filename>\t<# ab components>\n");
		cout << "Merging abinfo:\n";
		splitText(params.getValue(in_abinfo), filenames, ":");
		for (int i = 0; i < filenames.size(); i++) {
			if (i >= in_spec_sizes.size()) {
				cerr << "ERROR: Not enough .pklbin files loaded to determine spectrum index offsets for " << filenames[i] << "\n";
				return -1;
			}
			abinfo_t in_ab;
			if (!Load_abinfo(filenames[i].c_str(), in_ab)) {cerr<<"ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			cout << in_ab.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%lu\n", filenames[i].c_str(), (unsigned long)in_ab.size());
			int size_in = all_abinfo.size();
			for (abinfo_t::iterator compit = in_ab.begin(); compit != in_ab.end(); compit++) {
				pair<pair< vector<int>, vector<int> >, vector<pair<vector<int>, vector<double> > > > second = compit->second;
				for (int j = 0; j < second.first.first.size(); j++) {second.first.first[j] += in_spec_sizes[i];}
				for (int j = 0; j < second.second.size(); j++) {
					pair<vector<int>, vector<double> > vert = second.second[j];
					for (int k = 0; k < vert.first.size(); k++) second.second[j].first[k] += in_spec_sizes[i];
				}
				all_abinfo[compit->first + size_in] = second;
			}
		}
		cout << "\n\n";
	}
	
	filenames.resize(0);
	if (params.paramPresent(in_overlaps)) {
		fprintf(output, "\nOverlaps Merge Order:\n<filename>\t<# overlaps>\n");
		cout << "Merging overlaps:\n";
		splitText(params.getValue(in_overlaps), filenames, ":");
		for (int i = 0; i < filenames.size(); i++) {
			SpecSet in_specs;
			if (!in_specs.LoadSpecSet_pklbin(filenames[i].c_str())) {cerr<<"ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			cout << in_specs.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%d\n", filenames[i].c_str(), in_specs.size());
			overlaps.specs.resize(overlaps.size() + in_specs.size());
			for (int j = 0; j < in_specs.size(); j++) {
				overlaps[idxOv] = in_specs[j];
				idxOv ++;
			}
		}
		cout << "\n\n";
	}

	filenames.resize(0);
	if (params.paramPresent(in_prot_match)) {
		fprintf(output, "\nMatched Protein Order:\n<filename>\t<# matched spectra>\n");
		cout << "Merging matched proteins:\n";
		splitText(params.getValue(in_prot_match), filenames, ":");
		for (int i = 0; i < filenames.size(); i ++) {
			vector<vector<int> > prot_m;
			if (! Load_binArray(filenames[i].c_str(), prot_m)) {cerr << "ERROR: Failed to load <" << filenames[i] << ">\n"; return -1;}
			cout << prot_m.size() << " from " << filenames[i] << ", ";
			fprintf(output, "%s\t%lu\n", filenames[i].c_str(), (unsigned long)prot_m.size());
			prot_match.resize(prot_match.size() + prot_m.size());
			for (int j = 0; j < prot_m.size(); j++) {
				prot_match[idxPr] = prot_m[j];
				idxPr ++;
			}
		}
		cout << "\n\n";
	}
	
	if (params.paramPresent(out_mgf) || params.paramPresent(out_pklbin) || params.paramPresent(out_pkl) || params.paramPresent(out_ms2)) {
		fprintf(output, "\nOutput Spectra:\n<filename>\t<# spectra>\n");
	}

	if (params.paramPresent(out_mgf)) {
		cout << "Saving " << allSpecs.size() << " spectra to <" << params.getValue(out_mgf) << "> ...\n";
		const char* activation;
		if (params.paramPresent(activ)) {activation = params.getValue(activ);}
		else {activation = 0;}
		if (! allSpecs.SaveSpecSet_mgf(params.getValue(out_mgf), activation)) {cerr<<"ERROR: Failed to save to <" << params.getValue(out_mgf) << ">\n"; return -1;}
		fprintf(output, "%s\t%d\n", params.getValue(out_mgf), allSpecs.size());
	}

	if (params.paramPresent(out_pklbin)) {
		cout << "Saving " << allSpecs.size() << " spectra to <" << params.getValue(out_pklbin) << "> ...\n";
		if (! allSpecs.SaveSpecSet_pklbin(params.getValue(out_pklbin))) {cerr<<"ERROR: Failed to save to <" << params.getValue(out_pklbin) << ">\n"; return -1;}
		fprintf(output, "%s\t%d\n", params.getValue(out_pklbin), allSpecs.size());
	}

	if (params.paramPresent(out_pkl)) {
		cout << "Saving " << allSpecs.size() << " spectra to <" << params.getValue(out_pkl) << "> ...\n";
		if (! allSpecs.SaveSpecSet_pkl(params.getValue(out_pkl))) {cerr<<"ERROR: Failed to save to <" << params.getValue(out_pkl) << ">\n"; return -1;}
		fprintf(output, "%s\t%d\n", params.getValue(out_pkl), allSpecs.size());
	}

	if (params.paramPresent(out_ms2)) {
		cout << "Saving " << allSpecs.size() << " spectra to <" << params.getValue(out_ms2) << "> ...\n";
		if (! allSpecs.SaveSpecSet_ms2(params.getValue(out_ms2))) {cerr<<"ERROR: Failed to save to <" << params.getValue(out_ms2) << ">\n"; return -1;}
		fprintf(output, "%s\t%d\n", params.getValue(out_ms2), allSpecs.size());
	}
	
	if (params.paramPresent(out_abinfo)) {
		fprintf(output, "\nOutput Abinfo:\n<filename>\t<# ab components>\n");
		cout << "Saving " << all_abinfo.size() << " ab components to <" << params.getValue(out_abinfo) << "> ...\n";
		if (! Save_abinfo_v1_0(params.getValue(out_abinfo), all_abinfo)) {cerr<<"ERROR: Failed to save to <" << params.getValue(out_abinfo) << ">\n"; return -1;}
		fprintf(output, "%s\t%lu\n", params.getValue(out_abinfo), (unsigned long)all_abinfo.size());
	}

	if (params.paramPresent(out_overlaps)) {
		fprintf(output, "\nOutput Overlaps:\n<filename>\t<# overlaps>\n");
		cout << "Saving " << overlaps.size() << " overlaps to <" << params.getValue(out_overlaps) << "> ...\n";
		if (! overlaps.SaveSpecSet_pklbin(params.getValue(out_overlaps))) {cerr<<"ERROR: Failed to save to <" << params.getValue(out_overlaps) << ">\n"; return -1;}
		fprintf(output, "%s\t%d\n", params.getValue(out_overlaps), overlaps.size());
	}

	if (params.paramPresent(out_prot_match)) {
		fprintf(output, "\nOutput Matched Proteins:\n<filename>\t<# spectra matched>\n");
		cout << "Saving " << prot_match.size() << " protein matched spectra to <" << params.getValue(out_prot_match) << "> ...\n";
		if (! Save_binArray(params.getValue(out_prot_match), prot_match)) {cerr<<"ERROR: Failed to save to <" << params.getValue(out_prot_match) << ">\n"; return -1;}
		fprintf(output, "%s\t%lu\n", params.getValue(out_prot_match), (unsigned long)prot_match.size());
	}
	fclose(output);
	return 0;
}

