/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 7/29/09
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <sstream>
#include <cstring>
#include <list>
#include <vector>
#include <iomanip>

#include "PruneEnvelopes.h"
#include "../ion.h"
#include "../label.h"
#include "../range.h"
#include "../spectrum.h"
#include "../aminoacid.h"
#include "../dbg_print.h"
#include "../twovalues.h"
#include "../db_fasta.h"
#include "../utils.h"
#include "../batch.h"

using namespace std;

const char* spec_file = "INPUT_RAW_SPECS (*.pklbin)";
const char* spec_file_mgf = "INPUT_RAW_SPECS (*mgf)";
const char* dist_file = "INPUT_ISOTOPIC_DISTRIBUTION (*.bin)";
const char* out_file = "OUTPUT_SPECS (*.pklbin)";
const char* out_file_mgf = "OUTPUT_SPECS (*.mgf)";
const char* activ = "ACTIVATION";
const char* pep_file = "INPUT_PEPTIDES (*.txt)";

const char* peak_tol = "TOLERANCE_PEAK";
const char* res = "RESOLUTION";
const char* ppm = "TOLERANCE_PPM";
const char* thresh = "THRESHOLD (default: 0.58)";

TwoValues<int> countBY(SpecSet& specs, const char* pep_file, float tol) {
	BufferedLineReader orbcal_pep;
	if (! orbcal_pep.Load(pep_file)) {
		cerr << "ERROR: Failed to load " << pep_file << "\n";
	}
	int bcount = 0, ycount = 0;
	int specC = 0, peakC = 0;
	
	for (int i = 0; i < specs.size(); i++) {
		const char* next_contig = orbcal_pep.getline(i);
		if (next_contig[0] == '\0' || strlen(next_contig) == 0) continue;

		Spectrum orbcal = specs[i];

		if (orbcal.size() == 0) continue;
		specC ++;
		peakC += orbcal.peakList.size();
		vector<float> perfect_orbcal;
		AAJumps jumps;
		string pepseq = next_contig;
		pepseq = makeBracketsMods(pepseq);
		jumps.getPRMMasses(pepseq, perfect_orbcal, 0.0);
		
		for (int j = 0; j < perfect_orbcal.size(); j++) {
			float bpeak = perfect_orbcal[j] + AAJumps::massHion;
			float ypeak = perfect_orbcal[perfect_orbcal.size()-1] - perfect_orbcal[j] + AAJumps::massH2O + AAJumps::massHion;
			int bloc = orbcal.findClosest(bpeak);
			int yloc = orbcal.findClosest(ypeak);
			if (isEqual(orbcal[bloc][0], bpeak, tol)) {bcount ++;}
			if (isEqual(orbcal[yloc][0], ypeak, tol)) {ycount ++;}
		}
	}
	cout << specC << " spectra and " << peakC << " peaks\n";
	return TwoValues<int>(bcount, ycount);
}

int main(int argc, char *argv[], char **envp ) {
	InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("prunenv.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "prunenv.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs(1);
	paramStrs[0] = dist_file;
	
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"prunenv.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: " << paramStrs[0];
		for (int i = 1; i < paramStrs.size(); i++) cerr << ", " << paramStrs[i];
		cerr << endl << "see prunenv.params in sps/trunk/EulerOverlap\n";
		return -1;
	}
	
	InputParams::Resolution = getResolution(InputParams::PeakTol);
	
	bool usePPM = params.paramPresent(ppm);
	float threshold = (params.paramPresent(thresh)) ? (float)params.getValueDouble(thresh) : 0.58;
	
	PruneEnvelopes penv(params.getValue(dist_file));
	
	if (params.paramPresent(pep_file)) {
		SpecSet input;
		if (params.paramPresent(spec_file_mgf)) {
			input.LoadSpecSet_mgf(params.getValue(spec_file_mgf));
		} else if (params.paramPresent(spec_file)) {
			input.LoadSpecSet_pklbin(params.getValue(spec_file));
		} else {
			cerr<<"Input file not specified!!\n";
			return 0;
		}
		penv.pklSpecs = input;
		float minThresh = threshold - 0.1;
		float maxThresh = threshold + 0.1;
		float stepThresh = 0.01;
		float bestThresh;
		int maxPeaks = -1;
		SpecSet bestOut;
		short charge_c = 2;
		cout << "Min ParentCharge - #B Ions Matched - #Y Ions Matched\n";
		//for (float thresh = minThresh; thresh <= maxThresh; thresh += stepThresh) {
		for (short charge_c = 2; charge_c <= 7; charge_c += 1) {
		  penv.pklSpecs = input;

			for (int i = 0; i < penv.pklSpecs.size(); i++) {
				if (penv.pklSpecs[i].parentCharge != charge_c) {
				  penv.pklSpecs[i].resize(0);
				}

			}
			TwoValues<int> count = countBY(penv.pklSpecs, params.getValue(pep_file), InputParams::PeakTol);
			cout <<  " - " << charge_c << " - "<< count[0] << " - " << count[1] << " - " << "\n";

			penv.transferAllSpecs(InputParams::PeakTol, threshold, false, usePPM);

			count = countBY(penv.prunedSpecs, params.getValue(pep_file), InputParams::PeakTol);

			cout << " - " << charge_c << " - "<< count[0] << " - " << count[1] << " - " << "\n";

		}
		penv.prunedSpecs = bestOut;
		cout << "Selected output with " << maxPeaks << " ions and threshold " << bestThresh << "\n";
	} else {
		if (params.paramPresent(spec_file_mgf)) {
			penv.uploadMGF(params.getValue(spec_file_mgf), InputParams::PeakTol, threshold, false, usePPM);
		} else if (params.paramPresent(spec_file)) {
			penv.uploadPklBin(params.getValue(spec_file), InputParams::PeakTol, threshold, false, usePPM);
		} else {
			cerr<<"Input file not specified!!\n";
			return 0;
		}
	}
	
	cout << "Pruned " << penv.prunedSpecs.size() << " spectra, \n";
	if (params.paramPresent(out_file)) {
		cout << "Saving to " << params.getValue(out_file) << "\n";
		penv.prunedSpecs.SaveSpecSet_pklbin(params.getValue(out_file));
	}
	if (params.paramPresent(out_file_mgf)) {
		cout << "Saving to " << params.getValue(out_file_mgf) << "\n";
		if (params.paramPresent(activ)) {penv.prunedSpecs.SaveSpecSet_mgf(params.getValue(out_file_mgf), params.getValue(activ));}
		else {penv.prunedSpecs.SaveSpecSet_mgf(params.getValue(out_file_mgf), 0);}
	}
	
	if ((!params.paramPresent(out_file)) && (!params.paramPresent(out_file_mgf))) {
		cerr<<"Output file not specified!!\n";
	}
	return 0;
}
