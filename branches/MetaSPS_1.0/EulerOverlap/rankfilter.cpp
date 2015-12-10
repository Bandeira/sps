/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 3/16/10
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

#include "../inputParams.h"
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

const char* spec_file = "INPUT_SPECS (*.pklbin)";
const char* spec_file_mgf = "INPUT_SPECS (*.mgf)";
const char* out_file = "OUTPUT_SPECS (*.pklbin)";
const char* out_file_mgf = "OUTPUT_SPECS (*.mgf)";
const char* activ = "ACTIVATION";
const char* pep_file = "INPUT_PEPTIDES (*.txt)";

const char* peak_tol = "TOLERANCE_PEAK";
const char* res = "RESOLUTION";
const char* ppm = "TOLERANCE_PPM";
const char* rank = "RANK";
const char* wind = "WINDOW (default: 50)";

struct SortRank : public std::binary_function<vector<float>, vector<float>, bool> {
	  bool operator()(vector<float> left, vector<float> right) const
	{
		return left[2] > right[2];
	};
};

void countBY(SpecSet& specs, const char* pep_file, float tol, vector<int>& count) {
	count.resize(4);
	BufferedLineReader orbcal_pep;
	AAJumps jumps;
	if (! orbcal_pep.Load(pep_file)) {
		cerr << "ERROR: Failed to load " << pep_file << "\n";
	}
	int bcount = 0, ycount = 0, peakcount = 0, specount = 0;
	
	for (int i = 0; i < orbcal_pep.size(); i++) {
		const char* next_contig = orbcal_pep.getline(i);
		if (next_contig[0] == '\0' || strlen(next_contig) == 0) continue;
		Spectrum orbcal = specs[i];
		if (orbcal.size() == 0) continue;
		specount ++;
		peakcount += orbcal.size();
		vector<float> perfect_orbcal;
		jumps.getPRMMasses(string(next_contig), perfect_orbcal);
		
		for (int j = 0; j < perfect_orbcal.size(); j++) {
			float bpeak = perfect_orbcal[j];
			float ypeak = perfect_orbcal[perfect_orbcal.size()-1] - perfect_orbcal[j] + AAJumps::massH2O;
			int bloc = orbcal.findClosest(bpeak);
			int yloc = orbcal.findClosest(ypeak);
			if (isEqual(orbcal[bloc][0], bpeak, tol)) {bcount ++;}
			if (isEqual(orbcal[yloc][0], ypeak, tol)) {ycount ++;}
		}
	}
	count[0] = bcount; count[1] = ycount; count[2] = peakcount; count[3] = specount;
}

void rankFilter(SpecSet& inspecs, SpecSet& outspecs, int rank_thresh, float window = 50.0) {
	outspecs.specs.clear();
	outspecs.specs.resize(inspecs.size());
	list<vector<float> > ranklist;
	list<vector<float> >::iterator ranklistIt;
	vector<float> peakrank(3);
	TwoValues<float> peak;
	
	for (int idx = 0; idx < inspecs.size(); idx ++) {
		Spectrum inspec = inspecs[idx];
		outspecs[idx] = inspec;
		outspecs[idx].peakList.resize(0);
		int startIdx = 0;
		for (int i = 0; i < inspec.size(); i ++) {
			ranklist.clear();
			float mass = inspec[i][0];
			float intensity = inspec[i][1];
			bool foundRange = false;
			for (int j = startIdx; j < inspec.size(); j++) {
				if (! foundRange && inspec[j][0] >= mass - window) {
					foundRange = true;
					startIdx = j;
				}
				if (inspec[j][0] < mass - window) continue;
				if (inspec[j][0] > mass + window) break;
				peakrank[0] = (float)j;
				peakrank[1] = inspec[j][0];
				peakrank[2] = inspec[j][1];
				ranklist.push_back(peakrank);
			}
			ranklist.sort(SortRank());
			int count = 0;
			for (ranklistIt = ranklist.begin(); ranklistIt != ranklist.end(); ranklistIt ++) {
				if (count >= rank_thresh) break;
				count ++;
				peakrank = *ranklistIt;
				if (((int)(round(peakrank[0])+0.1)) == i) {
					peak[0] = peakrank[1];
					peak[1] = peakrank[2];
					outspecs[idx].peakList.push_back(peak);
					break;
				}
			}
		}
	}
}

int main(int argc, char *argv[], char **envp ) {

float testF = 12345.678910;
double testD = testF;
double testD2 = 12345.678910;

if (testF < testD2) {
printf("less\n");
} else {
printf("geq\n");
}


printf("double: %f\n", testD);
printf("float: %f\n", testF);

	InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("rankfilter.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "rankfilter.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs(1);
	paramStrs[0] = rank;
	
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"rankfilter.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: " << paramStrs[0];
		for (int i = 1; i < paramStrs.size(); i++) cerr << ", " << paramStrs[i];
		cerr << endl << "see rankfilter.params in sps/trunk/EulerOverlap\n";
		return -1;
	}
	
	InputParams::Resolution = getResolution(InputParams::PeakTol);
	
	bool usePPM = params.paramPresent(ppm);
	int rank_thresh = params.getValueInt(rank);
	float window = (params.paramPresent(wind)) ? (float)params.getValueDouble(wind) : 50.0;
	
	SpecSet outspecs;
	SpecSet input;
	if (params.paramPresent(spec_file_mgf)) {
		if (!input.LoadSpecSet_mgf(params.getValue(spec_file_mgf))) return -1;
	} else if (params.paramPresent(spec_file)) {
		if (!input.LoadSpecSet_pklbin(params.getValue(spec_file))) return -1;
	} else {
		cerr<<"Input file not specified!!\n";
		return -1;
	}
	
	if (params.paramPresent(pep_file)) {
		int bestThresh;
		float maxRatio = -1;
		cout << "RankFilter - #B Ions Matched - #Y Ions Matched - Total Ions Matched - Total Peaks - Peaks Per Spectrum - Matched peaks/Total peaks\n";
		for (int thresh = 0; thresh <= 5; thresh ++) {
			SpecSet output;
			if (thresh == 0) output = input;
			else rankFilter(input, output, thresh, window);
			vector<int> count(4);
			countBY(output, params.getValue(pep_file), InputParams::PeakTol, count);
			int allpeaks = count[2];
			int numspecs = count[3];
			int tot = count[0] + count[1];
			float ratio = ((float)tot)/((float)allpeaks);
			cout << thresh << " - " << count[0] << " - " << count[1] << " - " << tot << " - " << allpeaks << " - " << ((float)allpeaks)/((float)numspecs) << " - " << ratio << "\n";
			if (thresh > 0 && ratio > maxRatio) {
				maxRatio = ratio; bestThresh = thresh; outspecs = output;
			}
		}
		cout << "\nSelected output with " << maxRatio << " matched peak ratio and rank threshold " << bestThresh << "\n";
	} else rankFilter(input, outspecs, rank_thresh, window);
		
	
	if (params.paramPresent(out_file)) {
		cout << "Saving to " << params.getValue(out_file) << "\n";
		outspecs.SaveSpecSet_pklbin(params.getValue(out_file));
	}
	if (params.paramPresent(out_file_mgf)) {
		cout << "Saving to " << params.getValue(out_file_mgf) << "\n";
		char* activation = 0;
		if (params.paramPresent(activ)) {outspecs.SaveSpecSet_mgf(params.getValue(out_file_mgf), params.getValue(activ));}
		else {outspecs.SaveSpecSet_mgf(params.getValue(out_file_mgf));}
	}
	
	if ((!params.paramPresent(out_file)) && (!params.paramPresent(out_file_mgf))) {
		cerr<<"Output file not specified!!\n";
		return -1;
	}
	return 0;
}
