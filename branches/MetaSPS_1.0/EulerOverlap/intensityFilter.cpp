/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 4/5/10
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

char* spec_file = "INPUT_SPECS (*.pklbin)";
char* spec_file_mgf = "INPUT_SPECS (*.mgf)";
char* out_file = "OUTPUT_SPECS (*.pklbin)";
char* out_file_mgf = "OUTPUT_SPECS (*.mgf)";
char* activ = "ACTIVATION";
char* pep_file = "INPUT_PEPTIDES (*.txt)";
char* stats_out = "OUTPUT_STATS (*.csv)";

char* peak_tol = "TOLERANCE_PEAK";
char* res = "RESOLUTION";
char* ppm = "TOLERANCE_PPM";
char* rank = "MIN_INTENSITY";

char* aa_masses = "AMINO_ACID_MASSES";

char* BET = ",";

void getMasses(const char *sequence, vector<float> &masses, float offset, InputParams* aa_params) {
    int size = strlen(sequence);
    masses.resize(size+1);
    float f;
    int m;
    int float_i;
    char float_seq[25];
    int spec_index = 1;
    masses[0] = 0.0;
    float total = offset;
    for (unsigned int i = 0; sequence[i] != 0; i++) {
		
        m=0;
        while (m<20 and AAletters[m]!=sequence[i]) m++;
        if(m>=20) {
            if (sequence[i] == '[') {
                float_i = i+1;
                for(; sequence[i] != ']'; i++);
                strncpy(float_seq, &sequence[float_i], i-float_i);
                float_seq[i-float_i] = '\0';
                f = getFloat(float_seq);
                masses[spec_index-1] = masses[spec_index-1] + f;
                total += f;
            }
        }
        else {
			float aa_mass = AAmasses[m];
			if (aa_params != 0) {
				char aa = toupper(sequence[i]);
				aa_mass = (float)(*aa_params).getValueDouble(&aa);
			}
            masses[spec_index] = aa_mass + total;
            total += aa_mass;
            spec_index ++;
        }
    }
    masses.resize(spec_index);
}

void countBY(SpecSet& specs, char* pep_file, float tol, vector<int>& count, InputParams* aa_params) {
	count.resize(8);
	BufferedLineReader orbcal_pep;
	if (! orbcal_pep.Load(pep_file)) {
		cerr << "ERROR: Failed to load " << pep_file << "\n";
	}
	int bcount = 0, ycount = 0, peakcount = 0, specount = 0;
	int bcount2 = 0, ycount2 = 0, peakcount2 = 0, specount2 = 0;
	
	for (int i = 0; i < orbcal_pep.size(); i++) {
		char* next_contig = orbcal_pep.getline(i);
		if (next_contig[0] == '\0' || strlen(next_contig) == 0) continue;
		Spectrum orbcal = specs[i];
		if (orbcal.size() == 0) continue;
		vector<float> perfect_orbcal;
		getMasses(next_contig, perfect_orbcal, 0.0, aa_params);
		if (perfect_orbcal.size() <= 12) {
			specount ++;
			peakcount += orbcal.size();
		} else {
			specount2 ++;
			peakcount2 += orbcal.size();
		}
		
		for (int j = 0; j < perfect_orbcal.size(); j++) {
			float bpeak = perfect_orbcal[j];
			float ypeak = perfect_orbcal[perfect_orbcal.size()-1] - perfect_orbcal[j] + AAJumps::massH2O;
			//cout << bpeak << " " << ypeak << " ";
			int bloc = orbcal.findClosest(bpeak);
			int yloc = orbcal.findClosest(ypeak);
			if (perfect_orbcal.size() <= 12) {
				if (isEqual(orbcal[bloc][0], bpeak, tol)) {bcount ++;}
				if (isEqual(orbcal[yloc][0], ypeak, tol)) {ycount ++;}
			} else {
				if (isEqual(orbcal[bloc][0], bpeak, tol)) {bcount2 ++;}
				if (isEqual(orbcal[yloc][0], ypeak, tol)) {ycount2 ++;}
			}
		}
	}
	count[0] = bcount; count[1] = ycount; count[2] = peakcount; count[3] = specount;
	count[4] = bcount2; count[5] = ycount2; count[6] = peakcount2; count[7] = specount2;
}

void intensityFilter(SpecSet& inspecs, SpecSet& outspecs, float thresh) {
	outspecs.specs.clear();
	outspecs.specs.resize(inspecs.size());
	TwoValues<float> peak;
	for (int idx = 0; idx < inspecs.size(); idx ++) {
		Spectrum inspec = inspecs[idx];
		outspecs[idx] = inspec;
		outspecs[idx].peakList.resize(0);
		int startIdx = 0;
		for (int i = 0; i < inspec.size(); i ++) {
			peak = inspec[i];
			if (peak[1] >= thresh) outspecs[idx].peakList.push_back(peak);
		}
	}
}

int main(int argc, char *argv[], char **envp ) {
	InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("intensityFilter.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "intensityFilter.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	
	InputParams aa_params;
	if (params.paramPresent(aa_masses)) {
		paramsOk = aa_params.readParams(params.getValue(aa_masses));
		if(!paramsOk) {
			cerr << "Error opening parameters file " << params.getValue(aa_masses) << "\n";
			return -1;
		}
	}
	
	bool usePPM = params.paramPresent(ppm);
	float rank_thresh = (float)params.getValueDouble(rank);
	
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
	
	FILE* outputS;
	
	if (!params.paramPresent(rank) && params.paramPresent(pep_file)) {
		float minI = 0, maxI = 0;
		bool first = true;
		for (int idx = 0; idx < input.size(); idx++) {
			for (int i = 0; i < input[idx].size(); i++) {
				if (first) {minI = input[idx][i][1]; maxI = input[idx][i][1]; first = false; continue;}
				else {
					minI = min(minI, input[idx][i][1]);
					maxI = max(maxI, input[idx][i][1]);
				}
			}
		}
		float bestThresh;
		float maxRatio = -1;
		cout << "MinIntensity - #Prefix Masses - #Suffix Masses - #Total Masses Matched - Total Peaks - Total Spectra - Peaks Per Spectrum - Matched peaks/Total peaks\n";
		if (params.paramPresent(stats_out)) {
			outputS = fopen(params.getValue(stats_out), "w");
			fprintf(outputS, "Peptide Length%sMin Intensity%s#Prefix Masses%s#Suffix Masses%s#Total Masses Matched%sTotal Peaks%sTotal Spectra%sPeaks/Spectrum%sMatched Peaks/Total Peaks\n", BET, BET, BET, BET, BET, BET, BET, BET);
		}
		for (float thresh = minI; thresh <= maxI; thresh += (maxI-minI)/20.0) {
			SpecSet output;
			if (thresh == 0) output = input;
			else intensityFilter(input, output, thresh);
			vector<int> count(8);
			if (params.paramPresent(aa_masses)) {
				countBY(output, params.getValue(pep_file), InputParams::PeakTol, count, &aa_params);
			} else {
				countBY(output, params.getValue(pep_file), InputParams::PeakTol, count, 0);
			}
			int allpeaks = count[2] + count[6];
			int numspecs = count[3] + count[7];
			int tot = count[0] + count[1] + count[4] + count[5];
			float ratio = ((float)tot)/((float)allpeaks);
			float ratio1 = ((float)(count[0] + count[1]))/((float)count[2]);
			float ratio2 = ((float)(count[4] + count[5]))/((float)count[6]);
			cout << thresh << " - " << count[0] << " - " << count[1] << " - " << count[0] + count[1] << " - " << count[2] << " - " << count[3] << " - " << ((float)count[2])/((float)count[3]) << " - " << ratio1 << " (<= 12 length peptides)\n";
			cout << thresh << " - " << count[4] << " - " << count[5] << " - " << count[4] + count[5] << " - " << count[6] << " - " << count[7] << " - " << ((float)count[6])/((float)count[7]) << " - " << ratio2 << " (> 12 length peptides)\n";
			if (params.paramPresent(stats_out)) {
				fprintf(outputS, "<= 12%s%.2f%s%d%s%d%s%d%s%d%s%d%s%.2f%s%.2f\n", BET, thresh, BET, count[0], BET, count[1], BET, count[0] + count[1], BET, count[2], BET, count[3], BET, ((float)count[2])/((float)count[3]), BET, ratio1);
				fprintf(outputS, "> 12%s%.2f%s%d%s%d%s%d%s%d%s%d%s%.2f%s%.2f\n", BET, thresh, BET, count[4], BET, count[5], BET, count[4] + count[5], BET, count[6], BET, count[7], BET, ((float)count[6])/((float)count[7]), BET, ratio2);
			}
			if (thresh > 0 && ratio > maxRatio) {
				maxRatio = ratio; bestThresh = thresh; outspecs = output;
			}
		}
		if (params.paramPresent(stats_out)) fclose(outputS);
		cout << "\nSelected output with " << maxRatio << " matched peak ratio and minimum intensity  " << bestThresh << "\n";
	} else if (params.paramPresent(rank) && params.paramPresent(pep_file)) {
		vector<int> count(8);
		float rank_thresh = (float)params.getValueDouble(rank);
		intensityFilter(input, outspecs, rank_thresh);
		if (params.paramPresent(aa_masses)) {
			countBY(outspecs, params.getValue(pep_file), InputParams::PeakTol, count, &aa_params);
		} else {
			countBY(outspecs, params.getValue(pep_file), InputParams::PeakTol, count, 0);
		}
		int allpeaks = count[2] + count[6];
		int numspecs = count[3] + count[7];
		int tot = count[0] + count[1] + count[4] + count[5];
		float ratio = ((float)tot)/((float)allpeaks);
		float ratio1 = ((float)(count[0] + count[1]))/((float)count[2]);
		float ratio2 = ((float)(count[4] + count[5]))/((float)count[6]);
		if (params.paramPresent(stats_out)) {
			outputS = fopen(params.getValue(stats_out), "w");
			fprintf(outputS, "Peptide Length%sMin Intensity%s#Prefix Masses%s#Suffix Masses%s#Total Masses Matched%sTotal Peaks%sTotal Spectra%sPeaks/Spectrum%sMatched Peaks/Total Peaks\n", BET, BET, BET, BET, BET, BET, BET, BET);
		}
		cout << "MinIntensity - #Prefix Masses - #Suffix Masses - #Total Masses Matched - Total Peaks - Total Spectra - Peaks Per Spectrum - Matched peaks/Total peaks\n";
		cout << rank_thresh << " - " << count[0] << " - " << count[1] << " - " << count[0] + count[1] << " - " << count[2] << " - " << count[3] << " - " << ((float)count[2])/((float)count[3]) << " - " << ratio1 << " (<= 12 length peptides)\n";
		cout << rank_thresh << " - " << count[4] << " - " << count[5] << " - " << count[4] + count[5] << " - " << count[6] << " - " << count[7] << " - " << ((float)count[6])/((float)count[7]) << " - " << ratio2 << " (> 12 length peptides)\n";
		if (params.paramPresent(stats_out)) {
			fprintf(outputS, "<= 12%s%.2f%s%d%s%d%s%d%s%d%s%d%s%.2f%s%.2f\n", BET, rank_thresh, BET, count[0], BET, count[1], BET, count[0] + count[1], BET, count[2], BET, count[3], BET, ((float)count[2])/((float)count[3]), BET, ratio1);
			fprintf(outputS, "> 12%s%.2f%s%d%s%d%s%d%s%d%s%d%s%.2f%s%.2f\n", BET, rank_thresh, BET, count[4], BET, count[5], BET, count[4] + count[5], BET, count[6], BET, count[7], BET, ((float)count[6])/((float)count[7]), BET, ratio2);
		}
		if (params.paramPresent(stats_out)) fclose(outputS);
	} else if (params.paramPresent(rank) && !params.paramPresent(pep_file)) {
		float rank_thresh = (float)params.getValueDouble(rank);
		intensityFilter(input, outspecs, rank_thresh);
	} else {cerr<<"Must specify either peptides or threshold!\n"; return -1;}
	
	if (params.paramPresent(out_file)) {
		cout << "Saving to " << params.getValue(out_file) << "\n";
		outspecs.SaveSpecSet_pklbin(params.getValue(out_file));
	}
	if (params.paramPresent(out_file_mgf)) {
		cout << "Saving to " << params.getValue(out_file_mgf) << "\n";
		char* activation = 0;
		if (params.paramPresent(activ)) {activation = params.getValue(activ);}
		outspecs.SaveSpecSet_mgf(params.getValue(out_file_mgf), activation);
	}
	
	if ((!params.paramPresent(out_file)) && (!params.paramPresent(out_file_mgf))) {
		cerr<<"Output file not specified!!\n";
		return -1;
	}
	return 0;
}
