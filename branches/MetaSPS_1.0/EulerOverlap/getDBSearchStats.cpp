/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 6/29/09
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
#include "../uncertain.h"
#include "../dbg_print.h"
#include "../twovalues.h"
#include "../db_fasta.h"
#include "../utils.h"
#include "../batch.h"

using namespace std;

const char* SEP1 = ",";

const char* usage = "To output a .csv file of statistics on b/y peaks matched to peptide ids, run:\n\n./getDBSearchStats [input peptides] [input specs.pklbin] [output .csv] [peak tol]\n";

const int num_bins = 20;

TwoValues<float> getBYPeaks(vector<float>& perfect_orbcal, Spectrum spec, float peakTol) {
	float count = 0.0;
	float total = 0.0;
	float pm = perfect_orbcal[perfect_orbcal.size()-1];
	for (int i = 1; i < perfect_orbcal.size(); i++) {
		float bp = perfect_orbcal[i] + AAJumps::massHion;
		float yp = pm - perfect_orbcal[i] + AAJumps::massH2O + AAJumps::massHion;
		float closeb = spec[spec.findClosest(bp)][0];
		float closey = spec[spec.findClosest(yp)][0];
		total += 1.0;
		if (isEqual(closeb, bp, peakTol) || isEqual(closey, yp, peakTol)) {count += 1.0;}
	}
	return TwoValues<float>(count, total);
}

int main(int argc, char *argv[], char **envp ) {
	if (argc != 5) {cout << usage; return 0;}
	SpecSet specs;
	BufferedLineReader peptides;
	if (! specs.LoadSpecSet_pklbin(argv[2])) {cerr << "ERROR: Failed to load " << argv[2] << "\n"; return 0;}
	if (! peptides.Load(argv[1])) {cerr << "ERROR: Failed to load " << argv[1] << "\n"; return 0;}
	float peakTol = getFloat(argv[4]);
	AAJumps jumps;
	//vector<vector<float> > results(specs.specs.size());
	//vector<float> result(3);
	vector<float> perfect_orbcal;
	//for (int i = 0; i < results.size(); i++) results[i] = result;
	cout << "Found " << specs.size() << " spectra and " << peptides.size() << " lines in peptide text file.\n";
	TwoValues<float> res;
	int idxUse1 = 0, idxUse2 = 0, idxUse3 = 0;//, idxUse4 = 0;
	vector<float> data1(specs.specs.size()), data2(specs.specs.size()), data3(specs.specs.size());//, data4;
	vector<TwoValues<float> > binFreq1(num_bins), binFreq2(num_bins), binFreq3(num_bins);//, binFreq4(num_bins);
	int num_ident = 0;
	for (int i = 0; i < specs.size(); i++) {
		if (i >= peptides.size()) break;
		char* next_line = peptides.getline(i);
		if (next_line[0] == '\0' || strlen(next_line) == 0) continue;
		Spectrum spec = specs[i];
		perfect_orbcal.resize(0);
		getMassesCummulative(next_line, perfect_orbcal, 0);
		num_ident ++;
		int size = perfect_orbcal.size();
		
		/*
		if (i == 9616) {
			res = getBYPeaks(perfect_orbcal, spec, peakTol);
			printf("Spectrum %d of size %d has %.1f b/y breaks out of a possible %.1f\n", i, size, res[0], res[1]);
		}*/
		res = getBYPeaks(perfect_orbcal, spec, peakTol);
		
		if (size <= 14) {
			data1[idxUse1] = res[0]/res[1];
			//results[idxUse1][0] = res[0]/res[1];
			idxUse1 ++;
		}
		else if (size <= 25) {
			//res = getBYPeaks(perfect_orbcal, spec, peakTol);
			data2[idxUse2] = res[0]/res[1];
			//results[idxUse2][1] = res[0]/res[1];
			idxUse2 ++;
		} else {
			//res = getBYPeaks(perfect_orbcal, spec, peakTol);
			data3[idxUse3] = res[0]/res[1];
			//results[idxUse3][2] = res[0]/res[1];
			idxUse3 ++;
		}
		/*
		res = getBYPeaks(perfect_orbcal, spec, peakTol);
		data4.push_back(res[0]/res[1]);
		results[idxUse4][3] = res[0]/res[1];
		idxUse4 ++;
		*/
	}
	cout << "Found " << num_ident << " identified spectra\n";
	data1.resize(idxUse1);
	data2.resize(idxUse2);
	data3.resize(idxUse3);
	
	//float bin = 1.0;
	float step = 1.0/((float)num_bins);
	for (int i = 0; i < num_bins; i ++) {
		//binFreq1[i][0] = bin; binFreq2[i][0] = bin; binFreq3[i][0] = bin;
		//bin += 2.0;
		binFreq1[i][0] = step * ((float)i);
		binFreq2[i][0] = step * ((float)i);
		binFreq3[i][0] = step * ((float)i);
	}
	
	getHistogramInfo(data1, binFreq1);
	getHistogramInfo(data2, binFreq2);
	getHistogramInfo(data3, binFreq3);
	//getHistogramInfo(data4, binFreq4);
	int maxIdx = max(idxUse1, idxUse2);
	maxIdx = max(maxIdx, idxUse3);
	//maxIdx = max(maxIdx, idxUse4);
	maxIdx = max(maxIdx, num_bins);
	FILE* output = fopen(argv[3], "w");
	fprintf(output, "Percent Observed b/y Breaks:\n");
	fprintf(output, "0-14%sBins%sFrequency%s15-25%sBins%sFrequency%s>25%sBins%sFrequency\n", SEP1, SEP1, SEP1, SEP1, SEP1, SEP1, SEP1, SEP1);
	for (int i = 0; i < maxIdx; i++) {
		
		if ( i < idxUse1 ) fprintf(output, "%.2f", data1[i]);
		fprintf(output, "%s", SEP1);
		
		if (i < num_bins ) fprintf(output, "%.2f%s%.2f%s", binFreq1[i][0], SEP1, binFreq1[i][1], SEP1);
		else fprintf(output, "%s%s", SEP1, SEP1);
		
		if ( i < idxUse2 ) fprintf(output, "%.2f", data2[i]);
		fprintf(output, "%s", SEP1);
		
		if (i < num_bins ) fprintf(output, "%.2f%s%.0f%s", binFreq2[i][0], SEP1, binFreq2[i][1], SEP1);
		else fprintf(output, "%s%s", SEP1, SEP1);
		
		if ( i < idxUse3 ) fprintf(output, "%.2f", data3[i]);
		fprintf(output, "%s", SEP1);
		
		if (i < num_bins ) fprintf(output, "%.2f%s%.2f%s\n", binFreq3[i][0], SEP1, binFreq3[i][1], SEP1);
		else fprintf(output, "%s%s\n", SEP1, SEP1);

		/*
		if (i < num_bins ) fprintf(output, "%.2f%s%.0f%s", binFreq3[i][0], SEP1, binFreq3[i][1], SEP1);
		else fprintf(output, "%s%s", SEP1, SEP1);
		fprintf(output, "%s", SEP1);
		

		if ( i < idxUse4 ) fprintf(output, "%.2f", results[i][3]);
		fprintf(output, "%s", SEP1);
		
		if (i < num_bins ) fprintf(output, "%.2f%s%.0f\n", binFreq4[i][0], SEP1, binFreq4[i][1]);
		else fprintf(output, "%s\n", SEP1);
		*/
	}
	fclose(output);
	return 0;
}
