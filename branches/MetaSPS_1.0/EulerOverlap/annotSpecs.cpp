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
char* out_file_mgf = "OUTPUT_ANNOT_SPECS (*.mgf)";
char* pep_file = "INPUT_PEPTIDES (*.txt)";

int main(int argc, char *argv[], char **envp ) {
	InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("annotSpecs.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "annotSpecs.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	BufferedLineReader orbcal_pep;
	if (! orbcal_pep.Load(params.getValue(pep_file))) {
		cerr << "ERROR: Failed to load " << pep_file << "\n";
	}
	
	SpecSet input;
	char* infile;
	if (params.paramPresent(spec_file_mgf)) {
		infile = params.getValue(spec_file_mgf);
		if (!input.LoadSpecSet_mgf(infile)) return -1;
	} else if (params.paramPresent(spec_file)) {
		infile = params.getValue(spec_file);
		if (!input.LoadSpecSet_pklbin(infile)) return -1;
	} else {
		cerr<<"Input file not specified!!\n";
		return -1;
	}
	
	FILE* outputS = fopen(params.getValue(out_file_mgf), "w");
	for (int idx = 0; idx < input.size(); idx ++) {
		if (idx >= orbcal_pep.size()) break;
		char* next_contig = orbcal_pep.getline(idx);
		if (next_contig[0] == '\0' || strlen(next_contig) == 0) continue;
		Spectrum orbcal = input[idx];
		if (orbcal.size() == 0) continue;
		fprintf(outputS, "BEGIN IONS\nTITLE=%s_%d\nSEQ=%s\nPEPMASS=%.5f\nCHARGE=+%d\n", infile, idx, next_contig, orbcal.parentMass, orbcal.parentCharge);
		for (int i = 0; i < orbcal.size(); i++) {
			fprintf(outputS, "%.5f\t%.5f\n", orbcal[i][0], orbcal[i][1]);
		}
		fprintf(outputS, "END IONS\n\n");
	}
	fclose(outputS);
	return 0;
}
