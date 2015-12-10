/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 10/27/09
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

#include "OverlapOut.h"

using namespace std;

const char* contig_file = "INPUT_CONTIG_SPECTRA (*.pklbin)";
const char* overlaps_file = "INPUT_MATCHED_PEAKS_IDX (*midx.pklbin)";
const char* prot_match_file = "INPUT_MATCHED_PROTS (*mp.bin)";
const char* contig_peak_tol = "CONTIG_TOLERANCE_PEAK";

const char* contig_comps = "INPUT_CONTIG_COMPONENTS";

const char* connectors_file = "INPUT_CONNECTOR_SPECTRA (*.pklbin)";
const char* connector_peptides = "INPUT_CONNECTOR_PEPTIDES";
const char* prec_charge = "MINIMUM_CONNECTOR_CHARGE";
const char* connector_peak_tol = "CONNECTOR_TOLERANCE_PEAK";

const char* min_cont_orb_match = "MINIMUM_CONTIG_CONNECTOR_MP";

const char* fasta_file = "INPUT_FASTA";
const char* protein_idx = "PROTEIN_INDEX";
const char* line_length = "AA_PER_LINE";
const char* delimeter = "DELIMETER";
const char* outfile = "OUTPUT_FILE (*.csv)";

int main(int argc, char *argv[], char **envp ) {
	InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("spsoverlap.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "spsoverlap.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<const char *> paramStrs(5);
	paramStrs[0] = fasta_file;
	paramStrs[1] = outfile;
	paramStrs[2] = protein_idx;
	paramStrs[3] = line_length;
	paramStrs[4] = delimeter;
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"spsoverlap.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: " << paramStrs[0];
		for (int i = 1; i < paramStrs.size(); i++) cerr << ", " << paramStrs[i];
		cerr << endl << "see spsoverlap.params in sps/trunk/EulerOverlap\n";
		return -1;
	}
	
	const char* fasta = params.getValue(fasta_file);
	const char* output = params.getValue(outfile);
	const char* protein = params.getValue(protein_idx);
	int aa_per_line = (int) params.getValueInt(line_length);
	const char* contigs;
	const char* overlaps;
	const char* prot_match;
	const char* connectors;
	const char* peptides;
	const char* comps;
	const char* delim = params.getValue(delimeter);
	int precCharge, cont_orb_mp;
	float contig_pktol, connector_pktol;
	
	short overType;
	
	//OverlapOut(char* over_file, char* contigs_file, char* protein_match, char* fasta_file, char* delim = ",");
	//OverlapOut(char* over_file, char* contigs_file, char* protein_match, char* fasta_file, char* orbcal_seq, char* orbcal_spec, char* delim = ",");
	if (params.paramPresent(contig_file) && !params.paramPresent(contig_comps) && !params.paramPresent(connectors_file)) {
		contigs = params.getValue(contig_file);
		overlaps = params.getValue(overlaps_file);
		prot_match = params.getValue(prot_match_file);
		contig_pktol = (float) params.getValueDouble(contig_peak_tol);
		InputParams::Resolution = getResolution(contig_pktol);
		OverlapOut out(overlaps, contigs, prot_match, fasta, delim);
		ostringstream outs;
		char fileroot[strlen(output)];
		char* spot;
		strcpy(fileroot, output);
		spot = strstr(fileroot, ".csv");
		if (spot != NULL) fileroot[&spot[0] - &fileroot[0]] = 0;
		if (strcmp(protein, "all") == 0) {
			for (int i = 0; i < out.fasta.size(); i++) {
				outs << fileroot << i << ".csv";
				out.outputOverlap(outs.str().c_str(), i, contig_pktol, aa_per_line);
				outs.str("");
			}
		} else {
			out.outputOverlap(strcat(fileroot, ".csv"), getInt(protein), contig_pktol, aa_per_line);
		}
		overType = 0;
	} else if (params.paramPresent(contig_file) && params.paramPresent(contig_comps) && !params.paramPresent(connectors_file)) {
		contigs = params.getValue(contig_file);
		overlaps = params.getValue(overlaps_file);
		prot_match = params.getValue(prot_match_file);
		comps = params.getValue(contig_comps);
		contig_pktol = (float) params.getValueDouble(contig_peak_tol);
		InputParams::Resolution = getResolution(contig_pktol);
		OverlapOut out(overlaps, contigs, prot_match, fasta, delim);
		ostringstream outs;
		char fileroot[strlen(output)];
		char* spot;
		strcpy(fileroot, output);
		spot = strstr(fileroot, ".csv");
		if (spot != NULL) fileroot[&spot[0] - &fileroot[0]] = 0;
		if (strcmp(protein, "all") == 0) {
			for (int i = 0; i < out.fasta.size(); i++) {
				outs << fileroot << i << ".csv";
				out.outputOverlapComp(outs.str().c_str(), comps, i, contig_pktol, aa_per_line);
				outs.str("");
			}
		} else {
			out.outputOverlapComp(strcat(fileroot, ".csv"), comps, getInt(protein), contig_pktol, aa_per_line);
		}
		overType = 3;
	} else if (params.paramPresent(contig_file) && params.paramPresent(contig_comps) && params.paramPresent(connectors_file)) {
		contigs = params.getValue(contig_file);
		overlaps = params.getValue(overlaps_file);
		prot_match = params.getValue(prot_match_file);
		comps = params.getValue(contig_comps);
		connectors = params.getValue(connectors_file);
		peptides = params.getValue(connector_peptides);
		precCharge = params.getValueInt(prec_charge);
		connector_pktol = (float) params.getValueDouble(connector_peak_tol);
		contig_pktol = (float) params.getValueDouble(contig_peak_tol);
		InputParams::Resolution = getResolution(connector_pktol);
		OverlapOut out(overlaps, contigs, prot_match, fasta, peptides, connectors, delim);
		ostringstream outs;
		char fileroot[strlen(output)];
		char* spot;
		strcpy(fileroot, output);
		spot = strstr(fileroot, ".csv");
		if (spot != NULL) fileroot[&spot[0] - &fileroot[0]] = 0;
		if (strcmp(protein, "all") == 0) {
			for (int i = 0; i < out.fasta.size(); i++) {
				outs << fileroot << i << ".csv";
				
				out.outputOrbcalOverlapComp(outs.str().c_str(), comps, i, contig_pktol, (short)precCharge, connector_pktol, aa_per_line);
				outs.str("");
			}
		} else {
			out.outputOrbcalOverlapComp(strcat(fileroot, ".csv"), comps, getInt(protein), contig_pktol, (short)precCharge, connector_pktol, aa_per_line);
		}
		overType = 2;
	} else {
		contigs = params.getValue(contig_file);
		overlaps = params.getValue(overlaps_file);
		prot_match = params.getValue(prot_match_file);
		connectors = params.getValue(connectors_file);
		peptides = params.getValue(connector_peptides);
		precCharge = params.getValueInt(prec_charge);
		connector_pktol = (float) params.getValueDouble(connector_peak_tol);
		contig_pktol = (float) params.getValueDouble(contig_peak_tol);
		InputParams::Resolution = getResolution(connector_pktol);
		OverlapOut out(overlaps, contigs, prot_match, fasta, peptides, connectors, delim);
		ostringstream outs;
		char fileroot[strlen(output)];
		char* spot;
		strcpy(fileroot, output);
		spot = strstr(fileroot, ".csv");
		if (spot != NULL) fileroot[&spot[0] - &fileroot[0]] = 0;
		if (strcmp(protein, "all") == 0) {
			for (int i = 0; i < out.fasta.size(); i++) {
				outs << fileroot << i << ".csv";
				
				out.outputOrbcalOverlap(outs.str().c_str(), i, (short)precCharge, contig_pktol, connector_pktol, aa_per_line);
				outs.str("");
			}
		} else {
			out.outputOrbcalOverlap(strcat(fileroot, ".csv"), getInt(protein), (short)precCharge, contig_pktol, connector_pktol, aa_per_line);
		}
		overType = 1;
	}
	
	return 0;
}
