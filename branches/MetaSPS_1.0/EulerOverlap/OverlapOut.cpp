/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 6/29/09
 */

#include "OverlapOut.h"

using namespace std;

// comparison operator used to sort a list of twovalues by their second then first value.
struct SortOverlap : public std::binary_function<TwoValues<int>, TwoValues<int>, bool> {
	  bool operator()(TwoValues<int> left, TwoValues<int> right) const
	{
		if (left[1] != right[1]) return left[1] < right[1];
		else return left[0] < right[0];
	};
};

OverlapOut::OverlapOut() {}
OverlapOut::~OverlapOut() {}

OverlapOut::OverlapOut(const char* over_file, const char* contigs_file, const char* protein_match, const char* fasta_file, const char* delim) {
	if (! overlaps.LoadSpecSet_pklbin(over_file)) {
		cerr << "ERROR: Failed to load " << over_file << "\n";
	}
	if (! contigs.LoadSpecSet_pklbin(contigs_file)) {
		cerr << "ERROR: Failed to load " << contigs_file << "\n";
	}
	if (! Load_binArray(protein_match, prot_match)) {
		cerr << "ERROR: Failed to load " << protein_match << "\n";
	}
	if (! fasta.Load(fasta_file)) {
		cerr << "ERROR: Failed to load " << fasta_file << "\n";
	}
	DELIM = delim;
	contigSeq.resize(contigs.size());
	orbcalSeq.resize(orbcals.size());
	usesOrbcal = 0;
	checkComps = false;
}

OverlapOut::OverlapOut(const char* over_file, const char* contigs_file, const char* protein_match, const char* fasta_file, const char* orbcal_seq, const char* orbcal_spec, const char* delim) {
	if (! overlaps.LoadSpecSet_pklbin(over_file)) {
		cerr << "ERROR: Failed to load " << over_file << "\n";
	}
	//cout << "\noverlaps: " << overlaps.specs.size() << endl;
	if (! contigs.LoadSpecSet_pklbin(contigs_file)) {
		cerr << "ERROR: Failed to load " << contigs_file << "\n";
	}
	//cout << "contigs: " << contigs.specs.size() << endl;
	if (! Load_binArray(protein_match, prot_match)) {
		cerr << "ERROR: Failed to load " << protein_match << "\n";
	}
	//cout << "protein_match: " << prot_match.size() << endl;
	if (! fasta.Load(fasta_file)) {
		cerr << "ERROR: Failed to load " << fasta_file << "\n";
	}
	//cout << "fasta: " << fasta.size() << endl;
	if (! orbcal_pep.Load(orbcal_seq)) {
		cerr << "ERROR: Failed to load " << orbcal_seq << "\n";
	}
	//cout << "orbcal_pep: " << orbcal_pep.size() << endl;
	if (! orbcals.LoadSpecSet_pklbin(orbcal_spec)) {
		cerr << "ERROR: Failed to load " << orbcal_spec << "\n";
	}
	DELIM = delim;
	contigSeq.resize(contigs.size());
	orbcalSeq.resize(orbcals.size());
	//cout << "orbcals: " << orbcals.specs.size() << endl;
	usesOrbcal = 1;
	checkComps = false;
}

OverlapOut &OverlapOut::operator=(const OverlapOut& other) {
	overlaps = other.overlaps;
	contigs = other.contigs;
	orbcals = other.orbcals;
	fasta = other.fasta;
	prot_match = other.prot_match;
	usesOrbcal = other.usesOrbcal;
	proteinToContig = other.proteinToContig;
	proteinToOrbcal = other.proteinToOrbcal;
	mySet = other.mySet;
	DELIM = other.DELIM;
	contigSeq = other.contigSeq;
	orbcalSeq = other.orbcalSeq;
	checkComps = other.checkComps;
	return(*this);
}

bool OverlapOut::outputOverlap(const char* filename, int proteinIndex, float range, int line_length) {return outputOverlap(filename, proteinIndex, 0, range, 0.0, 0, 0, false, line_length);}

bool OverlapOut::outputOverlapComp(const char* filename, const char* componentFile, int proteinIndex, float range, int line_length) {
	BufferedLineReader comps;
	comps.Load(componentFile);
	contigComponents.clear();
	checkComps = true;
	char* next;
	char* result;
	vector<int> v(10);
	for (int i = 1; i < comps.size(); i++) {
		next = comps.getline(i);
		if (strlen( next ) == 0) continue;
		result = strtok( next, "," );
		int x = 0;
		while( result != NULL ) {
			v[x] = getInt(result);
			result = strtok( NULL, "," );
			x ++;
		}
		if (x >= 2) contigComponents[v[1]] = v[0];
	}
	return outputOverlap(filename, proteinIndex, 0, range, 0.0, 0, 0, false, line_length);
	checkComps = false;
}

bool OverlapOut::outputOrbcalOverlapComp(const char* filename, const char* componentFile, int proteinIndex, float range, short precCharge, float tolerance, int line_length) {
	BufferedLineReader comps;
	comps.Load(componentFile);
	contigComponents.clear();
	checkComps = true;
	char* next;
	char* result;
	vector<int> v(10);
	for (int i = 1; i < comps.size(); i++) {
		next = comps.getline(i);
		if (strlen( next ) == 0) continue;
		result = strtok( next, "," );
		int x = 0;
		while( result != NULL ) {
			v[x] = getInt(result);
			result = strtok( NULL, "," );
			x ++;
		}
		if (x >= 2) contigComponents[v[1]] = v[0];
	}
	return outputOverlap(filename, proteinIndex, 1, range, tolerance, precCharge, 0, false, line_length);
	checkComps = false;
}

bool OverlapOut::outputOrbcalOverlap(const char* filename, int proteinIndex, short precCharge, float range, float tolerance, int line_length) {return outputOverlap(filename, proteinIndex, 1, range, tolerance, precCharge, 0, false, line_length);}

bool OverlapOut::outputOrbcalContigFiltOverlap(const char* filename, int proteinIndex, short precCharge, float range, float tolerance, int match, bool checkConsPeaks, int line_length) {return outputOverlap(filename, proteinIndex, 1, range, tolerance, precCharge, match, checkConsPeaks, line_length);}

void OverlapOut::pruneContigs(int proteinIndex, float range, list<TwoValues<int> >& index_order, TwoValues<int>& curPair)
{
	Spectrum PepSpec = fasta.getMassesSpec(proteinIndex);
	proteinToContig.clear();
	proteinToContig.resize(PepSpec.size());
	proteinToConSt.clear();
	proteinToConSt.resize(PepSpec.size());
	
	//cout << proteinIndex << ", fasta size:" << fasta.size() << endl;
	
	for (int i = 0; i < overlaps.size(); i ++) {
		//if (checkComps && (contigComponents.count(i) == 0 || contigComponents[i] != 13)) continue;
		if (overlaps[i].size() == 0) continue;
		// only use contigs for this protein
		if (prot_match[i][0] != proteinIndex || contigs[i].size() == 0) continue;
		//reverse any specified contigs
		if (prot_match[i][2] == 1 && reversed.count(i) == 0) {
			contigs[i].reverse(0.0 - AAJumps::massH2O, 0);
			reversed.insert(i);
		}
		if (i == 46) {
		  cout << "Contig " << i << ": (size = " << contigs[i].size() << ", parent mass = " << contigs[i].parentMass << ")\n";
		  for (int j = 0; j < overlaps[i].size(); j++) {
		    int peakIdx = (int)(overlaps[i][j][0]+0.01);
		    cout << "idx = " << j << ", peakIdx = " << peakIdx << ", mass = " << contigs[i][peakIdx][0] << ", protIdx = " << overlaps[i][j][1] << "\n";
		  }
		  fflush(stdout);
		}
		//associate every contig index with its first matching index in the peptide
		contprot[i] = proteinIndex + 1;
		//cout << i << ":" <<  overlaps[i][0][1] << ":" << proteinToConSt.size() << ":" << PepSpec.size() << endl;
		proteinToConSt[(int)(overlaps[i][0][1]+0.01)].push_back(i);
		curPair.set(i, (int)(overlaps[i][0][1]+0.01));
		rememberContig(i, range, &PepSpec);
		index_order.push_back(curPair);
		char* peptide = fasta.getSequence(proteinIndex);
		char sequence[strlen(peptide)];
		strcpy(sequence, peptide);
		sequence[(int)(overlaps[i][overlaps[i].size()-1][1]+0.01)+1] = 0;
		contigSeq[i] = &sequence[(int)(overlaps[i][0][1] + 0.01)];
	}
	// sort contigs by their first matching index in the peptide
	index_order.sort(SortOverlap());
}

void OverlapOut::pruneOrbcals(int proteinIndex, short precCharge, float tolerance, list<TwoValues<int> >& orbcal_order, TwoValues<int>& curPair)
{
	Spectrum PepSpec = fasta.getMassesSpec(proteinIndex);
	Spectrum nextSpec;
	char* peptide = fasta.getSequence(proteinIndex);
	char* match;
	int orbcal_match;
	int prev;
	proteinToOrbcal.clear();
	proteinToOrbcal.resize(PepSpec.size());
	proteinToOrbSt.clear();
	proteinToOrbSt.resize(PepSpec.size());
	
	for (unsigned int i = 0; i < orbcal_pep.size() && i < orbcals.size(); i++) {
		
		const char* next_contig = orbcal_pep.getline(i);
		if (next_contig[0] == '\0' || strlen(next_contig) == 0) continue;
		
		nextSpec = orbcals[i];
		if (nextSpec.size() == 0) continue;
		if (nextSpec.parentCharge < precCharge) continue;
		char just_letters[strlen(next_contig)+1];
		getPepSeq(next_contig, just_letters);
		orbcal_match = -1; prev = -2;
		
		while (1) {
			//find all matches of orbcal in target peptide
			prev = orbcal_match;
			match = strstr(peptide, just_letters);
			if (match == NULL || &match[0] == &peptide[0]) break;
			orbcal_match = &match[0] - &peptide[0];
			if (orbcal_match != prev) {
				orbprot[i] = proteinIndex + 1;
				proteinToOrbSt[orbcal_match].push_back(i);
				curPair.set(i, orbcal_match);
				orbcal_order.push_back(curPair);
				rememberOrbcal(i, orbcal_match, tolerance);
				orbcalSeq[i] = next_contig;
			} else break;
		}
	}
	// sort orbcals by their first matching index in the peptide
	orbcal_order.sort(SortOverlap());
}

void OverlapOut::outputStats(char* outfile, list<Results_OCC>& alignment, short precCharge, float range, float tolerance, float parentMassTol, int startMinNumMatchedPeaks, int endMinNumMatchedPeaks, int stepMinNumMatchedPeaks, float startMinRatio, float endMinRatio, float stepMinRatio){
	
	orbcalContShifts.clear();
	orbcalContShifts.resize(orbcals.size());
	vector<float> curVec(5);
	list<TwoValues<int> > contig_order;
	list<TwoValues<int> > orbcal_order;
	list<vector<float> > ratioPeak;
	list<vector<float> >::iterator it;
	vector<float> curDat(3);
	TwoValues<int> curPair;
	list<TwoValues<int> >::iterator pivot_order;
	list<int>::iterator pivot_spec;
	list<Results_OCC>::iterator pivot_res;
	int contig_i;
	int orbcal_i;
	bool correct;
	bool thresh;
	bool haveRes = false;
	int tp, fn, fp, t;
	float minR;
	char* next_contig;
	contprot.clear(); contprot.resize(contigs.size());
	orbprot.clear(); orbprot.resize(orbcals.size());
	
	for (int i = 0; i < fasta.size(); i++) {
		contig_order.clear();
		orbcal_order.clear();
		pruneContigs(i, range, contig_order, curPair);
		pruneOrbcals(i, precCharge, tolerance, orbcal_order, curPair);
		for (int j = 0; j < proteinToConSt.size(); j ++) {
			if (proteinToConSt[j].size() == 0) continue;
			for (pivot_spec = proteinToConSt[j].begin(); pivot_spec != proteinToConSt[j].end(); pivot_spec ++) {
				contig_i = *pivot_spec;
				getCorrectShifts(contig_i, j, i, proteinToOrbSt, proteinToOrbcal, range, parentMassTol, true);
			}
		}
		for (int j = 0; j < proteinToOrbSt.size(); j ++) {
			if (proteinToOrbSt[j].size() == 0) continue;
			for (pivot_spec = proteinToOrbSt[j].begin(); pivot_spec != proteinToOrbSt[j].end(); pivot_spec ++) {
				orbcal_i = *pivot_spec;
				getCorrectShifts(orbcal_i, j, i, proteinToConSt, proteinToOrbcal, tolerance, parentMassTol, false);
			}
		}
	}
	
	FILE* output = fopen( outfile,  "w");
	//FILE* output2 = fopen( "debug_align_fp.txt",  "w");
	//FILE* output3 = fopen( "debug_align_tp.txt",  "w");
	fprintf(output, "MinNumMatchedPeaks%sminRatio%sTP%sFN%sFP%sNumAboveThreshhold%sTP/(TP+FP)%sTP/(TP+FN)\n", DELIM, DELIM, DELIM, DELIM, DELIM, DELIM, DELIM);
	int uneq = 0;
	int tot = 0;
	for (int minMatchedPeaks = startMinNumMatchedPeaks; minMatchedPeaks <= endMinNumMatchedPeaks; minMatchedPeaks += stepMinNumMatchedPeaks) {
		for (float minRatio = startMinRatio; minRatio < (endMinRatio + stepMinRatio); minRatio += stepMinRatio) {
			tp = 0; fn = 0; fp = 0; t = 0;
			for (pivot_res = alignment.begin(); pivot_res != alignment.end(); pivot_res ++) {
				Results_OCC result = *pivot_res;
				if (result.spec1 >= orbcal_pep.size()) continue;
				if (orbcals[result.spec1].parentCharge < precCharge) continue;
				next_contig = orbcal_pep.getline(result.spec1);
				
				if (strlen(next_contig) == 0) continue;
				minR = min(result.score1, result.score2);
				
				correct = false;
				thresh = false;
				if (orbcalContShifts[result.spec1].count(result.spec2) > 0) {
					float shift = orbcalContShifts[result.spec1][result.spec2];
					if (isEqual(shift, result.shift2, parentMassTol) || isEqual(shift, result.shift1, parentMassTol)) {
						correct = true;
					}
				}
				if (result.matchedPeaks >= minMatchedPeaks && minR >= minRatio) {
					thresh = true;
				}
				if (!haveRes && result.matchedPeaks >= startMinNumMatchedPeaks && minR >= startMinRatio) {
					curDat[0] = minR;
					curDat[1] = result.matchedPeaks;
					curDat[2] = (correct) ? 1.0 : 0.0;
					ratioPeak.push_back(curDat);
				}
				if (thresh) t ++;
				
				/*
				if (!haveRes && thresh && ! correct && result.matchedPeaks >= 7 && minR >= 0.2) {
					fprintf(output2, "orbcal: %d contig : %d numMatchedPeaks: %d minRatio: %.3f o_prot: ", result.spec1, result.spec2, result.matchedPeaks, minR);
					fprintf(output2, (orbprot[result.spec1] > 0) ? parseInt(orbprot[result.spec1] - 1).c_str() : "None");
					fprintf(output2, " c_prot: ");
					fprintf(output2, (contprot[result.spec2] > 0) ? parseInt(contprot[result.spec2] - 1).c_str() : "None");
					fprintf(output2, " shift1: %.1f shift2: %.1f", result.shift1, result.shift2);
					fprintf(output2, "\norbcalSeq: %s\ncontigSeq: %s\nCorrect shift: ", orbcalSeq[result.spec1].c_str(), contigSeq[result.spec2].c_str());
					if (orbcalContShifts[result.spec1].count(result.spec2) > 0) {
						float shift = orbcalContShifts[result.spec1][result.spec2];
						if (prot_match[result.spec2][2] == 1) fprintf(output2, "**");
						fprintf(output2, "%.1f", shift);
					} else fprintf(output2, "None");
					tot ++;
					if (orbprot[result.spec1] != contprot[result.spec2]) uneq ++;
					fprintf(output2, "\n");
				}
				
				if (!haveRes && thresh && correct && result.matchedPeaks >= 6){// && minR >= 0.225) {
					fprintf(output3, "orbcal: %d contig: %d numMatchedPeaks: %d minRatio: %.3f o_prot1: ", result.spec1, result.spec2, result.matchedPeaks, minR);
					fprintf(output3, (orbprot[result.spec1] > 0) ? parseInt(orbprot[result.spec1] - 1).c_str() : "None");
					fprintf(output3, " c_prot2: ");
					fprintf(output3, (contprot[result.spec2] > 0) ? parseInt(contprot[result.spec2] - 1).c_str() : "None");
					fprintf(output3, " shift1: %.1f shift2: %.1f\n",result.shift1, result.shift2);
				}
				*/
				
				if (correct && thresh) tp ++;
				else if (correct && ! thresh) fn ++;
				else if (! correct && thresh) fp ++;
				else continue;
			}
			haveRes = true;
			fprintf(output, "%d%s%.3f%s%d%s%d%s%d%s%d%s%.3f%s%.3f\n", minMatchedPeaks, DELIM, minRatio, DELIM, tp, DELIM, fn, DELIM, fp, DELIM, t, DELIM, ((float)tp)/((float)tp + (float)fp), DELIM, ((float)tp)/((float)tp + (float)fn));
		}
	}
	fprintf(output, "\n\nMinRatio%sNumPeaks%sCorrect?\n", DELIM, DELIM);
	for (it = ratioPeak.begin(); it != ratioPeak.end(); it ++) fprintf(output, "%.3f%s%.0f%s%.0f\n", (*it)[0], DELIM, (*it)[1], DELIM, (*it)[2]);
	//cout << uneq << " of " << tot << " do not match proteins.\n";
	fclose(output);
	//fclose(output2);
	//fclose(output3);
}

void OverlapOut::outputStats(char* outfile, char* inAlignmentFile, short precCharge, float range, float tolerance, float parentMassTol, int startMinNumMatchedPeaks, int endMinNumMatchedPeaks, int stepMinNumMatchedPeaks, float startMinRatio, float endMinRatio, float stepMinRatio){
	list<Results_OCC> alignment;
	if (! Load_resultsOCC(inAlignmentFile, alignment) ) cerr << "ERROR: Cannot load " << inAlignmentFile << ".\n";
	outputStats(outfile, alignment, precCharge, range, tolerance, parentMassTol, startMinNumMatchedPeaks, endMinNumMatchedPeaks, stepMinNumMatchedPeaks, startMinRatio, endMinRatio, stepMinRatio);
}


bool OverlapOut::outputOverlap(const char* filename, int proteinIndex, short showOrbcal, float range, float tolerance, short precCharge, int filterOrbcal, bool checkConsPeaks, int line_length)
{
	int LINE_LENGTH = line_length;
	Spectrum PepSpec = fasta.getMassesSpec(proteinIndex);
	char* peptide = fasta.getSequence(proteinIndex);
	FILE* output = fopen( filename,  "w");
	Spectrum nextSpec;
	list<TwoValues<int> > index_order;
	list<TwoValues<int> > orbcal_order;
	list<TwoValues<int> >::iterator pivot;
	list<TwoValues<int> >::iterator orbcal_pivot;
	list<TwoValues<int> > backwash_p;
	list<TwoValues<int> > backwash_n;
	list<TwoValues<int> >::iterator backpivot;
	list<vector<int> > orbcal_backwash_p;
	list<vector<int> > orbcal_backwash_n;
	list<vector<int> >::iterator orbcal_backpivot;
	list<int>::iterator pivot_spec;
	TwoValues<int> curPair;
	vector<int> curVec(3);
	vector<float> orb_masses;
	int contig_index;
	int backwash_index;
	int orbcal_backwash_index;
	contprot.clear();
	orbprot.clear();
	contprot.resize(contigs.size());
	orbprot.resize(orbcals.size());
	orbcalsFiltered.clear();
	orbcalsFiltered.resize(orbcals.size());
	contigsFiltered.clear();
	contigsFiltered.resize(contigs.size());
	orbcalContShifts.clear();
	orbcalContShifts.resize(orbcals.size());
	Bpeaks = 0;
	Ypeaks = 0;
	Totpeaks = 0;
	if (showOrbcal == 1) {
		if (usesOrbcal == 0 || precCharge == 0) {
			cerr << "ERROR: Must pass constructor orbcal sequence and specset files to display orbcal overlap. PrecCharge must also be > 0\n";
			return false;
		}
		pruneContigs(proteinIndex, range, index_order, curPair);
		pruneOrbcals(proteinIndex, precCharge, tolerance, orbcal_order, curPair);
	}
	else pruneContigs(proteinIndex, range, index_order, curPair);
	
	float percB = ((float)Bpeaks)/((float)Totpeaks);
	float percY = ((float)Ypeaks)/((float)Totpeaks);
	//cout << "Ratio B peaks present in orbcals: " << percB << endl;
	//cout << "Ratio Y peaks present in orbcals: " << percY << endl;
	if (filterOrbcal > 0) {
		for (int j = 0; j < proteinToConSt.size(); j ++) {
			if (proteinToConSt[j].size() == 0) continue;
			for (pivot_spec = proteinToConSt[j].begin(); pivot_spec != proteinToConSt[j].end(); pivot_spec ++) {
				int contig_i = *pivot_spec;
				getCorrectShifts(contig_i, j, proteinIndex, proteinToOrbSt, proteinToOrbcal, range, range, true);
			}
		}
		for (int j = 0; j < proteinToOrbSt.size(); j ++) {
			if (proteinToOrbSt[j].size() == 0) continue;
			for (pivot_spec = proteinToOrbSt[j].begin(); pivot_spec != proteinToOrbSt[j].end(); pivot_spec ++) {
				int orbcal_i = *pivot_spec;
				getCorrectShifts(orbcal_i, j, proteinIndex, proteinToConSt, proteinToOrbcal, tolerance, range, false);
			}
		}
		for (orbcal_pivot = orbcal_order.begin(); orbcal_pivot != orbcal_order.end(); orbcal_pivot ++) {
			isFilteredOrbcal((*orbcal_pivot)[0], (*orbcal_pivot)[1], tolerance, filterOrbcal, range, checkConsPeaks);
		}
		for (pivot = index_order.begin(); pivot != index_order.end(); pivot ++) {
			isFilteredContig(((*pivot)[0]), filterOrbcal, range, checkConsPeaks);
		}
	}
	
	if (checkComps) fprintf(output, "%sComponent\n", DELIM);
	else fprintf(output, "%sParent Charge\n", DELIM);
	pivot = index_order.begin();
	if (showOrbcal == 1) orbcal_pivot = orbcal_order.begin();
	for (int l = 0; l <= ((PepSpec.size() -1)/LINE_LENGTH); l ++) {
		fprintf(output, "%s", DELIM);
		for (int i = (l*LINE_LENGTH); (i < (l*LINE_LENGTH)+LINE_LENGTH) && (peptide[i] != 0); i++) {
			// write peptide segment
			fprintf(output, "%s%s", DELIM, DELIM);
			fputc(peptide[i], output);
		}
		fprintf(output, "%s%s\n", DELIM, DELIM);
		backwash_n.clear();
		for (backpivot = backwash_p.begin(); backpivot != backwash_p.end(); backpivot ++) {
			// write left-over contigs from previous line(s) and save their last written index
			backwash_index = outputContig(output, peptide, &PepSpec, (*backpivot)[0], (l*LINE_LENGTH), (*backpivot)[1], (l*LINE_LENGTH)+LINE_LENGTH, range);
			if (backwash_index >= 0) {
				// remember contig and last written index for next line
				curPair.set((*backpivot)[0], backwash_index);
				backwash_n.push_back(curPair);
			}
			fprintf(output, "\n");
		}
		
		while(pivot != index_order.end() && (*pivot)[1] < (l*LINE_LENGTH)+LINE_LENGTH) {
			if (filterOrbcal > 0 && ! contigsFiltered[(*pivot)[0]]) {
				pivot++;
				continue;
			}
			// write contigs in the order of their first matching protein index
			backwash_index = outputContig(output, peptide, &PepSpec, ((*pivot)[0]), (l*LINE_LENGTH), 0, (l*LINE_LENGTH)+LINE_LENGTH, range);
			
			if (backwash_index >= 0) {
				// remember contig and last written index for next line
				curPair.set(((*pivot)[0]), backwash_index);
				backwash_n.push_back(curPair);
			}
			fprintf(output, "\n");
			pivot++; // go to next contig
		}
		
		if (showOrbcal == 1) {
			orbcal_backwash_n.clear();
			for (orbcal_backpivot = orbcal_backwash_p.begin(); orbcal_backpivot != orbcal_backwash_p.end(); orbcal_backpivot ++) {
				// write left-over contigs from previous line(s) and save their last written index
				orbcal_backwash_index = outputOrbcalContig(output, peptide, (*orbcal_backpivot)[0], (l*LINE_LENGTH), (*orbcal_backpivot)[1], (*orbcal_backpivot)[2], (l*LINE_LENGTH)+LINE_LENGTH, tolerance);
				if (orbcal_backwash_index >= 0) {
					// remember contig and last written index for next line
					curVec[0] = (*orbcal_backpivot)[0];
					curVec[1] = (l*LINE_LENGTH)+LINE_LENGTH;
					curVec[2] = orbcal_backwash_index;
					orbcal_backwash_n.push_back(curVec);
				}
				fprintf(output, "\n");
			}
			
			while(orbcal_pivot != orbcal_order.end() && (*orbcal_pivot)[1] < (l*LINE_LENGTH)+LINE_LENGTH) {
				// write contigs in the order of their first matching protein index
				if (filterOrbcal > 0 && ! orbcalsFiltered[(*orbcal_pivot)[0]]) {
					orbcal_pivot++;
					continue;
				}
				orbcal_backwash_index = outputOrbcalContig(output, peptide, (*orbcal_pivot)[0], (l*LINE_LENGTH), (*orbcal_pivot)[1], 0, (l*LINE_LENGTH)+LINE_LENGTH, tolerance);
				if (orbcal_backwash_index >= 0) {
					// remember contig and last written index for next line
					curVec[0] = (*orbcal_pivot)[0];
					curVec[1] = (l*LINE_LENGTH)+LINE_LENGTH;
					curVec[2] = orbcal_backwash_index;
					orbcal_backwash_n.push_back(curVec);
				}
				fprintf(output, "\n");
				orbcal_pivot++; // go to next contig
			}
		}
		// switch overflow queues
		backwash_p = backwash_n;
		
		if (showOrbcal == 1) {
			orbcal_backwash_p = orbcal_backwash_n;
		}
		fprintf(output, "\n");
	}
	fclose(output);
	return true;
}

int OverlapOut::outputContig(FILE* output, char* peptide, Spectrum* PepSpec, int contig_ident, int start_pep_index, int prev_c, int end_pep_index, float range) {
	fprintf(output, "%d%s", contig_ident, DELIM);
	Spectrum contig = contigs[contig_ident];
	if (checkComps) {
		if (contigComponents.count(contig_ident) > 0) fprintf(output, "%d", contigComponents[contig_ident]);
		else fprintf(output, "-1");
	} else fprintf(output, "%d", contig.parentCharge);
	Spectrum overlap = overlaps[contig_ident];
	int start = start_pep_index;
	float pep_mass;
	float contig_mass;
	int backwash_index = -1;
	int i;
	
	for (i = start_pep_index; i < overlap[0][1]; i++) fprintf(output, "%s%s", DELIM, DELIM);
	
	if (i == overlap[prev_c][1]) fprintf(output, "%s|%s", DELIM, DELIM); // write very first break
	// loop until end of contig or end of line is reached
	for (; (i <= end_pep_index) && (i <= overlap[overlap.size()-1][1]) && peptide[i] != 0; i++) {
		if ((overlap[prev_c + 1][1] == i) || ((i == end_pep_index) && (i < overlap[overlap.size()-1][1]))) {
			pep_mass = (*PepSpec)[(int)(overlap[prev_c + 1][1] + 0.01)][0] - (*PepSpec)[(int)(overlap[prev_c][1] + 0.01)][0];
			contig_mass = contig[(int)(overlap[prev_c + 1][0] + 0.01)][0] - contig[(int)(overlap[prev_c][0]+0.01)][0];
			if ((pep_mass < contig_mass-range) || (pep_mass > contig_mass+range)) {
				// write mass difference if it is outside of our limits
				if (overlap[prev_c][1] >= start_pep_index) {
					if (contig_mass > pep_mass) {
						fprintf(output, "+");
					}
					fprintf(output, "%.1f", contig_mass - pep_mass);
				}
				for (int k = (int)(overlap[prev_c][1]+0.01)+1; k < i; k++) {
					if (k < start_pep_index) {continue;}
					// fill in unknown region
					fprintf(output, "%s%s", DELIM, DELIM);
				}
				if (overlap[prev_c + 1][1] == i) {fprintf(output, "%s|%s", DELIM, DELIM);} // write ending contig segment break
			} else {
				//write AA's for matched contig segment 
				if (overlap[prev_c][1] >= start_pep_index) {
					fputc(peptide[(int)(overlap[prev_c][1]+0.01)], output);
				}
				for (int k = (int)(overlap[prev_c][1]+0.01)+1; k < i; k++) {
					if (k < start_pep_index) {continue;}
					fprintf(output, "%s.%s", DELIM, DELIM);
					fputc(peptide[k], output);
				}
				// write ending contig segment break
				if (overlap[prev_c + 1][1] == i) {
					fprintf(output, "%s|%s", DELIM, DELIM);
				} else fprintf(output, "%s.%s", DELIM, DELIM);
			}
			prev_c ++; // go to next contig break
		}
	}
	
	if ((i-1 == end_pep_index) && (i-1 < overlap[overlap.size()-1][1])) {
		// return index of last AA written if more of the contig needs to be written
		backwash_index = prev_c - 1;
	}
	return backwash_index;
}

int OverlapOut::outputOrbcalContig(FILE* output, char* peptide, int orbcal_ident, int start_pep_index, int start_overlap, int prev_c, int end_pep_index, float tolerance) {
	Spectrum orbcal = orbcals[orbcal_ident];
	fputc('*', output);
	fprintf(output, "%d%s%d", orbcal_ident, DELIM, orbcal.parentCharge);
	vector<float> perfect_orbcal;
	vector<float> no_offset;
	int start_contig_index = prev_c;
	int backwash_index = -1;
	getMassesCummulative(orbcal_pep.getline(orbcal_ident), perfect_orbcal, 0.0);
	getMassesCummulativeNoShift(orbcal_pep.getline(orbcal_ident), no_offset);
	float mass;
	float shiftmass;
	int i;
	
	/*if (orbcal_ident == 6302 && prev_c == 0) {
		for (int jk = 0; jk < perfect_orbcal.size(); jk ++) {
			float yp = perfect_orbcal[perfect_orbcal.size()-1] - perfect_orbcal[jk] + AAJumps::massH2O + 1.007;
			cout << perfect_orbcal[jk] + 1.007 << " ; " << orbcal[orbcal.findClosest(perfect_orbcal[jk] + 1.007)][0] << " ; " << yp << " ; " << orbcal[orbcal.findClosest(yp)][0] << "\n";
		}
	}*/
	
	for (i = start_pep_index; i < start_overlap; i++) fprintf(output, "%s%s", DELIM, DELIM);
	// loop until end of contig or end of line is reached
	outputOrbcalBreak(output, &orbcal, perfect_orbcal[prev_c], perfect_orbcal[perfect_orbcal.size()-1] - perfect_orbcal[prev_c], tolerance);
	prev_c ++;
	
	for (; (i < end_pep_index) && (prev_c < perfect_orbcal.size()) && peptide[i] != 0; i++) {
		fputc(peptide[i], output);
		mass = no_offset[prev_c] - no_offset[prev_c - 1];
		shiftmass = perfect_orbcal[prev_c] - perfect_orbcal[prev_c - 1];
		
		if (mass <  shiftmass - 0.001 || mass > shiftmass + 0.001) {
			if (shiftmass > mass) {
				fputc('+', output);
			}

			fprintf(output, "%.0f", shiftmass - mass);
		}
		outputOrbcalBreak(output, &orbcal, perfect_orbcal[prev_c], perfect_orbcal[perfect_orbcal.size()-1] - perfect_orbcal[prev_c], tolerance);
		prev_c ++;
	}
	if (prev_c < perfect_orbcal.size()) {
		// return index of last AA written if more of the contig needs to be written
		backwash_index = prev_c - 1;
	}
	return backwash_index;
}

void OverlapOut::getCorrectShifts(int index, int start, int pepIndex, vector<list<int> >& start_order_other, vector<set<int> >& proteinToOrbcal, float tolrange, float pmTol, bool contig) {
	char* peptide = fasta.getSequence(pepIndex);
	vector<float> perfect_spec;
	Spectrum overlap;
	Spectrum contig_spec;
	if (contig) {
		overlap = overlaps[index];
		contig_spec = contigs[index];
	} else getMassesCummulative(orbcal_pep.getline(index), perfect_spec, 0.0);
	
	int end = (contig) ? (int)(overlap[overlap.size()-1][1] +0.01) : start + perfect_spec.size() - 1;
	list<int>::iterator start_other_pivot;
	map<int, float> shiftScore;
	map<int, float>::iterator shiftScoreIt;
	bool foundShift;
	float shift;
	int intShift;
	float peak;
	float score;
	float total_intensity;
	float intensity;
	float maxScore;
	int orb_index;
	int intPmTol = (int)((pmTol/InputParams::Resolution)+0.01);
	
	for (int i = start; i <= end && peptide[i] != 0; i++) {
		//while (contig && (int)overlap[next][1] < i && next < overlap.size()) next ++;
		if (start_order_other[i].size() == 0) continue;
		
		for (start_other_pivot = start_order_other[i].begin(); start_other_pivot != start_order_other[i].end(); start_other_pivot++) {
			
			shiftScore.clear();
			if (!contig) {
				overlap = overlaps[*start_other_pivot];
				contig_spec = contigs[*start_other_pivot];
			} else getMassesCummulative(orbcal_pep.getline(*start_other_pivot), perfect_spec, 0.0);
			
			total_intensity = 0.0;
			for (int x = 0; x < contig_spec.size(); x ++) total_intensity += contig_spec[x][1];
			
			for (int pep_next = i; pep_next <= end && peptide[pep_next] != 0; pep_next ++) {
				if (contig && proteinToOrbcal[pep_next].count(*start_other_pivot) == 0) continue;
				if (!contig && proteinToOrbcal[pep_next].count(index) == 0) continue;
				orb_index = (contig) ? pep_next - i : pep_next - start;
				
				for (int next = 0; next < overlap.size() && (int)(overlap[next][1]+0.01) <= pep_next; next ++) {
					if ((int)(overlap[next][1]+0.01) != pep_next) continue;
					intShift = (int)(((perfect_spec[orb_index] - contig_spec[(int)(overlap[next][0]+0.01)][0])/InputParams::Resolution)+0.01);
					intensity = contig_spec[(int)(overlap[next][0]+0.01)][1];
					for (int idx = intShift - intPmTol; idx <= intShift + intPmTol; idx ++) {
						if (shiftScore.count(idx) == 0) shiftScore[idx] = (intensity/total_intensity) - (0.01*(float)abs(intShift-idx));
						else shiftScore[idx] += (intensity/total_intensity) - (0.01*(float)abs(intShift-idx));
					}
				}
			}
			maxScore = 0.0;
			for (shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt ++) {
				if (shiftScoreIt->second > maxScore) {
					maxScore = shiftScoreIt->second;
					shift = ((float)shiftScoreIt->first)*InputParams::Resolution;
				}
			}
			if (maxScore > 0.0) {
				if (contig) {
					orbcalContShifts[*start_other_pivot][index] = shift;
				} else {
					orbcalContShifts[index][*start_other_pivot] = shift;
				}
			}
		}
	}
}

TwoValues<float> OverlapOut::getContigOverlapScore(int idx1, int idx2, float shiftFor, float shiftRev, float pmTol, bool debug) {
	Spectrum contig;
	if (reversed.count(idx1) > 0) contigs[idx1].reverse(0.0, &contig);
	else contig = contigs[idx1];
	Spectrum f_contig_o;
	if (reversed.count(idx2) > 0) contigs[idx2].reverse(0.0, &f_contig_o);
	else f_contig_o = contigs[idx2];
	return getContigOverlapScores(contig, f_contig_o, shiftFor, shiftRev, pmTol, 0.0, debug);
}

void OverlapOut::getCorrectConnShifts(int index, int start, int pepIndex, vector<list<int> >& start_order_orbcal, vector<set<int> >& proteinToOrbcal, float tolerance, float pmTol) {
	char* peptide = fasta.getSequence(pepIndex);
	vector<float> perfect_spec;
	
	Spectrum orbcal1 = orbcals[index];
	getMassesCummulative(orbcal_pep.getline(index), perfect_spec, 0.0);
	
	int end = start + perfect_spec.size() - 1;
	list<int>::iterator start_other_pivot;
	map<int, float> shiftScore;
	map<int, float>::iterator shiftScoreIt;
	bool foundShift;
	float shift;
	float score;
	float maxScore;
	int intPmTol = (int)((pmTol/InputParams::Resolution)+0.01);
	
	for (int i = start; i <= end && peptide[i] != 0; i++) {
		//while (contig && (int)overlap[next][1] < i && next < overlap.size()) next ++;
		if (start_order_orbcal[i].size() == 0) continue;
		
		for (start_other_pivot = start_order_orbcal[i].begin(); start_other_pivot != start_order_orbcal[i].end(); start_other_pivot++) {
			
			shiftScore.clear();
			vector<float> perfect_spec2;
			Spectrum orbcal2 = orbcals[*start_other_pivot];
			getMassesCummulative(orbcal_pep.getline(*start_other_pivot), perfect_spec2, 0.0);
			
			for (int pep_next = i; pep_next <= end && peptide[pep_next] != 0; pep_next ++) {
				if (proteinToOrbcal[pep_next].count(*start_other_pivot) == 0 || proteinToOrbcal[pep_next].count(index) == 0) continue;
				
				int orb1_index = orbcal1.findClosest(perfect_spec[pep_next - start]);
				int orb2_index = orbcal2.findClosest(perfect_spec[pep_next - i]);
				float intensity = orbcal1[orb1_index][1] + orbcal2[orb2_index][1];
				int intShift = (int)(((orbcal1[orb1_index][0] - orbcal2[orb2_index][0])/InputParams::Resolution)+0.01);
				for (int idx = intShift - intPmTol; idx <= intShift + intPmTol; idx ++) {
					if (shiftScore.count(idx) == 0) shiftScore[idx] = intensity - (0.01*(float)abs(intShift-idx));
					else shiftScore[idx] += intensity - (0.01*(float)abs(intShift-idx));
				}
			}
			maxScore = 0.0;
			for (shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt ++) {
				if (shiftScoreIt->second > maxScore) {
					maxScore = shiftScoreIt->second;
					shift = ((float)shiftScoreIt->first)*InputParams::Resolution;
				}
			}
			if (maxScore > 0.0) {
				orbcalOrbcalShifts[index][*start_other_pivot] = shift;
			}
		}
	}
}

void OverlapOut::getContigDistance(int index, vector<list<int> >& start_order_other, float pmTol) {
	unsigned int pepIndex = prot_match[index][0];
	char* peptide = fasta.getSequence(pepIndex);
	Spectrum masses = fasta.getMassesSpec(pepIndex);
	
	Spectrum overlap = overlaps[index];
	Spectrum contig_spec = contigs[index];
	
	int end = (int)(overlap[overlap.size()-1][1]+0.01);
	list<int>::iterator start_other_pivot;
	map<int, float> shiftScore;
	map<int, float>::iterator shiftScoreIt;
	float best_shift1;
	float best_shift2;
	float shift;
	float intensity;
	char sequence[strlen(peptide)];
	strcpy(sequence, peptide);
	int intShift;
	int intPmTol = (int)((pmTol/InputParams::Resolution)+0.01);
	
	for (int i = 0; i < overlap.size(); i++) {
		int cIdx = (int)(overlap[i][0] + 0.01);
		int pIdx = (int)(overlap[i][1] + 0.01);
		shift = masses[pIdx][0] - contig_spec[cIdx][0];
		//if (index == 98) cout << cIdx << " : " << contig_spec[cIdx][0] << ", " << pIdx << " : " << masses[pIdx][0] << " = " << shift << "\n";
		intShift = (int)((shift/InputParams::Resolution) + 0.01);
		intensity = contig_spec[cIdx][1];
		
		for (int idx = intShift - intPmTol; idx <= intShift + intPmTol; idx ++) {
			if (shiftScore.count(idx) == 0) shiftScore[idx] = intensity - (0.001*(float)abs(intShift-idx));
			else shiftScore[idx] += intensity - (0.001*(float)abs(intShift-idx));
		}
	}
	float maxScore = 0.0;
	for (shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt ++) {
		if (shiftScoreIt->second > maxScore) {
			maxScore = shiftScoreIt->second;
			best_shift1 = ((float)shiftScoreIt->first)*InputParams::Resolution;
		}
	}
	
	for (int i = (int)overlap[0][1]; peptide[i] != 0; i++) {
		
		if (start_order_other[i].size() == 0) continue;
		
		for (start_other_pivot = start_order_other[i].begin(); start_other_pivot != start_order_other[i].end(); start_other_pivot++) {
			shiftScore.clear();
			Spectrum other_overlap = overlaps[*start_other_pivot];
			Spectrum other_contig = contigs[*start_other_pivot];
			
			for (int j = 0; j < other_overlap.size(); j++) {
				int cIdx = (int)(other_overlap[j][0] + 0.01);
				int pIdx = (int)(other_overlap[j][1] + 0.01);
				shift = masses[pIdx][0] - other_contig[cIdx][0];
				intShift = (int)((shift/InputParams::Resolution)+0.01);
				intensity = other_contig[cIdx][1];
				//if (index == 98 && *start_other_pivot == 60) cout << cIdx << " : " << other_contig[cIdx][0] << ", " << pIdx << " : " << masses[pIdx][0] << " = " << shift << " = " << shift - best_shift1 << "\n";
				for (int idx = intShift - intPmTol; idx <= intShift + intPmTol; idx ++) {
					if (shiftScore.count(idx) == 0) shiftScore[idx] = intensity - (0.001*(float)abs(intShift-idx));
					else shiftScore[idx] += intensity - (0.001*(float)abs(intShift-idx));
				}
			}
			maxScore = 0.0;
			for (shiftScoreIt = shiftScore.begin(); shiftScoreIt != shiftScore.end(); shiftScoreIt ++) {
				if (shiftScoreIt->second > maxScore) {
					maxScore = shiftScoreIt->second;
					best_shift2 = ((float)shiftScoreIt->first)*InputParams::Resolution;
				}
			}
			if (index == 7 && *start_other_pivot == 277) cout << "shift: " << best_shift2 - best_shift1 << "\n";
			contContShifts[index][*start_other_pivot] = best_shift2 - best_shift1;
		}
	}
	
}

float OverlapOut::getMass(char *aa) {
	int m = 0;
	while (m<20 and AAletters[m]!=aa[0]) m++;
	if(m<20) {
		return AAmasses[m];
	} else return 0.0;
}

void OverlapOut::outputOrbcalBreak(FILE* output, Spectrum* orbcal, float bpeak_, float ypeak_, float tolerance) {
	float bpeak = bpeak_;
	float ypeak = ypeak_ + AAJumps::massH2O;
	int bloc = (*orbcal).findClosest(bpeak);
	int yloc = (*orbcal).findClosest(ypeak);
	bool isy = false;
	bool isb = false;
	isb = isEqual((*orbcal)[bloc][0], bpeak, tolerance);
	isy = isEqual((*orbcal)[yloc][0], ypeak, tolerance);
	if (isb && !isy ) fprintf(output, "%s\\%s", DELIM, DELIM);
	else if (!isb && isy ) fprintf(output, "%s/%s", DELIM, DELIM);
	else if (isb && isy ) fprintf(output, "%s|%s", DELIM, DELIM);
	else fprintf(output, "%s.%s", DELIM, DELIM);
}

void OverlapOut::rememberOrbcal(int index, int start, float tolerance) {
	Spectrum orbcal = orbcals[index];
	vector<float> perfect_orbcal;
	getMassesCummulative(orbcal_pep.getline(index), perfect_orbcal, 0.0);
	
	proteinToOrbcal[start].insert(index);
	for (int i = 1; i < perfect_orbcal.size(); i++) {
		if (! hasPeak(&orbcal, perfect_orbcal[i], perfect_orbcal[perfect_orbcal.size()-1] - perfect_orbcal[i], tolerance)) continue;
		proteinToOrbcal[start+i].insert(index);
	}
}

void OverlapOut::rememberContig(int index, float range, Spectrum* PepSpec) {
	Spectrum contig = contigs[index];
	Spectrum overlap = overlaps[index];
	float pep_mass;
	float contig_mass;
	int pepI;
	int contI;
	
	proteinToContig[(int)(overlap[0][1]+0.01)].insert(index);
	
	for (int i = 1; i < overlap.size(); i++) {
		pepI = (int)(overlap[i][1]+0.01);
		contI = (int)(overlap[i][0]+0.01);
		pep_mass = (*PepSpec)[pepI][0] - (*PepSpec)[(int)(overlap[i-1][1]+0.01)][0];
		contig_mass = contig[contI][0] - contig[contI-1][0];
		if ((pep_mass >= contig_mass-range) && (pep_mass <= contig_mass+range)) {
			proteinToContig[pepI].insert(index);
		}
	}
}

bool OverlapOut::isFilteredOrbcal(int index, int start, float tolerance, int filterOrbcal, float peakTol, bool countConsec) {
	Spectrum orbcal = orbcals[index];
	vector<float> perfect_orbcal;
	set<int> contigOverlap;
	set<int>::iterator contigIt;
	list<int> filtered;
	
	getMassesCummulative(orbcal_pep.getline(index), perfect_orbcal, 0.0);
	
	for (int i = 0; i < perfect_orbcal.size(); i++) {
		if (proteinToContig[start + i].size() == 0) continue;
		mySet = proteinToContig[start + i];
		for (contigIt = mySet.begin(); contigIt != mySet.end(); contigIt ++) {
			if (contigOverlap.count(*contigIt) == 0) contigOverlap.insert(*contigIt);
		}
	}
	for (contigIt = contigOverlap.begin(); contigIt != contigOverlap.end(); contigIt ++) {
		if (numMatchingPeaks(contigs[*contigIt], orbcals[index], orbcalContShifts[index][*contigIt], peakTol, countConsec) < filterOrbcal) continue;
		filtered.push_back(*contigIt);
	}
	if (filtered.size() >= 2) {
		orbcalsFiltered[index] = true;
		return true;
	}
	return false;
}

bool OverlapOut::isFilteredContig(int index, int filterOrbcal, float peakTol, bool countConsec) {
	Spectrum overlap = overlaps[index];
	set<int> orbcalOverlap;
	set<int>::iterator orbcalIt;
	int protI;
	
	for (int i = 0; i < overlap.size(); i++) {
		protI = (int)(overlap[i][1]+0.01);
		if (proteinToOrbcal[protI].size() == 0) continue;
		mySet = proteinToOrbcal[protI];
		for (orbcalIt = mySet.begin(); orbcalIt != mySet.end(); orbcalIt ++) {
			if (! orbcalsFiltered[*orbcalIt]) continue;
			if (orbcalOverlap.count(*orbcalIt) == 0) orbcalOverlap.insert(*orbcalIt);
		}
	}
	for (orbcalIt = orbcalOverlap.begin(); orbcalIt != orbcalOverlap.end(); orbcalIt ++) {
		if (numMatchingPeaks(contigs[index], orbcals[*orbcalIt], orbcalContShifts[*orbcalIt][index], peakTol, countConsec) < filterOrbcal) continue;
		contigsFiltered[index] = true;
		return true;
	}
	return false;
}

unsigned int OverlapOut::numMatchingPeaks(Spectrum& contig, Spectrum& orbcal, float contigOffset, float peakTol, bool countConsec) {
	vector<bool> peaksMatchContig(contig.size());
	vector<bool> peaksMatchOrbcal(orbcal.size());
	int i = 0;
	for (int j = 0; j < contig.size(); j++) {
		float peak2 = contig[j][0] + contigOffset;
		while (i < orbcal.size() && orbcal[i][0] < peak2 - peakTol) i ++;
		for (int k = i; k < orbcal.size() && orbcal[k][0] >= peak2 - peakTol && orbcal[k][0] <= peak2 + peakTol; k ++) {
			peaksMatchContig[j] = true;
			peaksMatchOrbcal[k] = true;
		}
	}
	if (countConsec) {
		unsigned int max1 = 0;
		unsigned int num1 = 0;
		for (int x = 0; x < peaksMatchContig.size(); x++) {
			if (peaksMatchContig[x]) {
				num1 ++;
				if (num1 > max1) max1 = num1;
			} else num1 = 0;
		}
		
		return max1;
	} else {
		unsigned int num1 = 0;
		unsigned int num2 = 0;
		for (int x = 0; x < peaksMatchContig.size(); x++) {
			if (peaksMatchContig[x]) num1 ++;
		}
		for (int x = 0; x < peaksMatchOrbcal.size(); x++) {
			if (peaksMatchOrbcal[x]) num2 ++;
		}
		return min(num1, num2);
	}
}

bool OverlapOut::hasPeak(Spectrum* orbcal, float bpeak_, float ypeak_, float tolerance) {
	float bpeak = bpeak_;
	float ypeak = ypeak_ + AAJumps::massH2O;
	int bloc = (*orbcal).findClosest(bpeak);
	int yloc = (*orbcal).findClosest(ypeak);
	bool isy = false;
	bool isb = false;
	isb = isEqual((*orbcal)[bloc][0], bpeak, tolerance);
	isy = isEqual((*orbcal)[yloc][0], ypeak, tolerance);
	return isy || isb;
}