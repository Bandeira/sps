/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 6/29/09
 */

#include "PruneEnvelopes.h"

using namespace std;

PruneEnvelopes::PruneEnvelopes() {}
PruneEnvelopes::~PruneEnvelopes() {}

PruneEnvelopes::PruneEnvelopes(const char* filename) {
	if (! env.LoadModel(filename) ) cerr << "ERROR: Failed to load " << filename << "\n";
}

PruneEnvelopes::PruneEnvelopes(IsoEnvelope& env2) {
	env = env2;
}

void PruneEnvelopes::getPrunedSpectrum(Spectrum& spectrum, Spectrum* putHere, vector<TwoValues<float> >* putChargeScore, float peakTol, float threshold, bool strictMode, bool usePPM) {
	vector<TwoValues<float> > chargeScore(spectrum.size());
	Spectrum spec = spectrum;
	Spectrum temp = spec;
	vector<list<unsigned int> > holdIndicies(spec.size());
	vector<vector<float> > holdEnvelopes(spec.size());
	vector<list<unsigned int> > indicesUsed;
	vector<list<unsigned int> > indicesUsed2;
	list<unsigned int> indicesUsedList;
	list<unsigned int>::iterator indiciesIt;
	vector<float> massEnvelope;
	vector<float> massEnvelope2;
	float mass, score, monoisotopicmass, intensity;
	set< int > dontuse;
	//FILE* output = fopen("debug_envelope.txt", "w");
	(*putHere) = spec;
	int putIdxUsed = 0;
	(*putHere).peakList.resize(spec.size());
	(*putChargeScore).resize(spec.size());
	int index = 0;
	while (true) {
		if (index >= spec.size()) break;
		mass = spec[index][0];
		intensity = spec[index][1];
		chargeScore[index][0] = -1.0;
		chargeScore[index][1] = -1.0;
		for (unsigned short charge = spectrum.parentCharge; charge > 1; charge --) {
			monoisotopicmass = getMonoisotopicmass(mass, charge);
			if (monoisotopicmass >= spec.parentMass) continue;
			massEnvelope.clear();
			float pkTolUse = (usePPM) ? ((mass * InputParams::PPM)/1000000.0) : peakTol;
			env.ExtractEnvelope(mass, charge, spec, pkTolUse, massEnvelope, indicesUsed, strictMode);
			if (massEnvelope[0] == 0 || massEnvelope[1] == 0) continue;
			for (int i = massEnvelope.size(); i > 1; i--) {
				if (i < massEnvelope.size() && indicesUsed[i].size() == 0) continue;
				if (i < massEnvelope.size()) {
					massEnvelope[i] = 0.0;
					indicesUsed[i].clear();
				}
				
				bool badcharge = false;
				for (int j = 0; j < indicesUsed.size(); j++) {
					for (indiciesIt = indicesUsed[j].begin(); indiciesIt != indicesUsed[j].end(); indiciesIt ++) {
						if (spec[*indiciesIt][1] > 12.0*intensity) {
								badcharge = true;
								break;
						}
					}
					if (badcharge) break;
				}
				if (badcharge) continue;
				vector<float> normalizedEnvelope(massEnvelope);
				env.normalize(normalizedEnvelope, true);
				score = env.ScoreEnvelope(monoisotopicmass, normalizedEnvelope, strictMode);
				/*
				if (isEqual(mass, 733.039, 0.001)) {
					ostringstream out;
					float inten = intensity;
					for (int j = 0; j < indicesUsed.size(); j++) {
						for (indiciesIt = indicesUsed[j].begin(); indiciesIt != indicesUsed[j].end(); indiciesIt ++) {
							inten += spec[*indiciesIt][1];
						}
					}
					out << "mass: " << parseFloat(mass, 3) << "\tmonoisotopicmass: " << parseFloat(monoisotopicmass, 3) << "\tintensity: " << parseFloat(inten, 3) << "\tcharge: " << charge << "\tscore: " << score << "\tenvelope idx:";
					for (int j = 0; j < indicesUsed.size(); j++) {
						for (indiciesIt = indicesUsed[j].begin(); indiciesIt != indicesUsed[j].end(); indiciesIt ++) {
							out << " " << *indiciesIt << "(" << spec[*indiciesIt][0] << ")";
						}
					}
					out << endl;
					fprintf(output, out.str().c_str());
				}
				*/
				if (1.0/score > chargeScore[index][1] && 1.0/score >= (1.0/threshold) - 0.00001) {
					chargeScore[index][0] = (float)charge;
					chargeScore[index][1] = 1.0/score;
					indicesUsedList.clear();
					for (int j = 0; j < indicesUsed.size(); j++) {
						for (indiciesIt = indicesUsed[j].begin(); indiciesIt != indicesUsed[j].end(); indiciesIt ++) {
							indicesUsedList.push_back(*indiciesIt);
						}
					}
					holdIndicies[index] = indicesUsedList;
					holdEnvelopes[index] = massEnvelope;
				}
			}
		}
		if (chargeScore[index][1] < -0.0001) {
			(*putHere)[putIdxUsed] = spec[index]; putIdxUsed ++; index ++;
			continue;
		}

		unsigned short charge = (unsigned short)round(chargeScore[index][0]);
		cout << "Charge at " << mass << " is " << charge << " with score " << 1.0/chargeScore[index][1] << " : ";

		twovalues[0] = chargeScore[index][0];
		twovalues[1] = 1.0/chargeScore[index][1];
		(*putChargeScore)[putIdxUsed] = twovalues;
		twovalues[0] = getMonoisotopicmass(mass, charge);
		intensity = 0.0;
		for (int j = 0; j < holdEnvelopes[index].size(); j++) intensity += holdEnvelopes[index][j];
		twovalues[1] = intensity;
		dontuse.clear();
		for (indiciesIt = holdIndicies[index].begin(); indiciesIt != holdIndicies[index].end(); indiciesIt ++) {
		  cout << spec[(int)(*indiciesIt)][0] << ", ";
		  dontuse.insert((int)(*indiciesIt));
		}
		cout << "\n";
		/*
		ostringstream out;
		out << "mass: " << parseFloat(mass, 3) << "\tmonoisotopicmass: " << parseFloat(twovalues[0], 3) << "\tintensity: " << parseFloat(intensity, 3) << "\tcharge: " << chargeScore[index][0] << "\tscore: " << 1.0/chargeScore[index][1] << "\tenvelope idx:";
		for (indiciesIt = holdIndicies[index].begin(); indiciesIt != holdIndicies[index].end(); indiciesIt ++) {
			out << " " << *indiciesIt << "(" << spec[*indiciesIt][0] << ")";
		}
		out << endl;
		fprintf(output, out.str().c_str());
		*/
		(*putHere)[putIdxUsed] = twovalues; putIdxUsed++;
		temp.peakList.clear();
		temp.peakList.resize(spec.size());
		int idxUsed = 0;
		int allIdx = 0;
		for (; allIdx < index; allIdx++) {temp[idxUsed] = spec[allIdx]; idxUsed++;}
		allIdx++;
		for (; allIdx < spec.size(); allIdx++) {
			if (dontuse.count(allIdx) == 0) {temp[idxUsed] = spec[allIdx]; idxUsed++;}
		}
		temp.peakList.resize(idxUsed);
		spec = temp;
	}
	(*putHere).peakList.resize(putIdxUsed);
	(*putChargeScore).resize(putIdxUsed);
	sort((*putHere).peakList.begin(), (*putHere).peakList.end());
	//fclose(output);
}

bool PruneEnvelopes::uploadPklBin(const char* filename, float peakTol, float threshold, bool strictMode, bool usePPM) {
	if (! pklSpecs.LoadSpecSet_pklbin(filename)) {
		cerr << "ERROR: Cannot open file " << filename << endl;
		return false;
	}
	transferAllSpecs(peakTol, threshold, strictMode, usePPM);
	return true;
}

bool PruneEnvelopes::uploadMGF(const char* filename, float peakTol, float threshold, bool strictMode, bool usePPM) {
	if (! pklSpecs.LoadSpecSet_mgf(filename)) {
		cerr << "ERROR: Cannot open file " << filename << endl;
		return false;
	}
	transferAllSpecs(peakTol, threshold, strictMode, usePPM);
	return true;
}

void PruneEnvelopes::transferAllSpecs(float peakTol, float threshold, bool strictMode, bool usePPM) {
	//cout << peakTol << ", " << threshold << ", " << usePPM << "\n";
	vector<TwoValues<float> > putChargeScore;
	Spectrum prunedSpec;
	prunedSpecs.specs.clear();
	specChargeScore.clear();
	prunedSpecs.specs.resize(pklSpecs.size());
	specChargeScore.resize(pklSpecs.size());
	ProgressDisplay prog(cout, pklSpecs.size()-1, "Pruning Envelopes: % ");
	for (int idx = 0; idx < pklSpecs.size(); idx ++) {
		prog.showProgress(cout, idx);
		getPrunedSpectrum(pklSpecs[idx], &prunedSpec, &putChargeScore, peakTol, threshold, strictMode, usePPM);
		prunedSpecs.specs[idx] = prunedSpec;
		specChargeScore[idx] = putChargeScore;
		prunedSpec.output(cout);
		break;
	}
	prog.clear(cout);
}

float PruneEnvelopes::getMonoisotopicmass(float mass, unsigned short charge) {
	return (mass*((float)((short)charge))) - ((float)((short)charge))+AAJumps::massHion;
}

