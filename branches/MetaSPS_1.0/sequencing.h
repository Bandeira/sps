#ifndef SEQUENCING_H
#define SEQUENCING_H

#include "alignment_scoring.h"
#include "spectral_pairs.h"
#include "batch.h"
#include "filters.h"
#include "abruijn.h"
#include "graph.h"
#include "setmerger.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>


struct masab_parameters {
	SpecSet INPUT_SPECS;
	vector<Results_PA> INPUT_ALIGNS;
	
	float PENALTY_PTM;
	float PENALTY_SAME_VERTEX;
	int GRAPH_TYPE;

	vector<SpectrumPeakLabels>* LABELS;
	int MAX_AA_JUMP;
	float MAX_MOD_MASS;
	float TOLERANCE_PEAK;
	float TOLERANCE_PM;
	short EDGE_SCORE_TYPE;
	unsigned int MIN_MATCHED_PEAKS;
	unsigned int MIN_EDGES_TO_COMPONENT;
	unsigned int PATH_MIN_SPECS;
	short PATH_MIN_PEAKS;
	int SPEC_TYPE_MSMS;

	bool NO_SEQUENCING;
	bool ADD_ENDPOINTS;
	bool OUTPUT_COMPLETE_ABRUIJN;

	SetMerger COMPONENTS;
	vector<vector<float> > COMPONENT_STATS;
	vector<list<int> > COMPONENT_SPECTRA;
	abinfo_t COMPONENT_INFO;
	abinfo_t COMPLETE_ABRUIJN;
	Clusters PATH_SPECTRA_AS_CLUSTER;
	vector<vector<short> > ABCOUNTS;
	SpecSet OUTPUT_SPECS;
};

void SplitASPPA(vector<Results_PA> &pairs, float pmTol, vector<Results_ASP> &pairsASP, vector<Results_PA> &pairsPA);

void Masab(masab_parameters& m_parameters, ostream& STDOUT = cout);

#endif
