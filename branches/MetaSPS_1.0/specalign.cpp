#include "inputParams.h"
#include "alignment_scoring.h"
#include "batch.h"
#include "spectral_alignment.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;

int main(int argc, char **argv){
    
    InputParams params; bool paramsOk;
	if(argc<=1) paramsOk=params.readParams("specalign.params");	else paramsOk=params.readParams(argv[1]);
	if(!paramsOk) {
		cerr << "Error opening parameters file ";
		if(argc<=1) cerr << "specalign.params\n"; else cerr << argv[1] << endl;
		return -1;
	}
	vector<char *> paramStrs;   paramStrs.resize(5);
	paramStrs[0] = "INPUT_ALIGNS";
	paramStrs[1] = "OUTPUT_SPECS";
	paramStrs[2] = "OUTPUT_MODPOS";
	paramStrs[3] = "PENALTY_PTM";
	paramStrs[4] = "PENALTY_SAME_VERTEX";
	if(!params.confirmParams(paramStrs)) {
		cerr << "ERROR: Parameters file ";
		if(argc==1) cerr<<"specalign.params"; else cerr<<argv[1];
		cerr << " is incomplete. One of the following is missing: INPUT_ALIGNS, OUTPUT_SPECS, OUTPUT_MODPOS, PENALTY_PTM, PENALTY_SAME_VERTEX\n";
		return -1;
	}
	
	char *alignsFN = params.getValue("INPUT_ALIGNS");
	char *resultsFN = params.getValue("OUTPUT_SPECS");
	char *modPosFN = params.getValue("OUTPUT_MODPOS");
	float penalty_ptm = (float) params.getValueDouble("PENALTY_PTM");
	float penalty_sameVert = (float) params.getValueDouble("PENALTY_SAME_VERTEX");

	int startBaseIdx = params.getValue("IDX_START")?params.getValueInt("IDX_START"):0;
	int endBaseIdx = params.getValue("IDX_END")?params.getValueInt("IDX_END"):0;
	int maxAAjump = params.getValue("MAX_AA_JUMP")?params.getValueInt("MAX_AA_JUMP"):0;
	float peakTol = params.getValue("TOLERANCE_PEAK")?(float) params.getValueDouble("TOLERANCE_PEAK"):0.5;
	float resolution = params.getValue("RESOLUTION")?(float) params.getValueDouble("RESOLUTION"):0.1;
    int specType = params.paramPresent("SPEC_TYPE_MSMS")?((int) params.getValueInt("SPEC_TYPE_MSMS")?1:0):0;
	float ionOffset = specType?AAJumps::massHion:0;
    
    SpecSet specSet;   short loadOk=0;   unsigned int specIdx;
    if(params.paramPresent("INPUT_SPECS")) loadOk=specSet.LoadSpecSet_pkl(params.getValue("INPUT_SPECS"));
    else if(params.paramPresent("INPUT_SPECS_PKLBIN")) loadOk=specSet.LoadSpecSet_pklbin(params.getValue("INPUT_SPECS_PKLBIN"));
    if (loadOk<=0 or specSet.size()==0) return -1;

	vector<Results_ASP> aligns;   
	if (!Load_resultsASPbin(alignsFN, aligns)) { cerr << "Error reading "<<alignsFN<<"!\n"; return -1; }
    if(endBaseIdx==0) endBaseIdx = aligns.size()-1;

	SpecSet resultSpecsSet;   resultSpecsSet.specs.resize(2*aligns.size());  // 2 results per pair of spectra
	vector<float> modPositions;   modPositions.resize(aligns.size());        // Keeps track of the mass where the mod was placed
	int idxResultSpecs=0, spec1, spec2, idxPair;
	for (idxPair=startBaseIdx; idxPair<=endBaseIdx; idxPair++) {
		spec1 = aligns[idxPair].spec1;    spec2 = aligns[idxPair].spec2;
		resultSpecsSet[idxResultSpecs].copyNP(specSet[spec1]);
		resultSpecsSet[idxResultSpecs+1].copyNP(specSet[spec2]);

		modPositions[idxPair]=spec_align(&specSet[spec1],&specSet[spec2],peakTol,&resultSpecsSet[idxResultSpecs],&resultSpecsSet[idxResultSpecs+1],maxAAjump,penalty_sameVert,penalty_ptm);
		idxResultSpecs+=2;
	}

	resultSpecsSet.SaveSpecSet_pklbin(resultsFN);
	Save_binArray(modPosFN, modPositions);
	return(0);
}
