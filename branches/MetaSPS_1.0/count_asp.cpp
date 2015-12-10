#include "aminoacid.h"
#include "spectrum.h"
//#include "batch.h"

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv){
    vector<PairAlign> aligns(0);
    float score;

    if (argc<5) { 
        cerr << "Usage: count_asp <startBaseIdx> <endBaseIdx> <prmSpecs.pkl> <output.txt>\n";
        cerr << "       <startBaseIdx>, <endBaseIdx> are numbers from 1 to number of spectra in <prmSpecs.pkl>\n";
        return(-1);
    }
    
    int startBaseIdx=(int)atof(argv[1]), endBaseIdx=(int)atof(argv[2]);
    vector<int> baseSpecIdx(endBaseIdx-startBaseIdx+1);  for(int i=0; i<endBaseIdx-startBaseIdx+1; i++) baseSpecIdx[i]=startBaseIdx+i-1;
    
    SpecSet specSet;
    specSet.LoadSpecSet_pkl(argv[3]);
    cout << "Loading specs complete. Num specs: " << specSet.size() << "\n";

//    countPairAlignsASP(specSet, baseSpecIdx, 3, 1, 0.5, results);  // Original version
    
    int i, j, massIdx;
    double mass;
    vector<list<int> > hist;   specSet.getMassesHistogram(hist,resolution);
    
    // Determine parent mass differences for valid matches
    AAJumps jumps(aaDiff);    jumps.forceJump(0);   jumps.forceTolerance(pmTol);    jumps.forceDoubleSided();
    jumps.forceUnique();
    
    // Count number of pairs
    int numPairs = 0;
    for(i=0; i<baseSpectraIdx.size(); i++) {
        mass = specSet[baseSpectraIdx[i]].parentMass;
        for(j=0, massIdx=(int)round((mass+jumps[j])/resolution); j<jumps.size() && massIdx<0 ; j++, massIdx=(int)round((mass+jumps[j])/resolution));  // Skips masses<0
        for(; j<jumps.size() && massIdx<hist.size(); j++, massIdx=(int)round((mass+jumps[j])/resolution)) 
            for(curPair = hist[massIdx].begin(); curPair!=hist[massIdx].end(); curPair++) if (*curPair>baseSpectraIdx[i]) numPairs++;
    }
cout << "Total number of pairs to match: " << numPairs << endl;
    
    return 1;
}