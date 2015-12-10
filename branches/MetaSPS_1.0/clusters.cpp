#include "clusters.h"
#include "aminoacid.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

void Clusters::changeResolution(float newResolution, bool enforceRounding) {
	unsigned int s,i;
	if (enforceRounding)
		for(s=0; s<consensus.size(); s++) {
			consensus[s].changeResolution(newResolution,enforceRounding);
			for(i=0; i<shifts[s].size(); i++) shifts[s][i][0]=round(shifts[s][i][0]/newResolution);
		}
	else
		for(s=0; s<consensus.size(); s++) {
			consensus[s].changeResolution(newResolution,enforceRounding);
			for(i=0; i<shifts[s].size(); i++) shifts[s][i][0]=shifts[s][i][0]/newResolution;
		}
}

void Clusters::set_endpoints(int idx, SpecSet &specs, vector<list<TwoValues<int> > > &abVertices,
		float peakTol, float pmTol) {
	float peakMass, scoreBep, scoreYep;
	list<TwoValues<int> >::iterator vertexIter;
	unsigned int setIdx, specIdx, peakIdx, endpointIdx;

	if(consensus[idx].size()!=abVertices.size()) {
		cerr<<"[Warning in Clusters::set_endpoints()] Cannot set endpoint masses because non-matching number of consensus peaks and ABruijn vertices.";
		endpoints[idx].resize(0);
		return;
	}

	endpoints[idx].resize(consensus[idx].size());   endpointIdx = 0;
	for(setIdx=0; setIdx<abVertices.size(); setIdx++) {
		scoreBep = 0;   scoreYep = 0;
		for(vertexIter=abVertices[setIdx].begin(); vertexIter!=abVertices[setIdx].end(); vertexIter++) {
			specIdx = (*vertexIter)[0];   peakIdx = (*vertexIter)[1];   peakMass = specs[specIdx][peakIdx][0];
			if(abs(peakMass)<=peakTol or abs(peakMass+AAJumps::massMH-specs[specIdx].parentMass)<2*pmTol) scoreBep+=specs[specIdx][peakIdx][1];
			if(abs(peakMass-AAJumps::massH2O)<=peakTol or abs(peakMass+AAJumps::massHion-specs[specIdx].parentMass)<2*pmTol) scoreYep+=specs[specIdx][peakIdx][1];
		}
		if(scoreBep+scoreYep<0.00001) endpoints[idx][setIdx]=0;
		else { if(scoreBep>=scoreYep) endpoints[idx][setIdx]=1; else endpoints[idx][setIdx]=2; }
	}
}

void Clusters::getSpecIfB(int idx, Spectrum &putHere, float peakTol, float ionOffset) {
	float deltaY0=0;  // Mass added to all peaks when the first peak is y0 (shift right by AAJumps::massH2O)

	// Adjust masses for peaks that are endpoint masses in the assembled spectra
	putHere = consensus[idx];
	if(not endpoints[idx].empty())
		for(unsigned int peakIdx=0; peakIdx<putHere.size(); peakIdx++) {
			putHere[peakIdx][0]+=deltaY0;
			if(endpoints[idx][peakIdx]==2) {
				if(peakIdx==0) {
					deltaY0 = AAJumps::massH2O;
					putHere.parentMass+=AAJumps::massH2O;  // First mass was y0, need to adjust parent mass (for addZPMpeaks)
				} else putHere[peakIdx][0]-=AAJumps::massH2O;
				if(peakIdx==putHere.size()-1) putHere.parentMass-=AAJumps::massH2O;  // Last mass was yN, need to adjust parent mass (for addZPMpeaks)
			}
		}

	putHere.addZPMpeaks(peakTol,ionOffset,false); // Guarantee peaks at b0/bN with default intensity
}

void Clusters::getSpecIfY(int idx, Spectrum &putHere, float peakTol, float ionOffset) {
	Spectrum tmp = consensus[idx];
	float aaMass=tmp.parentMass-AAJumps::massMH;

	// Adjust masses for peaks that are endpoint masses in the assembled spectra
	if(not endpoints[idx].empty())
		for(unsigned int peakIdx=0; peakIdx<tmp.size(); peakIdx++) {
			if(endpoints[idx][peakIdx]==1) {
				tmp[peakIdx][0]+=AAJumps::massH2O;
				if(peakIdx==0) tmp.parentMass-=AAJumps::massH2O;  // First mass was b0, need to adjust parent mass (for addZPMpeaks)
				if(peakIdx==putHere.size()-1) {
					tmp.parentMass+=AAJumps::massH2O;  // Last mass was bN, need to adjust parent mass (for addZPMpeaks)
					aaMass+=AAJumps::massH2O;  // Last mass was bN, need to adjust symmetry
				}
			}
		}

	// Reverse merged spectrum
	int topIdx=tmp.size()-1;  while(topIdx>=0 && tmp[topIdx][0]>aaMass) topIdx--;
	putHere.peakList.resize(topIdx+1);    putHere.parentMass = tmp.parentMass;
	for(int i=0; i<=topIdx; i++) putHere[topIdx-i].set(aaMass-tmp[i][0],tmp[i][1]);

	putHere.addZPMpeaks(peakTol,ionOffset,false); // Guarantee peaks at b0/bN
}


//
//  File format is:
//  <numClusters>
//  Repeat per cluster:
//    <numEntries>
//    Repeat per spectrum
//        <specIdx> <shift1> <shift2>
//    <num consensus PRMs> <mhMass>
//    Repeat per consensus PRM
//        <mass> <score>
//    <numEndpoints>
//    Repeat per endpoint
//        <endpoint types> (0=internal, 1=b0/bN, 2=y0/yN)
//        <mass> <score>   <<--- discontinued old format

int Clusters::Load(const char *contigFilename,
        const char *spectrumPklbinFilename /* = 0x0 */,
        const char *psmFilename /* = 0x0 */){
	ifstream input(contigFilename);
			if (!input) {
	       cerr << "Can not open  [" << contigFilename << "]\n";
			  specIdx.resize(0);
			  return -1;
			}

			unsigned int count, cluster, i;

			input >> count;  specIdx.resize(count);  shifts.resize(count);  consensus.resize(count);  endpoints.resize(count);
			for(cluster=0; cluster<specIdx.size(); cluster++) {
					// Read spectrum indices and shifts
					input >> count;   specIdx[cluster].resize(count);   shifts[cluster].resize(count);
					for(i=0; i<count; i++) input >> specIdx[cluster][i] >> shifts[cluster][i][0] >> shifts[cluster][i][1];

					// Read consensus spectrum PRMs
					input >> count >> consensus[cluster].parentMass;   consensus[cluster].resize(count);
					for(i=0; i<count; i++) input >> consensus[cluster][i][0] >> consensus[cluster][i][1];

					// Read endpoints
					input >> count;   endpoints[cluster].resize(count);
	//        for(i=0; i<count; i++) input >> endpoints[cluster][i][0] >> endpoints[cluster][i][1];
					for(i=0; i<count; i++) input >> endpoints[cluster][i];
			}
			input.close();

	    if (spectrumPklbinFilename) {
	      if (consensus.LoadSpecSet_pklbin("assembly/sps_seqs.pklbin") <= 0) {
	        return -2;
	      }
	    }
	    if (psmFilename) {/*
	      PeptideSpectrumMatchSet psmSetTemp;
	      if (psmSetTemp.loadFromFile("assembly/tagsearchpsm.txt")) {
		      psmSetTemp.addSpectra(&consensus, true);
		    } else {
		    	cerr << "Problem loading PSM file [" << psmFilename << "]\n";
	        return -3;
	      }*/
	    }

			return(specIdx.size());
}

int Clusters::Load_pklbin(const char *filename){
    unsigned int count, cluster, i;

    if(not consensus.LoadSpecSet_pklbin(filename)) { specIdx.resize(0); shifts.resize(0); endpoints.resize(0); return -1; }

    specIdx.resize(consensus.size());  shifts.resize(consensus.size());  endpoints.resize(consensus.size());
    for(cluster=0; cluster<consensus.size(); cluster++) {
    	specIdx[cluster].resize(0);
    	shifts[cluster].resize(0);
    	endpoints[cluster].resize(0);
    }

	return(specIdx.size());
}

int Clusters::Save(const char *filename){
    ofstream output(filename);
    if (!output) { cerr << "ERROR: cannot open " << filename << "\n";   return -1; }

    unsigned int count, cluster, i;

    output << specIdx.size() << endl;
    for(cluster=0; cluster<specIdx.size(); cluster++) {
        output << specIdx[cluster].size() << endl;
        for(i=0; i<specIdx[cluster].size(); i++) output << specIdx[cluster][i] << " " << shifts[cluster][i][0] << " " << shifts[cluster][i][1] << endl;
        output << consensus[cluster].size() << " " << consensus[cluster].parentMass << endl;
        for(i=0; i<consensus[cluster].size(); i++) output << consensus[cluster][i][0] << " " << consensus[cluster][i][1] << endl;
        output << endpoints[cluster].size() << endl;
//        for(i=0; i<endpoints[cluster].size(); i++) output << endpoints[cluster][i][0] << " " << endpoints[cluster][i][1] << endl;
        for(i=0; i<endpoints[cluster].size(); i++) output << endpoints[cluster][i] << endl;
    }

    output.close();  return(0);
}
