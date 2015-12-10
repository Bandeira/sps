/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 10/26/09
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

const char* usage = "To merge contig components and contigs, run:\n./mergeContigComponentSets [orig contigs] [orig overlaps] [orig prot_match] [component .bla] [component contigs] [component overlaps] [component prot_match] [out contigs] [out overlaps] [out prot_match]\n";

int main(int argc, char *argv[], char **envp ) {
	if (argc != 11) {cout << usage; return 0;}
	
	SpecSet contigs;
	SpecSet contigsComp;
	SpecSet contigsOut;
	SpecSet overlaps;
	SpecSet overlapsComp;
	SpecSet overlapsOut;
	vector<vector<int> > prot_match;
	vector<vector<int> > prot_matchComp;
	
	vector<list<int> > comps;
	
	
	if (! contigs.LoadSpecSet_pklbin(argv[1])) {cerr << "ERROR: Failed to load " << argv[1] << "\n"; return 1;}
	if (! contigsComp.LoadSpecSet_pklbin(argv[5])) {cerr << "ERROR: Failed to load " << argv[5] << "\n"; return 1;}
	if (! overlaps.LoadSpecSet_pklbin(argv[2])) {cerr << "ERROR: Failed to load " << argv[2] << "\n"; return 1;}
	if (! overlapsComp.LoadSpecSet_pklbin(argv[6])) {cerr << "ERROR: Failed to load " << argv[6] << "\n"; return 1;}
	if (! Load_binArray(argv[3], prot_match)) {cerr << "ERROR: Failed to load " << argv[3] << "\n"; return 1;}
	if (! Load_binArray(argv[7], prot_matchComp)) {cerr << "ERROR: Failed to load " << argv[7] << "\n"; return 1;}
	if (! Load_binListArray<int,list<int>,list<int>::iterator>(argv[4],comps)) {cerr << "ERROR: Failed to load " << argv[4] << "\n"; return 1;}
	
	vector<vector<int> > prot_match_out(prot_matchComp);
	contigsOut = contigsComp;
	overlapsOut = overlapsComp;
	
	int size = contigs.specs.size();
	//cout << size << endl;
	prot_match_out.resize(size);
	contigsOut.specs.resize(size);
	overlapsOut.specs.resize(size);
	
	set<int> used;
	for (int i = 0; i < comps.size(); i++) {
		//cout << i << ": ";
		for (list<int>::iterator node = comps[i].begin(); node != comps[i].end(); node ++) {
			used.insert(*node);
			//cout << *node << " ";
		}
		//cout << endl;
	}
	
	
	int idxUse = prot_matchComp.size();
	//cout << idxUse << endl;
	for (int i = 0; i < contigs.size(); i++) {
		if (used.count(i) > 0) continue;
		//cout << idxUse << " ";
		//cout.flush();
		contigsOut[idxUse] = contigs[i];
		overlapsOut[idxUse] = overlaps[i];
		prot_match_out[idxUse] = prot_match[i];
		idxUse ++;
	}
	
	contigsOut.specs.resize(idxUse);
	overlapsOut.specs.resize(idxUse);
	prot_match_out.resize(idxUse);
	contigsOut.SaveSpecSet_pklbin(argv[8]);
	overlapsOut.SaveSpecSet_pklbin(argv[9]);
	Save_binArray(argv[10], prot_match_out);
	
	return 0;
}
