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

const char* usage = "To filter matchma results, run:\n./filterMatchma [contig indicies] [in prot match] [out prot match]\n";

int main(int argc, char *argv[], char **envp ) {
	if (argc != 4) {cout << usage; return 0;}
	
	vector<vector<int> > prot_match;
	vector<vector<int> > contig_idx;
	
	if (! Load_binArray(argv[1], contig_idx)) cerr << "ERROR: Failed to load " << argv[1] << "\n";
	if (! Load_binArray(argv[2], prot_match)) cerr << "ERROR: Failed to load " << argv[2] << "\n";
	
	vector<vector<int> > prot_matchOut(prot_match);
	set<int> filtered;
	for (int i = 0; i < contig_idx.size(); i++) filtered.insert(contig_idx[i][1] - 1);
	for (int i = 0; i < prot_matchOut.size(); i++) {
		if (filtered.count(i) == 0) prot_matchOut[i][0] = -1;
	}
	if (! Save_binArray(argv[3], prot_matchOut)) cerr << "ERROR: Failed to save " << argv[3] << "\n";
	return 0;
}
