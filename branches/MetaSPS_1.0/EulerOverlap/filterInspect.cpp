/**
 * Author: Adrian Guthals
 * Email: adrian.guthals@gmail.com
 * Date: 11/5/09
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
#include "../inspect_parse.h"

using namespace std;

const char* usage = "\n\nTo annotate spectra and output matched peptide to text file indexed by line number, run:\n\n./filterInspect [inspect results] [searched spectra] [annotation model] [output .txt file]\n\n";

int main(int argc, char *argv[], char **envp ) {
	if (argc != 5) {cout << usage; return 0;}
    SpecSet in_specs;
    if (! in_specs.LoadSpecSet_pklbin(argv[2]) ) {
        cerr << "ERROR: Failed to load " << argv[2] << "\n";
        return -1;
    }
    
    MS2ScoringModel model;
    if (! model.LoadModel(argv[3]) ) {
		cerr << "ERROR: Failed to load " << argv[3] << "\n";
        return -1;
    }
    InspectResultsSet in_results;
	if (! in_results.loadInspectResultsFile(argv[1]) ) {
		cerr << "ERROR: Failed to load " << argv[1] << "\n";
        return -1;
	}
    string ionNamesInclude("all");
	int num_annotated = 0;
	for (int i = 0; i < in_specs.size(); i++) {
		if (in_results.results.count(i) == 0) continue;
		string specnet_annot;
		in_results.inspectToSpecNets(in_results.results[i].Annotation, specnet_annot);
		in_specs[i].annotate(specnet_annot, ionNamesInclude, model, 0, 0, 0.5);
		num_annotated ++;
	}
	cout << "annotated " << num_annotated << " of " << in_specs.size() << " spectra\n";
    if (! in_specs.SaveAnnotations(argv[4]) ) return -1;
	return 0;
}
