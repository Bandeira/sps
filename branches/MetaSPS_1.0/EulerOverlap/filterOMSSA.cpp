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
#include "../batch.h"

using namespace std;

const char* usage = "\n\nTo filter out all omssa search results by a given p-value and output each matched peptide to text file indexed by line number, run:\n\n./filterOMSSA [omssa results] [p-value] [output .txt file]\n\noptimum p-value is 0.001\n";

int main(int argc, char *argv[], char **envp ) {
	if (argc != 4) {cout << usage; return 0;}
	loadCsvToTxt(argv[1], argv[3], getFloat(argv[2]));
	return 0;
}
