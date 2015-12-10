////////////////////////////////////////////////////////////////////////////////
#include "spectrum.h"
#include "PWizInterface.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstring>

////////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;
////////////////////////////////////////////////////////////////////////////////
// Dump specset data (Debug)
////////////////////////////////////////////////////////////////////////////////
void dumpSpecset(SpecSet &specSet)
{
  if(!specSet.size()) {
    cout << "Specset size is 0" << endl;
    return;
  }

  specnets::Spectrum &spectrum = specSet[0];

  if(!spectrum.size()) {
    cout << "spectrum size is 0" << endl;
    return;
  }


  cout << "parentMassTol :       "  << spectrum.parentMassTol       << endl;
  cout << "parentMZ :            "  << spectrum.parentMZ            << endl;
  cout << "parentMass :          "  << spectrum.parentMass          << endl;
  cout << "parentCharge :        "  << spectrum.parentCharge        << endl;
  cout << "scan :                "  << spectrum.scan                << endl;
  cout << "msLevel :             "  << spectrum.msLevel             << endl;
  cout << "msFragType :          "  << spectrum.msFragType          << endl;
  cout << "fileName :            "  << spectrum.fileName            << endl;
  cout << "resolution :          "  << spectrum.resolution          << endl;
  cout << "instrument_name :     "  << spectrum.instrument_name     << endl;
  cout << "ITOL :                "  << spectrum.ITOL                << endl;
  cout << "ITOLU :               "  << spectrum.ITOLU               << endl;
  cout << "TOL :                 "  << spectrum.TOL                 << endl;
  cout << "TOLU :                "  << spectrum.TOLU                << endl;
  cout << "spectrum_quality :    "  << spectrum.spectrum_quality    << endl;
  cout << "idDist :              "  << spectrum.idDist              << endl;
  cout << "retention_time :      "  << spectrum.retention_time      << endl;
  cout << "precursor_intensity : "  << spectrum.precursor_intensity << endl;
}
////////////////////////////////////////////////////////////////////////////////
// Entry point
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  SpecSet specs;
  string  outName;
  FilenameManager fm;
  vector<TwoValues<unsigned int> > specsInfo; // Scan num (pos.0), msLevel (pos.1)


  if (argc < 3) {
    // if less than 2 args, that means args are missing
    if (argc < 2) {

      cerr << "Syntax: " << argv[0] << " <file_type> <input_file_name> <output_file_name_prefix>" << endl;
      cerr
          << "        file_type: mgf, ms2, pkl (multiple spectra per file), prms (pepnovo_prm output format)\n";
      return -1;

      // if 2 args, let's see if only the filename was used. In that case, let's get the extension from the filename
    } else {
      fm.filenameFull   = argv[1];
      fm.splitFilename();
    }

  // all
  } else {
    fm.filenameFull  = argv[2];
    fm.splitFilename();
    fm.extension = argv[1];
  }

  // Generate output filename
  if (argc <= 3) {
    size_t found   = fm.filenameFull.find_last_of(".");
    outName        = fm.filenameFull.substr(0, found);
    outName       += ".pklbin";
  } else {
    outName = argv[3];
  }

  // load the specset
  bool loadOk = specs.Load(fm.filenameFull.c_str());


  if (loadOk) {
    specs.savePklBin(outName.c_str());
  } else {
    cerr << "ERROR parsing file_type; valid options are ms2, mgf, mzml, mzxml, pkl or prms.\n";
    return -1;
  }

  //dumpSpecset(specs);

  return (0);
}
////////////////////////////////////////////////////////////////////////////////
