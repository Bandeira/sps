////////////////////////////////////////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include <string>

#include "PWizInterface.h"

#if defined(__linux__)
#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"

using namespace pwiz::msdata;
using namespace pwiz::cv;

#endif

////////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace specnets;

////////////////////////////////////////////////////////////////////////////////
// Load data using proteowizard section
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
#if defined(__linux__)
////////////////////////////////////////////////////////////////////////////////
// Find the index of a given ID in pwiz data structure for controlled vocabulary
////////////////////////////////////////////////////////////////////////////////
int findIndex(SpectrumPtr &spectrum, int id)
{
  for(int i = 0 ; i < spectrum->cvParams.size() ; i++) {
    if(spectrum->cvParams[i].cvid == id)
      return i;
  }
  return -1;
}
////////////////////////////////////////////////////////////////////////////////
// Get the scan from the ID sting
////////////////////////////////////////////////////////////////////////////////
bool getScan(SpectrumPtr &spectrum, int &scan)
{
  // test if id contains zomething
  if(!spectrum->id.length())
    return -1;

  // hold strings
  vector<string> aux;
  stringSplit(spectrum->id, aux, " \t");

  for(int i = 0 ; i < aux.size() ; i++) {

    vector<string> aux2;
    stringSplit(aux[i], aux2, "=");
    if(aux2.size() == 2) {
      if(aux2[0].compare("scan") == 0) {
        scan = getInt(aux2[1].c_str());
        return true;
      }
    }
  }

  // if ID exists, get
  return false;
}
////////////////////////////////////////////////////////////////////////////////
// Convert MSData (PWiz) to specset (SPS)
////////////////////////////////////////////////////////////////////////////////
int convertMSData2SpecSet(MSData &msd, SpecSet &specSet, string & fn)
{
  SpectrumList& spectrumList = *msd.run.spectrumListPtr;
  SpectrumPtr spectrum;
  //Read m/z and intensity values from the spectra
  const bool getBinaryData = true;
  size_t numSpectra = spectrumList.size();

  for (int i = 0 ; i < numSpectra ; ++i) {
    // define specnets spectrum object
    specnets::Spectrum spect;

    ///////// Get the peak list ///////
    // Get spectrum binary data
  	spectrum = spectrumList.spectrum(i, getBinaryData);
  	// Define pairs data holder
  	vector<MZIntensityPair> pairs;
  	// Get MZ intensity pairs
  	spectrum->getMZIntensityPairs(pairs);
    // Add the peaks
  	for (vector<MZIntensityPair>::const_iterator it = pairs.begin(), end = pairs.end(); it!=end; ++it)
      // add a peak
      spect.push_back(TwoValues<float>(it->mz, it->intensity));


    ///////// Data ///////
    int index;

    // MS Level
    if((index = findIndex(spectrum, MS_ms_level)) >= 0)
      spect.msLevel             = getInt(spectrum->cvParams[index].value.c_str());

    // parent M/Z
    if((index = findIndex(spectrum, MS_base_peak_m_z)) >= 0)
      spect.parentMZ            = getFloat(spectrum->cvParams[index].value.c_str());
    else
    if((index = findIndex(spectrum, MS_base_peak)) >= 0)
      spect.parentMZ            = getFloat(spectrum->cvParams[index].value.c_str());


    // parent charge
    if(spect.msLevel == 1)
      spect.parentCharge = 1;
    else
    if((index = findIndex(spectrum, MS_charge_state)) >= 0)
      spect.parentCharge        = getInt(spectrum->cvParams[index].value.c_str());

    // base peak intensity
    if((index = findIndex(spectrum, MS_base_peak_intensity)) >= 0)
      spect.precursor_intensity = getInt(spectrum->cvParams[index].value.c_str());

    // scan number
    int scan;
    if(getScan(spectrum, scan))
      spect.scan = scan;

    // File name
    size_t  found = fn.find_last_of("/\\");
    spect.fileName  = fn.substr(found + 1);

    // parent mass
    if(spect.parentCharge > 0)
      spect.parentMass          =  spect.parentMZ * spect.parentCharge - AAJumps::massHion * (spect.parentCharge - 1.0);
    else
      spect.parentMass          =  spect.parentMZ;

    // MS Fragmentation model
    spect.msFragType            = specnets::Spectrum::FragType_CID;
    if((index = findIndex(spectrum, MS_electron_transfer_dissociation)) >= 0)
      spect.msFragType          = specnets::Spectrum::FragType_ETD;
    if((index = findIndex(spectrum, MS_high_energy_collision_induced_dissociation)) >= 0)
      spect.msFragType          = specnets::Spectrum::FragType_HCD;




    // store the spectrum in the specset
    specSet.push_back(spect);
  }
}
////////////////////////////////////////////////////////////////////////////////
// Convert MSData (PWiz) to specset (SPS)
////////////////////////////////////////////////////////////////////////////////
int convertMSData2psmSet(MSData &msd, PeptideSpectrumMatch &psmSet, string & fn)
{
}
////////////////////////////////////////////////////////////////////////////////
bool loadDataUsingPWiz(string &inName, SpecSet &specs)
{
	FullReaderList readers;
	//populate an MSData object from an MS data filepath
  try {
		MSDataFile msd(inName, &readers);
    convertMSData2SpecSet(msd, specs, inName);
	} catch (exception& e) {
    return false;
  }
  return true;
}
////////////////////////////////////////////////////////////////////////////////
bool loadDataUsingPWiz(string &inName, PeptideSpectrumMatch &psmSet)
{
	FullReaderList readers;
	//populate an MSData object from an MS data filepath
  try {
		MSDataFile msd(inName, &readers);
    convertMSData2psmSet(msd, psmSet, inName);
	} catch (exception& e) {
    return false;
  }
  return true;
}
////////////////////////////////////////////////////////////////////////////////
#else
////////////////////////////////////////////////////////////////////////////////
bool loadDataUsingPWiz(string &inName, SpecSet &specs)
{
  return false;
}
////////////////////////////////////////////////////////////////////////////////
#endif
////////////////////////////////////////////////////////////////////////////////
