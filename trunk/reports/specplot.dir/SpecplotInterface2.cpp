///////////////////////////////////////////////////////////////////////////////
#include "SpecplotInterface2.h"
#include "mzxml.h"
#include "PWizInterface.h"

///////////////////////////////////////////////////////////////////////////////
namespace specnets {
///////////////////////////////////////////////////////////////////////////////
SpecplotInterface2::SpecplotInterface2()
{
}

SpecplotInterface2::~SpecplotInterface2()
{
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int SpecplotInterface2::load(string &fn, string &ext, bool spectrumPresent)
{
  bool fileLoaded = false;
  // load te file
  if(!specSet) {
    // create object
    specSet = new SpecSet();
    // load mzxml file
    if(ext.compare("mzxml") == 0) {
      // auxilizary vector needed
      vector<short> msLevel;
      // hold the return value
      int ret;
      // load method depends on the presence of spectrumscan specification
      try {
        if(m_spectrumScan.empty()) {
          fileLoaded = LoadMzxml( (char * const)(fn.c_str()), *specSet, & msLevel, 2);
        } else {
          fileLoaded = LoadMzxml(fn.c_str(), *specSet, m_spectrumScan.c_str(), & msLevel, 2);
        }
      } catch (...) {
        fileLoaded = false;
      }
    }

    if(!fileLoaded) {

      // load other formats
      fileLoaded = specSet->Load(fn.c_str());

      // load using pwiz
      if(!fileLoaded)
        fileLoaded = loadDataUsingPWiz(fn, *specSet);
    }
  }
  
  if(!fileLoaded) {
    stringstream err;
    err << "Error loading file: " << fn;
    return error(err.str());
    } else if(!spectrumPresent) {
        // get spectrumscan from index
        bool found = false;
        for(int i = 0 ; i < specSet->size() ; i++){
        if((*specSet)[i].scan == getInt(m_spectrumScan.c_str())) {
            m_spectrumIndex = i+1;
            plotSpectrum.setSpectrumIndex(m_spectrumIndex);
            found = true;
        }
        }
        

        if(!found) {
        stringstream err; err << "ERROR: Spectrumscan not found.";
        return error(err.str());
        }
    }
  
  return fileLoaded;
}
///////////////////////////////////////////////////////////////////////////////
}; // namespace specnets
///////////////////////////////////////////////////////////////////////////////
