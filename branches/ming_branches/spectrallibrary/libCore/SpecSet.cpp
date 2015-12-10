/*
 * SpecSet.cpp
 *
 *  Created on: Aug 30, 2011
 *      Author: jsnedecor
 */

//
// **************************************************************
//    SpecSet methods
// **************************************************************
//
#include "SpecSet.h"

namespace specnets
{
  // -------------------------------------------------------------------------
  SpecSet::SpecSet(unsigned int sz)
  {
    specs.resize(sz);
  }
  // -------------------------------------------------------------------------
  SpecSet::SpecSet(char *filename)
  {
    LoadSpecSet_pkl(filename);
  }
  // -------------------------------------------------------------------------
  Spectrum &SpecSet::operator[](unsigned int i)
  {
    return specs[i];
  }
  // -------------------------------------------------------------------------
  const Spectrum &SpecSet::operator[](unsigned int i) const
  {
    return specs[i];
  }
  // -------------------------------------------------------------------------
  void SpecSet::insert(vector<Spectrum>::iterator position,
                       vector<Spectrum>::iterator start,
                       vector<Spectrum>::iterator end)
  {
    specs.insert(position, start, end);
  }
  // -------------------------------------------------------------------------
  SpecSet &SpecSet::operator=(const SpecSet &other)
  {
    specs.resize(other.specs.size());
    for (unsigned int i = 0; i < other.specs.size(); i++)
      specs[i] = other[i];
  }
  // -------------------------------------------------------------------------
  void SpecSet::push_back(const Spectrum & x)
  {
    specs.push_back(x);
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::resize(unsigned int newSize)
  {
    specs.resize(newSize);
    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::size() const
  {
    return specs.size();
  }
  // -------------------------------------------------------------------------
  vector<Spectrum>::iterator SpecSet::begin(void)
  {
    return specs.begin();
  }
  // -------------------------------------------------------------------------
  vector<Spectrum>::iterator SpecSet::end(void)
  {
    return specs.end();
  }

  void SpecSet::setPeakTolerance(float tolerance, bool applyPPM) {
	for (int i = 0; i < size(); i++) {
		if (applyPPM) {
			specs[i].setPeakTolerancePPM(tolerance);
		} else {
			specs[i].setPeakTolerance(tolerance);
		}
	}
  }

  void SpecSet::setParentMassTolerance(float tolerance, bool applyPPM) {
	for (int i = 0; i < size(); i++) {
		if (applyPPM) {
			specs[i].setParentMassTolPPM(tolerance);
		} else {
			specs[i].setParentMassTol(tolerance);
		}
	}
  }


  // -------------------------------------------------------------------------
  void SpecSet::appendSpecSet(SpecSet other)
  {
    unsigned int prevSize = specs.size();
    specs.resize(prevSize + other.specs.size());
    for (unsigned int i = 0; i < other.specs.size(); i++)
    {
      specs[prevSize + i] = other[i];
      specs[prevSize + i].scan = prevSize + other[i].scan;
    }
  }
  // -------------------------------------------------------------------------
  vector<list<int> > &SpecSet::getMassesHistogram(vector<list<int> > &results,
                                                  float resolution)
  {
    int i;
    double t, massIdx = 0; // Value is always integer but max is (double,double)

    for (i = 0; i < specs.size(); i++)
    {
      t = round(specs[i].parentMass / resolution);
      if (massIdx < t)
        massIdx = t;
    }
    results.resize((int)massIdx + 1);

    for (i = 0; i < specs.size(); i++)
      results[(int)(round(specs[i].parentMass / resolution))].push_front(i);
    return results;
  }
  // -------------------------------------------------------------------------
  void SpecSet::setResolution(float newResolution, bool enforceRounding)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].setResolution(newResolution, enforceRounding);
  }
  // -------------------------------------------------------------------------
  void SpecSet::addZPMpeaks(float tolerance, float ionOffset, bool includeY0k)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].addZPMpeaks(tolerance, ionOffset, includeY0k);
  }
  // -------------------------------------------------------------------------
  void SpecSet::normalize(float newTotalInt)
  {
    for (unsigned int i = 0; i < specs.size(); i++)
      specs[i].normalize(newTotalInt);
  }
  // -------------------------------------------------------------------------
  void SpecSet::clearPsms(void)
  {
    for (int i = 0; i < specs.size(); i++)
    {
      specs[i].psmList.clear();
    }
  }
  // -------------------------------------------------------------------------
  template<class T> unsigned int SpecSet::extract(vector<T> &features,
                                                  T featureValue,
                                                  SpecSet &output)
  {
    output.resize(max(features.size(), specs.size()));
    unsigned int keptSpectra = 0, pivot;

    for (pivot = 0; pivot < output.size(); pivot++)
      if (features[pivot] == featureValue)
        output[keptSpectra++] = specs[pivot];

    output.resize(keptSpectra);
    return keptSpectra;
  }
  // -------------------------------------------------------------------------
  void instantiate_SpecSet_extract()
  {
    SpecSet specs, specsOut;

    vector<int> features;
    specs.extract(features, (int)0, specsOut);
  }
  // -------------------------------------------------------------------------
  bool SpecSet::maximumParsimony(void)
  {
    // Initialize assignment vector
    vector<psmPtr> assignments(specs.size());
    for (int pepIdx = 0; pepIdx < specs.size(); pepIdx++)
    {
      psmPtr p(new PeptideSpectrumMatch);
      assignments[pepIdx] = p;
      assignments[pepIdx]->m_dbIndex = -1;
    }

    DEBUG_TRACE;
    int protIdx = 0;
    // Find maximum index of a matched protein
    list<psmPtr>::iterator iter;
    for (int pepIdx = 0; pepIdx < specs.size(); pepIdx++)
    {
      for (iter = specs[pepIdx].psmList.begin(); iter
          != specs[pepIdx].psmList.end(); iter++)
      {
        protIdx = max((int)protIdx, (*iter)->m_dbIndex);
      }
    }
    vector < list<psmPtr> > pepsPerProt(protIdx + 1);
    DEBUG_VAR(pepsPerProt.size());

    unsigned int maxHits, maxHitsIdx;
    float maxHitsMass;
    do
    {
      maxHits = 0;
      maxHitsIdx = 0;

      // Clear the peptides per protein vector
      for (int protIdx = 0; protIdx < pepsPerProt.size(); protIdx++)
        pepsPerProt[protIdx].clear();

      // Create the lists of peptides per protein and compute the max
      for (int pepIdx = 0; pepIdx < specs.size(); pepIdx++)
      {
        for (iter = specs[pepIdx].psmList.begin(); iter
            != specs[pepIdx].psmList.end(); iter++)
        {
          int dbIndex = (*iter)->m_dbIndex;
          pepsPerProt[dbIndex].push_back(*iter);
          if (pepsPerProt[dbIndex].size() > maxHits)
          {
            maxHits = pepsPerProt[dbIndex].size();
            maxHitsIdx = dbIndex;
          }
          //cout << (*iter)->m_scanNum << ", " << dbIndex << endl;
        }
      }
      DEBUG_VAR(maxHits);
      DEBUG_VAR(maxHitsIdx);

      if (maxHits > 0)
      {
        for (iter = pepsPerProt[maxHitsIdx].begin(); iter
            != pepsPerProt[maxHitsIdx].end(); iter++)
        {
          //DEBUG_VAR((*iter)->m_scanNum);
          Spectrum * spectrum = getScan((*iter)->m_scanNum);
          if (spectrum == 0)
          {
            ERROR_MSG("Unable to get scan [" << (*iter)->m_scanNum << "]");
            return false;
          }
          spectrum->psmList.clear();
          assignments[(*iter)->m_scanNum - 1] = (*iter);
        }
      }

    } while (maxHits > 0);

    DEBUG_TRACE;
    // Copy the assignments back to the spectrum
    for (int pepIdx = 0; pepIdx < specs.size(); pepIdx++)
    {
      if (assignments[pepIdx]->m_dbIndex != -1)
      {
        specs[pepIdx].psmList.clear();
        specs[pepIdx].psmList.push_back(assignments[pepIdx]);
      }
    }

    DEBUG_TRACE;
    return true;
  }
  // -------------------------------------------------------------------------
  int SpecSet::saveMatchedProts(const char *filename)
  {
    vector < vector<int> > tempMatchedProts(specs.size());
    for (unsigned int i = 0; i < specs.size(); i++)
    {
      tempMatchedProts[i].resize(3);
      if (specs[i].psmList.size() == 0)
      {
        tempMatchedProts[i][0] = -1;
        tempMatchedProts[i][1] = 0;
        tempMatchedProts[i][2] = 0;
      }
      else
      {
        tempMatchedProts[i][0] = specs[i].psmList.front()->m_dbIndex;
        tempMatchedProts[i][1] = specs[i].psmList.front()->m_numMods;
        tempMatchedProts[i][2] = specs[i].psmList.front()->m_matchOrientation;
      }
    }
    return Save_binArray(filename, tempMatchedProts);
  }
  // -------------------------------------------------------------------------
  bool SpecSet::LoadSpecSet_pkl_mic(const char* filename)
  {
    BufferedLineReader blr;
    TwoValues<float> peak;
    Spectrum loadSpec;
    specs.resize(0);
    if (blr.Load(filename) <= 0 || blr.size() < 2)
    {
      cerr << "ERROR: Not enough lines loaded\n";
      return false;
    }
    bool inspec = false;
    for (unsigned int i = 0; i < blr.size(); i++)
    {

      vector < string > peak_vals;
      if (blr.getline(i) == NULL)
        break;
      const char* next_line = blr.getline(i);
      if (!splitText(next_line, peak_vals, "\t "))
      {
        cerr << "ATTENTION: stopped reading peaks before eof in " << filename
            << "\n";
        break;
      }
      if (peak_vals.size() == 0)
      {
        if (inspec)
        {
          loadSpec.sortPeaks();
          specs.push_back(loadSpec);
        }
        inspec = false;
        continue;
      }
      if (peak_vals.size() < 2)
      {
        cerr << "ATTENTION: found wierd line after header: <" << next_line
            << "> ---- skippping ...\n";
        continue;
      }

      if (!inspec && peak_vals.size() == 3)
      {
        loadSpec.resize(0);
        loadSpec.psmList.resize(0);
        loadSpec.parentMass = (float)atof(peak_vals[0].c_str());
        loadSpec.parentMZ = (double)strtod(peak_vals[0].c_str(), NULL);
        loadSpec.parentCharge = (short)atoi(peak_vals[2].c_str());
        if (loadSpec.parentCharge > 0)
        {
          loadSpec.parentMass = loadSpec.parentMass
              * ((float)loadSpec.parentCharge) - ((float)loadSpec.parentCharge)
              + AAJumps::massHion;
        }
        inspec = true;
        continue;
      }

      if (inspec && peak_vals.size() < 2)
      {
        cerr << "ATTENTION: not enough peak info after header: <" << next_line
            << "> ---- skippping ...\n";
        continue;
      }

      if (inspec)
      {
        peak[0] = (float)atof(peak_vals[0].c_str());
        peak[1] = (float)atof(peak_vals[1].c_str());
        loadSpec.insertPeak(peak[0], peak[1], 0);
      }
    }
    if (inspec)
    {
      loadSpec.sortPeaks();
      specs.push_back(loadSpec);
    }
    return true;
  }
  // -------------------------------------------------------------------------
  Spectrum* SpecSet::getScan(int scan_num)
  {
    //cheating a bit, check to see whether scans are in order and we can just jump to the correct
    //scan
    if (scan_num < specs.size() && scan_num > 0)
    { // only try this if scan number is within range.
      Spectrum *curr_scan = &(specs[scan_num - 1]);

      if (curr_scan->scan == scan_num)
      {
        return curr_scan;
      }
    }

    vector<Spectrum>::iterator spec_iter;
    for (spec_iter = specs.begin(); spec_iter != specs.end(); spec_iter++)
    {
      if (scan_num == spec_iter->scan)
      {
        return &(*spec_iter);
      }
    }
    return (Spectrum*)NULL;
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_mgf(const char *filename)
  {
    BufferedLineReader blr;
    resize(0);
    if (blr.Load(filename) <= 0 or blr.size() < 3)
      return 0; // A valid file must have at least 3 lines: BEGIN_IONS, END_IONS and one mass/intensity peak

    unsigned int lineIdx, specIdx, peakIdx;

    // Counts number of spectra and number of peaks per spectrum
    list<int> peaksPerSpec;
    int first;
    int numPeaks = 0; // Counts number of peaks in the spectrum
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      if (strncmp("END IONS", blr.getline(lineIdx), 8) == 0)
      {
        peaksPerSpec.push_back(numPeaks);
        numPeaks = -1;
        continue;
      }
      if (strncmp("BEGIN IONS", blr.getline(lineIdx), 10) == 0)
      {
        numPeaks = 0;
        continue;
      }
      first = (int)blr.getline(lineIdx)[0];
      if (numPeaks >= 0 and first >= 48 and first <= 57)
        numPeaks++;
    }

    // Parse spectra
    resize(peaksPerSpec.size());
    lineIdx = 0;
    char *token, *line;
    for (specIdx = 0; specIdx < specs.size(); specIdx++)
    {
      // Skip empty lines
      while (lineIdx < blr.size() and (blr.getline(lineIdx)[0] == 0
          or strcmp("BEGIN IONS", blr.getline(lineIdx)) != 0))
        lineIdx++;
      if (lineIdx == blr.size())
      {
        cerr << "Error loading " << filename << " - " << specIdx
            << " spectra instead of " << specs.size() << "?\n";
        resize(0);
        return 0;
      }

      // Start of spectrum
      if (strcmp("BEGIN IONS", blr.getline(lineIdx)) != 0)
      {
        cerr << "ERROR: Expected BEGIN IONS, found '" << blr.getline(lineIdx)
            << "' (line " << lineIdx + 1 << ")\n";
        return 0;
      }
      else
        lineIdx++;

      // Read peaks/charge/parent mass
      specs[specIdx].resize(peaksPerSpec.front());
      peaksPerSpec.pop_front();
      peakIdx = 0;
      specs[specIdx].parentCharge = 0;
      specs[specIdx].parentMassTol = 0;
      specs[specIdx].scan = 0;
      while (lineIdx < blr.size() and strcmp("END IONS", blr.getline(lineIdx))
          != 0)
      {
        line = blr.getline(lineIdx++);
        if (line[0] >= 48 and line[0] <= 57)
        {
          token = strtok(line, " \t");
          if (!token)
          {
            cerr << "Error loading " << filename
                << " - could not parse peak mass on line " << lineIdx << "!\n";
            resize(0);
            return 0;
          }
          specs[specIdx][peakIdx][0] = atof(token);
          token = strtok(NULL, " \t");
          if (!token)
          {
            cerr << "Error loading " << filename
                << " - could not parse peak intensity on line " << lineIdx
                << "!\n";
            resize(0);
            return 0;
          }
          specs[specIdx][peakIdx][1] = atof(token);
          peakIdx++;
          continue;
        }

        if (strncmp("CHARGE=+", line, 8) == 0)
        {
          specs[specIdx].parentCharge = (short)atof(&line[8]);
          continue;
        }
        else if (strncmp("CHARGE=", line, 7) == 0)
        {
          specs[specIdx].parentCharge = (short)atof(&line[7]);
          continue;
        }

        if (strncmp("TITLE=Scan Number: ", line, 19) == 0)
          specs[specIdx].scan = (unsigned int)atof(&line[19]);
        if (strncmp("SCANS=", line, sizeof("SCANS=") - 1) == 0)
          specs[specIdx].scan = (unsigned int)atof(line + sizeof("SCANS=") - 1);
        if (strncmp("PEPMASS=", line, 8) == 0)
        {
          specs[specIdx].parentMass = (double)strtod(&line[8], NULL);
          specs[specIdx].parentMZ = (double)strtod(&line[8], NULL);
        }
        if (strncmp("PRECURSOR=", line, 10) == 0)
        {
          specs[specIdx].parentMZ = (double)strtod(&line[10], NULL);
        }
        if (strncmp("PEPTIDE=", line, sizeof("PEPTIDE=") - 1) == 0)
        {
          psmPtr psm(new PeptideSpectrumMatch());
          psm->m_annotation = line + sizeof("PEPTIDE=") - 1;
          specs[specIdx].psmList.push_back(psm);
        }
        if (strncmp("SEQ=", line, sizeof("SEQ=") - 1) == 0)
        {
          psmPtr psm(new PeptideSpectrumMatch());
          psm->m_annotation = line + sizeof("SEQ=") - 1;
          specs[specIdx].psmList.push_back(psm);
        }
        if (strncmp("MSLEVEL=", line, sizeof("MSLEVEL=") - 1) == 0){
            specs[specIdx].msLevel = atoi(line + sizeof("MSLEVEL=") - 1);
        }
        if (strncmp("SOURCE_INSTRUMENT=", line, sizeof("SOURCE_INSTRUMENT=") - 1) == 0){
            specs[specIdx].instrument_name = line + sizeof("SOURCE_INSTRUMENT=") - 1;
        }
        if (strncmp("ITOL=", line, sizeof("ITOL=") - 1) == 0){
            specs[specIdx].ITOL = atof(line + sizeof("ITOL=") - 1);
        }
        if (strncmp("ITOLU=", line, sizeof("ITOLU=") - 1) == 0){
            specs[specIdx].ITOLU = line + sizeof("ITOLU=") - 1;
        }
        if (strncmp("TOL=", line, sizeof("TOL=") - 1) == 0){
            specs[specIdx].TOL = atof(line + sizeof("TOL=") - 1);
        }
        if (strncmp("TOLU=", line, sizeof("TOLU=") - 1) == 0){
            specs[specIdx].TOLU = line + sizeof("TOLU=") - 1;
        }
        if (strncmp("SPECTRUMQUALITY=", line, sizeof("SPECTRUMQUALITY=") - 1) == 0){
            specs[specIdx].spectrum_quality = atoi(line + sizeof("SPECTRUMQUALITY=") - 1);
        }
        
        if (strncmp("SUBMISSION_USER=", line, sizeof("SUBMISSION_USER=") - 1) == 0){
            if(specs[specIdx].psmList.size() > 0){
                specs[specIdx].psmList.front()->m_submission_user = line + sizeof("SUBMISSION_USER=") - 1;
            }
        }
        if (strncmp("SUBMISSION_ID=", line, sizeof("SUBMISSION_ID=") - 1) == 0){
            cout<<"HERE"<<endl;
            if(specs[specIdx].psmList.size() > 0){
                specs[specIdx].psmList.front()->m_submission_id = line + sizeof("SUBMISSION_ID=") - 1;
            }
        }
        if (strncmp("SUBMISSION_DATE=", line, sizeof("SUBMISSION_DATE=") - 1) == 0){
            if(specs[specIdx].psmList.size() > 0){
                specs[specIdx].psmList.front()->m_submission_date = line + sizeof("SUBMISSION_DATE=") - 1;
            }
        }
        if (strncmp("ORGANISM=", line, sizeof("ORGANISM=") - 1) == 0){
            if(specs[specIdx].psmList.size() > 0){
                specs[specIdx].psmList.front()->m_organism = line + sizeof("ORGANISM=") - 1;
            }
        }
        if (strncmp("COMPOUNDNAME=", line, sizeof("COMPOUNDNAME=") - 1) == 0){
            if(specs[specIdx].psmList.size() > 0){
                specs[specIdx].psmList.front()->m_compound_name = line + sizeof("COMPOUNDNAME=") - 1;
            }
        }
        if (strncmp("MOLECULEMASS=", line, sizeof("MOLECULEMASS=") - 1) == 0){
            if(specs[specIdx].psmList.size() > 0){
                specs[specIdx].psmList.front()->m_molecule_mass = atof(line + sizeof("MOLECULEMASS=") - 1);
            }
        }
        if (strncmp("ACTIVATION=", line, sizeof("ACTIVATION=") - 1) == 0)
        {

            string activ = line + sizeof("ACTIVATION=") - 1;
            std::transform(activ.begin(), activ.end(), activ.begin(), ::toupper);
            if (activ == "CID")
            {
              specs[specIdx].msFragType = Spectrum::FragType_CID;
            }
            else if (activ == "HCD")
            {
              specs[specIdx].msFragType = Spectrum::FragType_HCD;
            }
            else if (activ == "ETD")
            {
              specs[specIdx].msFragType = Spectrum::FragType_ETD;
            }
            else
            {
              WARN_MSG("Unrecognized ACTIVATION field \'" << activ << "\'");
            }

         }
      }
      
      
      if (specs[specIdx].parentCharge > 0)
        specs[specIdx].parentMass = (specs[specIdx].parentMass
            * specs[specIdx].parentCharge) - (AAJumps::massHion
            * (specs[specIdx].parentCharge - 1));
      if (specs[specIdx].scan <= 0)
      {
        specs[specIdx].scan = specIdx + 1;
      }
      lineIdx++;
    }

    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_ms2(const char *filename)
  {
    BufferedLineReader blr;
    resize(0);
    if (blr.Load(filename) <= 0 or blr.size() < 3)
      return 0; // A valid file must have at least 3 lines: 1 header (2 lines) + 1 (m/z,intensity) pair

    unsigned int lineIdx, specIdx, peakIdx;

    // Counts number of spectra and number of peaks per spectrum
    list<int> peaksPerSpec;
    int numPeaks = 1; // Counts number of peaks in the spectrum
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
      if (blr.getline(lineIdx)[0] == ':')
      {
        if (numPeaks > 0)
          peaksPerSpec.push_back(numPeaks);
        numPeaks = -1;
      }
      else if (blr.getline(lineIdx)[0] != 0)
        numPeaks++;
    peaksPerSpec.push_back(numPeaks); // Number of peaks in the last spectrum
    peaksPerSpec.pop_front(); // First element is just the '1' used to initialize numPeaks

    // Parse spectra
    resize(peaksPerSpec.size());
    lineIdx = 0;
    char *token;
    for (specIdx = 0; specIdx < specs.size(); specIdx++)
    {
      // Skip empty lines
      while (blr.getline(lineIdx)[0] == 0 and lineIdx < blr.size())
        lineIdx++;
      if (lineIdx == blr.size())
      {
        cerr << "Error loading " << filename << " - " << specIdx
            << " spectra instead of " << specs.size() << "?\n";
        resize(0);
        return 0;
      }

      // Parse header(s)
      while (blr.getline(lineIdx)[0] == ':')
      {
        lineIdx++;
        token = strtok(blr.getline(lineIdx++), " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse parent mass for spectrum " << specIdx + 1
              << "?\n";
          resize(0);
          return 0;
        }
        specs[specIdx].parentMass = (float)atof(token);
        specs[specIdx].parentMZ = (double)strtod(token, NULL);
        token = strtok(NULL, " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse parent charge for spectrum " << specIdx
              + 1 << "?\n";
          resize(0);
          return 0;
        }
        specs[specIdx].parentCharge = (short)atof(token);
        if (lineIdx == blr.size())
        {
          cerr << "Error loading " << filename << " - " << specIdx
              << " spectra instead of " << specs.size() << "?\n";
          resize(0);
          return 0;
        }
      }

      // Read spectrum peaks
      specs[specIdx].resize(peaksPerSpec.front());
      peaksPerSpec.pop_front();
      for (peakIdx = 0; peakIdx < specs[specIdx].size(); peakIdx++)
      {
        token = strtok(blr.getline(lineIdx++), " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse peak mass for spectrum " << specIdx + 1
              << ", peak " << peakIdx + 1 << "?\n";
          resize(0);
          return 0;
        }
        specs[specIdx][peakIdx][0] = (float)atof(token);
        token = strtok(NULL, " \t");
        if (!token)
        {
          cerr << "Error loading " << filename
              << " - could not parse peak intensity for spectrum " << specIdx
              + 1 << ", peak " << peakIdx + 1 << "?\n";
          resize(0);
          return 0;
        }
        specs[specIdx][peakIdx][1] = (float)atof(token);
        if (lineIdx == blr.size() and peakIdx < specs[specIdx].size() - 1)
        {
          cerr << "Error loading " << filename
              << " - end of file before end of spectrum " << specIdx + 1
              << " ended, got only " << peakIdx + 1 << " peaks instead of "
              << specs.size() << "?\n";
          resize(0);
          return 0;
        }
      }
    }

    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_prms(const char *filename)
  {
    unsigned int specIdx, lineIdx, peakIdx;
    list<unsigned int> specSizes;
    char *line;
    BufferedLineReader blr;

    // Load whole file into memory
    if (blr.Load(filename) <= 0)
      return 0;

    // Count # spectra in the file
    bool inSpectrum = false;
    peakIdx = 0;
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      line = blr.getline(lineIdx);
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      {
        inSpectrum = true;
        peakIdx = 0;
        continue;
      }
      if (line[0] == 0)
      {
        if (inSpectrum)
          specSizes.push_back(peakIdx);
        inSpectrum = false;
      }
      else if (inSpectrum)
        peakIdx++;
    }
    if (inSpectrum)
      specSizes.push_back(peakIdx); // In case there is no empty line after the last spectrum
    specs.resize(specSizes.size());

    // Parse the text
    peakIdx = 0;
    specIdx = 0;
    list<unsigned int>::iterator sizesIter = specSizes.begin();
    char *token;
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      line = blr.getline(lineIdx);

      if (line[0] == 0)
      {
        if (inSpectrum)
        { // End of current spectrum
          if (peakIdx > 0 and peakIdx <= specs[specIdx].size())
            specs[specIdx].parentMass = specs[specIdx][peakIdx - 1][0]
                + AAJumps::massMH;
          specIdx++;
          inSpectrum = false;
        }
        continue;
      }

      // Check for start of a new spectrum
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      {
        specs[specIdx].info = (char *)malloc(strlen(line) - 1);
        strcpy(specs[specIdx].info, &line[2]);
        char *tok = strtok(line, " \t");
        tok = strtok(NULL, " \t"); // Skip ">> " and "<file_index> "
        tok = strtok(NULL, " \t");
        specs[specIdx].scan = (unsigned int)strtoul(tok, NULL, 10);
        specs[specIdx].resize(*sizesIter);
        specs[specIdx].psmList.resize(0);
        sizesIter++;
        peakIdx = 0;
        inSpectrum = true;
        continue;
      }

      // Peak <mass> <intensity> pair
      token = strtok(line, " \t");
      if (token[0] == 0)
      {
        cerr << "ERROR reading peak mass for peak " << peakIdx
            << " for the spectrum entitled (" << specs[specIdx].info
            << ") in file " << filename << "!\n";
        specs.resize(0);
        return 0;
      }
      specs[specIdx][peakIdx][0] = atof(token);
      token = strtok(NULL, " \t");
      if (token[0] == 0)
      {
        cerr << "ERROR reading peak intensity for peak " << peakIdx
            << " for the spectrum entitled (" << specs[specIdx].info
            << ") in file " << filename << "!\n";
        specs.resize(0);
        return 0;
      }
      specs[specIdx][peakIdx][1] = atof(token);
      peakIdx++;
    }

    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_prmsv3(const char *filename)
  {
    unsigned int specIdx, lineIdx, peakIdx;
    list<unsigned int> specSizes;
    char *line;
    BufferedLineReader blr;

    // Load whole file into memory
    if (blr.Load(filename) <= 0)
      return 0;

    // Count # spectra in the file
    bool inSpectrum = false;
    peakIdx = 0;
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    {
      line = blr.getline(lineIdx);
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      {
        inSpectrum = true;
        peakIdx = 0;
        continue;
      }
      if (line[0] != 0 and line[0] == '#' and inSpectrum)
        inSpectrum = false;
      if (line[0] == 0)
      {
        if (inSpectrum)
          specSizes.push_back(peakIdx);
        inSpectrum = false;
      }
      else if (inSpectrum and line[0] != 'C')
        peakIdx++;
    }
    if (inSpectrum)
      specSizes.push_back(peakIdx); // In case there is no empty line after the last spectrum
    specs.resize(specSizes.size());

    // Parse the text
    peakIdx = 0;
    specIdx = 0;
    list<unsigned int>::iterator sizesIter = specSizes.begin();
    char *token;
    for (lineIdx = 0; lineIdx < blr.size(); lineIdx++)
    { // Skip header lines
      line = blr.getline(lineIdx);
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
        break;
    }

    for (; lineIdx < blr.size(); lineIdx++)
    {
      line = blr.getline(lineIdx);

      if (line[0] == 0)
      {
        if (inSpectrum)
        { // End of current spectrum
          //  -- Not valid for Pepnovo v3
          //        if(peakIdx>0 and peakIdx<=specs[specIdx].size()) specs[specIdx].parentMass = specs[specIdx][peakIdx-1][0]+AAJumps::massMH;
          specIdx++;
          inSpectrum = false;
        }
        continue;
      }

      // Check for start of a new spectrum
      if (line[0] != 0 and line[0] == '>' and line[1] != 0 and line[1] == '>')
      {
        specs[specIdx].info = (char *)malloc(strlen(line) - 1);
        strcpy(specs[specIdx].info, &line[2]);
        char *tok = strtok(line, " \t");
        tok = strtok(NULL, " \t"); // Skip ">> " and "<file_index> "
        tok = strtok(NULL, " \t");
        specs[specIdx].scan = (unsigned int)strtoul(tok, NULL, 10);

        // Read charge/mass header
        line = blr.getline(++lineIdx);
        if (line[0] != 0 and line[0] == '#')
        {
          inSpectrum = false;
          continue;
        }
        specs[specIdx].parentCharge = (short)atoi(&line[8]);
        specs[specIdx].parentMass = (float)atof(&line[19]);
        peakIdx = 0;
        inSpectrum = true;

        specs[specIdx].resize(*sizesIter);
        specs[specIdx].psmList.resize(0);
        sizesIter++;
        continue;
      }

      // Peak <mass> <intensity> pair
      token = strtok(line, " \t");
      if (token[0] == 0)
      {
        cerr << "ERROR reading peak mass for peak " << peakIdx
            << " for the spectrum entitled (" << specs[specIdx].info
            << ") in file " << filename << "!\n";
        specs.resize(0);
        return 0;
      }
      specs[specIdx][peakIdx][0] = atof(token);
      token = strtok(NULL, " \t");
      if (token[0] == 0)
      {
        cerr << "ERROR reading peak intensity for peak " << peakIdx
            << " for the spectrum entitled (" << specs[specIdx].info
            << ") in file " << filename << "!\n";
        specs.resize(0);
        return 0;
      }
      specs[specIdx][peakIdx][1] = atof(token);
      peakIdx++;
    }

    return specs.size();
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_pklbin_with_annotation(const char * spectra_filename,
                                                           const char * annotation_filename)
  {
    SpecSet temp_specs;
    temp_specs.LoadSpecSet_pklbin(spectra_filename);

    PeptideSpectrumMatchSetSpectralLibraryLoader psm_set;
    psm_set.loadSpecnetsResultsFile(annotation_filename);

    int original_specs_size = specs.size();
    cout << "Original Size: " << original_specs_size << endl;
    cout << temp_specs.size() << endl;
    specs.insert(specs.end(), temp_specs.specs.begin(), temp_specs.specs.end());

    for (int psm_idx = 0; psm_idx < psm_set.size(); psm_idx++)
    {
      psm_set[psm_idx]->m_spectrum = &(specs[psm_set[psm_idx]->m_scanNum - 1
          + original_specs_size]);
      psm_set[psm_idx]->m_spectrum->psmList.push_back(psm_set[psm_idx]);
    }

    return 0;
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_mgf_with_annotation(const char * spectra_filename,
                                                        const char * annotation_filename)
  {
    SpecSet temp_specs;
    temp_specs.LoadSpecSet_mgf(spectra_filename);

    PeptideSpectrumMatchSetSpectralLibraryLoader psm_set;
    psm_set.loadSpecnetsResultsFile(annotation_filename);

    int original_specs_size = specs.size();
    cout << "Original Size: " << original_specs_size << endl;
    cout << temp_specs.size() << endl;
    specs.insert(specs.end(), temp_specs.specs.begin(), temp_specs.specs.end());

    for (int psm_idx = 0; psm_idx < psm_set.size(); psm_idx++)
    {
      psm_set[psm_idx]->m_spectrum = &(specs[psm_set[psm_idx]->m_scanNum - 1
          + original_specs_size]);
      psm_set[psm_idx]->m_spectrum->psmList.push_back(psm_set[psm_idx]);
    }

    return 0;
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_pkl(const char *filename)
  {
    int numSpecs = 1; // Assume file contains at least one spectrum and count one additional spectrum per blank line
    ifstream input(filename);
    if (!input.is_open() || !input.good())
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      specs.resize(0);
      return 0;
    }

    list<int> specSizes; // Register number of peaks per spectrum

    char *lineBuffer = (char *)malloc(1025);
    //    char *lineBuffer = new char [1024];
    int numPeaks = -1; // zero only after reading first tuple of (parent mass, intensity, charge)
    while (!input.eof() && !input.fail())
    {
      input.getline(lineBuffer, 1024, '\n');
      if (lineBuffer[0] == 0 or lineBuffer[0] == '\n' or lineBuffer[0] == '\r')
      {
        if (numPeaks >= 0)
        {
          numSpecs++;
          specSizes.push_back(numPeaks);
          numPeaks = -1;
        }
      }
      else
      {
        if (((int)lineBuffer[0] >= (int)'0' and (int)lineBuffer[0] <= (int)'9')
            or lineBuffer[0] == '-')
          numPeaks++;
      }
    }
    if (numPeaks >= 0)
    {
      specSizes.push_back(numPeaks);
    }
    else
    {
      numSpecs--; // In case there was no empty line at the end of the pkl file
    }

    specs.resize(numSpecs);
    //    input.seekg(0,ios_base::beg);
    //    input.close();  input.open(filename);
    //    Neither of the above worked, so ...
    input.close();
    ifstream input2(filename);

    float foo;
    for (int i = 0; i < numSpecs; i++)
    {
      input2 >> specs[i].parentMass >> foo >> specs[i].parentCharge; // Read precursor mass, total intensity, charge state
      if (specs[i].parentCharge > 0)
        specs[i].parentMass = specs[i].parentMass * specs[i].parentCharge
            - specs[i].parentCharge + AAJumps::massHion; // Convert to PM+19

      //        specs[i].numPeaks = specSizes.front();   specSizes.pop_front();
      //        specs[i].peakList.resize(specs[i].numPeaks);
      specs[i].resize(specSizes.front());
      specs[i].psmList.resize(0);
      specSizes.pop_front();

      for (int j = 0; j < specs[i].size(); j++)
        input2 >> specs[i][j][0] >> specs[i][j][1];
      //            input2 >> specs[i].peakList[j][0] >> specs[i].peakList[j][1];
      specs[i].scan = i + 1;
      input2.getline(lineBuffer, 1024, '\n'); // Read intermediate newline
    }

    //    delete lineBuffer;
    input2.close();
    free(lineBuffer);
    return 1;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_info(const char *filename)
  {
    ofstream output(filename);
    if (!output)
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      return -1;
    }

    for (int i = 0; i < specs.size(); i++)
      output << specs[i].info << "\n";

    output.close();
    return 1;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_pkl(const char *filename)
  {
    ofstream output(filename);
    if (!output)
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      return -1;
    }

    for (int i = 0; i < specs.size(); i++)
    {
      output << specs[i].parentMass << " -1 " << specs[i].parentCharge << "\n";
      for (int j = 0; j < specs[i].size(); j++)
      {
        output << specs[i][j][0] << " " << specs[i][j][1] << endl;
        //            output << specs[i].peakList[j][0] << " " << specs[i].peakList[j][1] << endl;
      }
      if (i < specs.size() - 1)
      {
        output << endl;
      }
    }

    output.close();
    return 1;
  }
  // -------------------------------------------------------------------------
  bool SpecSet::SaveAnnotations(const char* filename)
  {
    FILE* out_buf = fopen(filename, "w");

    if (out_buf == NULL)
      return false;

    for (int i = 0; i < size(); i++)
    {
      if (specs[i].psmList.size() == 0)
      {
        fprintf(out_buf, "\n");
      }
      else
      {
        list<psmPtr>::iterator it;
        for (it = specs[i].psmList.begin(); it != specs[i].psmList.end(); it++)
        {
          fprintf(out_buf,
                  "%s,",
                  (*it)->m_annotation.substr(2, (*it)->m_annotation.length()
                      - 4).c_str());
        }
        fprintf(out_buf, "\n");
      }
    }
    fclose(out_buf);
    return true;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_ms2(const char* filename)
  {
    for (int i = 0; i < specs.size(); i++)
    {
      ostringstream outs;
      outs << filename << i << ".ms2";
      ofstream output(outs.str().c_str());
      if (!output)
      {
        cerr << "ERROR: cannot open " << filename << "\n";
        return -1;
      }
      specs[i].output_ms2(output);
      output.close();
    }
    return 1;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_mgf(const char* filename, const char* activation)
  {
    FILE* output = fopen(filename, "w");
    if (!output)
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      return -1;
    }

    float mz, z, pepmass;
    for (int i = 0; i < specs.size(); i++)
    {
      z = (float)specs[i].parentCharge;
      mz = specs[i].parentMZ;
      pepmass = specs[i].parentMass;
      if (z > 0 && mz == 0)
        mz = (pepmass + (z - 1) * AAJumps::massHion) / z;
      fprintf(output, "BEGIN IONS\nPEPMASS=%.5f\nCHARGE=%.0f+\n", mz, z);
      //output << "BEGIN IONS\nPEPMASS="<<mz<<"\nCHARGE=+"<<z<<"\n";
      
      
      if (specs[i].msLevel > 0)
          fprintf(output, "MSLEVEL=%d\n", specs[i].msLevel);
      if (specs[i].instrument_name.length() > 0)
          fprintf(output, "SOURCE_INSTRUMENT=%s\n", specs[i].instrument_name.c_str());
      if (specs[i].ITOL > 0)
          fprintf(output, "ITOL=%f\n", specs[i].ITOL);
      if (specs[i].ITOLU.length() > 0)
          fprintf(output, "ITOLU=%s\n", specs[i].ITOLU.c_str());
      if (specs[i].TOL > 0)
          fprintf(output, "TOL=%f\n", specs[i].TOL);
      if (specs[i].TOLU.length() > 0)
          fprintf(output, "TOLU=%s\n", specs[i].TOLU.c_str());
      if (specs[i].spectrum_quality > 0)
          fprintf(output, "SPECTRUMQUALITY=%d\n", specs[i].spectrum_quality);
      
      
      if (specs[i].psmList.size() > 0)
      {
        list<psmPtr>::iterator it;
        for (it = specs[i].psmList.begin(); it != specs[i].psmList.end(); it++)
        {
          fprintf(output, "SEQ=%s\n", (*it)->m_annotation.c_str());
          
          if((*it)->m_submission_user.length() > 0)
              fprintf(output, "SUBMISSION_USER=%s\n", (*it)->m_submission_user.c_str());
          if((*it)->m_submission_id.length() > 0)
              fprintf(output, "SUBMISSION_ID=%s\n", (*it)->m_submission_id.c_str());
          if((*it)->m_submission_date.length() > 0)
              fprintf(output, "SUBMISSION_DATE=%s\n", (*it)->m_submission_date.c_str());
          if((*it)->m_organism.length() > 0)
              fprintf(output, "ORGANISM=%s\n", (*it)->m_organism.c_str());
          
          if((*it)->m_compound_name.length() > 0)
              fprintf(output, "COMPOUNDNAME=%s\n", (*it)->m_compound_name.c_str());
          if((*it)->m_molecule_mass> 0)
              fprintf(output, "MOLECULEMASS=%f\n", (*it)->m_molecule_mass);
          
        }
      }

      if (activation != 0)
      {
        fprintf(output, "ACTIVATION=%s\n", activation);
        //output<<"ACTIVATION="<<activation<<"\n";
      }
      else
      {
        if (specs[i].msFragType == Spectrum::FragType_CID)
        {
          fprintf(output, "ACTIVATION=CID\n");
        }
        else if (specs[i].msFragType == Spectrum::FragType_HCD)
        {
          fprintf(output, "ACTIVATION=HCD\n");
        }
        else if (specs[i].msFragType == Spectrum::FragType_ETD)
        {
          fprintf(output, "ACTIVATION=ETD\n");
        }
      }
      if (specs[i].scan > 0)
      {
        fprintf(output, "TITLE=Scan Number: %d\n", specs[i].scan);
        //output<<"TITLE=Scan Number: "<<specs[i].scan<<"\n";
      }
      
          
          
      for (int j = 0; j < specs[i].size(); j++)
      {
        fprintf(output, "%.6f %.6f\n", specs[i][j][0], specs[i][j][1]);
        //output << specs[i][j][0] << " " << specs[i][j][1] << endl;
        //            output << specs[i].peakList[j][0] << " " << specs[i].peakList[j][1] << endl;
      }
      fprintf(output, "END IONS\n");
      //output << "END IONS\n";
      if (i < specs.size() - 1)
      {
        fprintf(output, "\n");
        //output << endl;
      }
    }
    fclose(output);
    //output.close();
    return 1;
  }
  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_bin(const char *filename)
  {
    WARN_MSG("SaveSpecSet_bin has been DEPRECATED!");
    WARN_MSG("Data previously saved in bin file is now save in pklbin");
    WARN_MSG("Use savePklBin() and loadPklBin().");
    // Matrix to store data
    vector < vector<int> > data;
    // get # of spectra
    unsigned int numSpecs = specs.size();
    // cycle thru all spectra
    for (int i = 0; i < numSpecs; i++)
    {
      // auxiliary vector to store a pair (scan #, msLevel)
      vector<int> aux;
      // store scan #
      aux.push_back(specs[i].scan);
      // store msLevel
      aux.push_back(specs[i].msLevel);
      // add spectrum data to matrix
      data.push_back(aux);
    }
    // call procedure to save .bin file and return it's return value
    return Save_binArray(filename, data);
  }
  //-------------------------------------------------------------------------
  //  save - saves a SpecSet to a file in binary format.
  int SpecSet::savePklBin(const char * spectrumFilename,
                          const char * psmFilename /* = 0x0 */,
                          const char * peaksFilename /* = 0x0 */)
  {
    FILE *fp;
    unsigned int numSpecs = specs.size();
    unsigned int i, p;

    fp = fopen(spectrumFilename, "w");
    if (fp == 0) {
      cerr << "ERROR: cannot open " << spectrumFilename << "\n";
      return -1;
    }

    unsigned int dummy = 0;
    unsigned int count = fwrite(&dummy, sizeof(int), 1, fp); // Empty 4 bytes

    // Latest version
    char version = 2;
    count = fwrite(&version, sizeof(char), 1, fp); // Version of file
    char subversion = 0;
    count = fwrite(&subversion, sizeof(char), 1, fp); // Sub-version of file

    count = fwrite(&numSpecs, sizeof(int), 1, fp); // Number of spectra in the file


    // Write the scan numbers
    unsigned int *scanNums = (unsigned int *) malloc(sizeof(unsigned int)
        * numSpecs);
    for (unsigned int i = 0; i < numSpecs; i++) {
      scanNums[i] = specs[i].scan;
    }
    count = fwrite(scanNums, sizeof(unsigned int), numSpecs, fp);
    free(scanNums);

    // Write the MS levels
    short *msLevels = (short *) malloc(sizeof(short) * numSpecs);
    for (unsigned int i = 0; i < numSpecs; i++) {
      msLevels[i] = specs[i].msLevel;
    }
    count = fwrite(msLevels, sizeof(short), numSpecs, fp);
    free(msLevels);

    // New in Version 2: Write the fragmentation types
    short *fragTypes = (short *) malloc(sizeof(short) * numSpecs);
    for (unsigned int i = 0; i < numSpecs; i++) {
      fragTypes[i] = (short) specs[i].msFragType;
    }
    count = fwrite(fragTypes, sizeof(short), numSpecs, fp);
    free(fragTypes);

    // New in Version 2: Write the parent masses here
    float *PMs = (float *) malloc(sizeof(float) * numSpecs);
    for (unsigned int i = 0; i < numSpecs; i++) {
      PMs[i] = specs[i].parentMass;
    }
    count = fwrite(PMs, sizeof(float), numSpecs, fp);
    free(PMs);

    // New in Version 2: Write the precursor charges here
    short *precCharges = (short *) malloc(sizeof(short) * numSpecs);
    for (unsigned int i = 0; i < numSpecs; i++) {
      precCharges[i] = specs[i].parentCharge;
    }
    count = fwrite(precCharges, sizeof(short), numSpecs, fp);
    free(precCharges);

    // New in Version 2: Write the parent mass tolerances
    float *pmTols = (float *) malloc(sizeof(float) * numSpecs);
    for (unsigned int i = 0; i < numSpecs; i++) {
      pmTols[i] = specs[i].parentMassTol;
    }
    count = fwrite(pmTols, sizeof(float), numSpecs, fp);
    free(pmTols);

    unsigned short *numPeaks = (unsigned short *) malloc(sizeof(short)
        * numSpecs);
    unsigned short maxNumPeaks = 0;
    for (i = 0; i < numSpecs; i++) {
      numPeaks[i] = specs[i].size();
      maxNumPeaks = max(maxNumPeaks, numPeaks[i]);
    }
    fwrite(numPeaks, sizeof(short), numSpecs, fp); // Number of peaks per spectrum in the file

    // New in Version 2: Write the tolerance of each peak
    //                   Also, parent masses and charges were written earlier, so they are not written here any more
    float *peaksBuffer = (float *) malloc(3 * maxNumPeaks * sizeof(float));
    unsigned int pbIdx;
    for (i = 0; i < numSpecs; i++) {
      for (pbIdx = 0, p = 0; p < numPeaks[i]; p++) {
        peaksBuffer[pbIdx++] = specs[i][p][0];
        peaksBuffer[pbIdx++] = specs[i][p][1];
        peaksBuffer[pbIdx++] = specs[i].getTolerance(p);
      }
      fwrite(peaksBuffer, sizeof(float), 3 * numPeaks[i], fp); // [parentMass charge] followed by [masses intensities]
    }

    free(peaksBuffer);
    free(numPeaks);
    fclose(fp);

    if (psmFilename != 0x0) {
      PeptideSpectrumMatchSet psmSetTemp;
      psmSetTemp.getPSMSet(this);
      psmSetTemp.saveToFile(psmFilename);
    }

    if (peaksFilename) {
      SpecSet tempMatchedPeaks(specs.size());
      for (int i = 0; i < specs.size(); i++) {
        list<psmPtr>::iterator litr = specs[i].psmList.begin();
        int peakListSize = (*litr)->m_matchedPeaks.size();
        tempMatchedPeaks[i].resize(peakListSize);
        for (int j = 0; j < peakListSize; j++) {
          tempMatchedPeaks[i][j].set((*litr)->m_matchedPeaks[j][0],
                                     (*litr)->m_matchedPeaks[j][1]);
        }
      }
      tempMatchedPeaks.SaveSpecSet_pklbin(peaksFilename);
    }

    return 1;
  }
  bool SpecSet::loadPklBin_1(FILE* fp,
                             int numSpecs,
                             bool oldVersion,
                             char subversion)
  {
    float *data; // Pointer to array containing spectrum data
    unsigned int i, p, dataIdx, count, numValues;
    unsigned short *numPeaks; // Pointer to array containing 1+number of peaks per spectrum
    unsigned int dataSize = 0; // Size of the data buffer

    specs.resize(numSpecs);

    // Read in the scan numbers and MS Levels if this is a second generation file
    // First generation files didn't have the dummy field nor this data
    if (!oldVersion) {

      // Read the scan numbers
      unsigned int *scanNums = (unsigned int *) malloc(sizeof(unsigned int)
          * numSpecs);
      if (scanNums == (unsigned int *) 0) {
        ERROR_MSG("Not enough memory for " << numSpecs << " scan numbers");
        fclose(fp);
        return false;
      }

      count = fread(scanNums, sizeof(unsigned int), numSpecs, fp);
      if (count != numSpecs) {
        free(scanNums);
        return false;
      }
      for (unsigned int i = 0; i < numSpecs; i++) {
        specs[i].scan = scanNums[i];
        if (specs[i].scan <= 0) {
          specs[i].scan = i + 1;
        }
      }
      free(scanNums);
      // Read the MS levels
      short *msLevels = (short *) malloc(sizeof(short) * numSpecs);
      if (msLevels == (short *) 0) {
        ERROR_MSG("Not enough memory for " << numSpecs << " MS Levels.");
        fclose(fp);
        return false;
      }
      count = fread(msLevels, sizeof(short), numSpecs, fp);
      if (count != numSpecs) {
        free(msLevels);
        fclose(fp);
        return false;
      }
      for (unsigned int i = 0; i < numSpecs; i++) {
        specs[i].msLevel = msLevels[i];
      }
      free(msLevels);
    }
    else {
      for (unsigned int i = 0; i < numSpecs; i++) {
        specs[i].scan = i + 1;
      }
    }

    // Number of peaks per spectrum
    numPeaks = (unsigned short *) malloc(sizeof(unsigned short) * numSpecs);
    if (numPeaks == (unsigned short *) 0) {
      ERROR_MSG("Not enough memory for " << numSpecs << " spectra");
      fclose(fp);
      return false;
    }
    count = fread(numPeaks, sizeof(unsigned short), numSpecs, fp);
    if (count != numSpecs) {
      free(numPeaks);
      fclose(fp);
      return false;
    }
    for (i = 0; i < numSpecs; i++) {
      if (numPeaks[i] > dataSize)
        dataSize = numPeaks[i];

      //cout << "numPeaks[" << i << "] = " << numPeaks[i] << endl;
    }
    dataSize = 2 * dataSize + 2; // Each spectrum 'line' has 2 values and there's one additional 'line' with parent mass/charge
    data = (float *) malloc(sizeof(float) * dataSize);
    if (data == (float *) 0) {
      ERROR_MSG("Not enough memory for " << dataSize << " floats");
      fclose(fp);
      return false;
    }

    for (i = 0; i < numSpecs; i++) {
      unsigned int numValues = 2 * numPeaks[i] + 2;
      count = fread(data, sizeof(float), numValues, fp);
      if (count != numValues) {
        ERROR_MSG("Not enough memory for " << numSpecs << " numValues");
        free(numPeaks);
        free(data);
        fclose(fp);
        return false;
      }

      specs[i].parentMass = data[0];
      specs[i].parentCharge = (int) data[1];

      if (specs[i].parentCharge > 0) {
        specs[i].parentMZ = (specs[i].parentMass + ((specs[i].parentCharge
            - 1.0) * AAJumps::massHion)) / specs[i].parentCharge;
      }
      else {
        specs[i].parentMZ = specs[i].parentMass;
      }

      //specs[i].scan = i + 1;
      specs[i].resize(numPeaks[i]);
      specs[i].psmList.resize(0);
      for (unsigned int p = 0, dataIdx = 2; dataIdx < numValues; dataIdx += 2, p++)
        specs[i][p].set(data[dataIdx], data[dataIdx + 1]);
    }

    free(numPeaks);
    free(data);
    fclose(fp);

    return true;
  }

  bool SpecSet::loadPklBin_2(FILE* fp, int numSpecs, char subversion)
  {
    float *data; // Pointer to array containing spectrum data
    unsigned int i, p, dataIdx, count, numValues;
    unsigned short *numPeaks; // Pointer to array containing 1+number of peaks per spectrum
    unsigned int dataSize = 0; // Size of the data buffer
    float *PMTols, *PMs;
    short *msLevels, *fragTypes, *precCharges;

    specs.resize(numSpecs);

    // Read in the scan numbers and MS Levels if this is a second generation file
    // First generation files didn't have the dummy field nor this data

    // Read the scan numbers
    unsigned int *scanNums = (unsigned int *) malloc(sizeof(unsigned int)
        * numSpecs);
    if (scanNums == (unsigned int *) 0) {
      ERROR_MSG("Not enough memory for " << numSpecs << " scan numbers.");
      goto load_fail;
    }

    count = fread(scanNums, sizeof(unsigned int), numSpecs, fp);
    if (count != numSpecs) {
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++) {
      specs[i].scan = scanNums[i];
      if (specs[i].scan <= 0) {
        specs[i].scan = i + 1;
      }
    }
    free(scanNums);

    // Read the MS levels
    msLevels = (short *) malloc(sizeof(short) * numSpecs);
    if (msLevels == (short *) 0) {
      ERROR_MSG("Not enough memory for " << numSpecs << " MS Levels.");
      goto load_fail;
    }
    count = fread(msLevels, sizeof(short), numSpecs, fp);
    if (count != numSpecs) {
      free(msLevels);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++) {
      specs[i].msLevel = msLevels[i];
    }
    free(msLevels);

    // Read the fragmentation types
    fragTypes = (short *) malloc(sizeof(short) * numSpecs);
    if (fragTypes == (short *) 0) {
      ERROR_MSG("Not enough memory for " << numSpecs << " MS frag types.");
      goto load_fail;
    }
    count = fread(fragTypes, sizeof(short), numSpecs, fp);
    if (count != numSpecs) {
      free(fragTypes);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++) {
      if (fragTypes[i] == (short) Spectrum::FragType_CID) {
        specs[i].msFragType = Spectrum::FragType_CID;
      }
      else if (fragTypes[i] == (short) Spectrum::FragType_HCD) {
        specs[i].msFragType = Spectrum::FragType_HCD;
      }
      else if (fragTypes[i] == (short) Spectrum::FragType_ETD) {
        specs[i].msFragType = Spectrum::FragType_ETD;
      }
      else {
        ERROR_MSG("Found unsupported fragmentation ID " << fragTypes[i]);
        free(fragTypes);
        goto load_fail;
      }
    }
    free(fragTypes);

    // Read the parent masses
    PMs = (float *) malloc(sizeof(float) * numSpecs);
    if (PMs == (float *) 0) {
      ERROR_MSG("Not enough memory for " << numSpecs << " parent masses.");
      goto load_fail;
    }
    count = fread(PMs, sizeof(float), numSpecs, fp);
    if (count != numSpecs) {
      free(PMs);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++) {
      specs[i].parentMass = PMs[i];
    }
    free(PMs);

    // Read the MS precursor charges
    precCharges = (short *) malloc(sizeof(short) * numSpecs);
    if (precCharges == (short *) 0) {
      ERROR_MSG("Not enough memory for " << numSpecs << " MS Charges.");
      goto load_fail;
    }
    count = fread(precCharges, sizeof(short), numSpecs, fp);
    if (count != numSpecs) {
      free(precCharges);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++) {
      specs[i].parentCharge = precCharges[i];
      if (specs[i].parentCharge > 0) {
        specs[i].parentMZ = (specs[i].parentMass + ((specs[i].parentCharge
            - 1.0) * AAJumps::massHion)) / specs[i].parentCharge;
      }
      else {
        specs[i].parentMZ = specs[i].parentMass;
      }
    }
    free(precCharges);

    // Read the parent masses tolerances
    PMTols = (float *) malloc(sizeof(float) * numSpecs);
    if (PMTols == (float *) 0) {
      ERROR_MSG("Not enough memory for " << numSpecs
          << " parent mass tolerances.");
      goto load_fail;
    }
    count = fread(PMTols, sizeof(float), numSpecs, fp);
    if (count != numSpecs) {
      free(PMTols);
      goto load_fail;
    }
    for (unsigned int i = 0; i < numSpecs; i++) {
      specs[i].parentMassTol = PMTols[i];
    }
    free(PMTols);

    // Number of peaks per spectrum
    numPeaks = (unsigned short *) malloc(sizeof(unsigned short) * numSpecs);
    if (numPeaks == (unsigned short *) 0) {
      ERROR_MSG("Not enough memory for " << numSpecs << " spectra");
      goto load_fail;
    }
    count = fread(numPeaks, sizeof(unsigned short), numSpecs, fp);
    if (count != numSpecs) {
      free(numPeaks);
      goto load_fail;
    }
    for (i = 0; i < numSpecs; i++) {
      if (numPeaks[i] > dataSize) {
        dataSize = numPeaks[i];
      }

      //cout << "numPeaks[" << i << "] = " << numPeaks[i] << endl;
    }
    dataSize = 3 * dataSize; // Each spectrum 'line' has 3 values
    data = (float *) malloc(sizeof(float) * dataSize);
    if (data == (float *) 0) {
      ERROR_MSG("Not enough memory for " << dataSize << " floats.");
      goto load_fail;
    }

    for (i = 0; i < numSpecs; i++) {
      unsigned int numValues = 3 * numPeaks[i];
      count = fread(data, sizeof(float), numValues, fp);
      if (count != numValues) {
        free(numPeaks);
        free(data);
        goto load_fail;
      }

      specs[i].resize(numPeaks[i]);
      specs[i].psmList.resize(0);
      for (unsigned int p = 0, dataIdx = 0; dataIdx < numValues; dataIdx += 3, p++) {
        specs[i][p].set(data[dataIdx], data[dataIdx + 1]);
        specs[i].setTolerance(p, data[dataIdx + 2]);
      }
    }

    free(numPeaks);
    free(data);
    fclose(fp);

    return true;

    load_fail:

    fclose(fp);

    return false;

  }

  //-------------------------------------------------------------------------
  //  load - loads a SpecSet from a binary format file.
  int SpecSet::loadPklBin(const char * filename,
                          const char * psmFilename /* = 0x0 */,
                          const char * peaksFilename /* = 0x0 */)
  {
    FILE *fp;
    unsigned int numSpecs = 0;
    unsigned int i, p, dataIdx, count, numValues;
    char version = 0;
    char subversion = 0;

    fp = fopen(filename, "r");
    if (fp == 0) {
      ERROR_MSG("Opening " << filename);
      return 0;
    }

    count = fread(&numSpecs, sizeof(unsigned int), 1, fp); // Number of spectra in the file
    if (count != 1) {
      ERROR_MSG("Reading number of spectra from " << filename);
      fclose(fp);
      return 0;
    }

    // Read header information, including which version the file was saved in
    bool oldVersion = true;
    if (numSpecs != 0) {
      WARN_MSG("PKLBIN file non-zero first value.");
      WARN_MSG("Assuming first value is number of spectra [" << numSpecs << "]");
    }
    else {
      oldVersion = false;
      version = 0;
      count = fread(&version, sizeof(char), 1, fp); // Version of file
      if (count != 1) {
        ERROR_MSG("Reading version number from " << filename);
        fclose(fp);
        return 0;
      }
      subversion = 0;
      count = fread(&subversion, sizeof(char), 1, fp); // Sub-version of file
      if (count != 1) {
        ERROR_MSG("Reading sub-version number from " << filename);
        fclose(fp);
        return 0;
      }
      count = fread(&numSpecs, sizeof(unsigned int), 1, fp); // Number of spectra in the file
      if (count != 1) {
        ERROR_MSG("Reading number of spectra from " << filename);
        fclose(fp);
        return 0;
      }
    }

    // Load the SpecSet and Spectrum data fields here
    // Each version has its own load method that should be called here
    if (oldVersion or version == 1) {
      // Load version 1
      if (!loadPklBin_1(fp, numSpecs, oldVersion, subversion)) {
        ERROR_MSG("Reading " << filename);
        return 0;
      }
    }
    else if (version == 2) {
      // Load version 2
      if (!loadPklBin_2(fp, numSpecs, subversion)) {
        ERROR_MSG("Reading " << filename);
        return 0;
      }
    }
    else {
      ERROR_MSG("Found unsupported version " << (int) version);
      ERROR_MSG("Reading " << filename);
      return 0;
    }
    /*
     if (countSpectraOnly)
     return numSpecs;
     */

    if (psmFilename) {
      PeptideSpectrumMatchSet psmSetTemp;
      psmSetTemp.loadFromFile(psmFilename);
      DEBUG_VAR(psmSetTemp.size());
      psmSetTemp.addSpectra(this, true);
    }

    if (peaksFilename) {
      DEBUG_TRACE;
      SpecSet tempMatchedPeaks;
      if (tempMatchedPeaks.LoadSpecSet_pklbin(peaksFilename) <= 0) {
        ERROR_MSG("Problem loading matched peaks from [" << peaksFilename
            << "]");
        return 0;
      }

      DEBUG_VAR(tempMatchedPeaks.size());
      for (int i = 0; i < tempMatchedPeaks.size(); i++) {
        if (specs[i].psmList.size() == 0) {
          continue;
        }
        specs[i].psmList.front()->m_matchedPeaks.resize(tempMatchedPeaks[i].size());
        for (int j = 0; j < tempMatchedPeaks[i].size(); j++) {
          specs[i].psmList.front()->m_matchedPeaks[j].set(tempMatchedPeaks[i][j][0],
                                                          tempMatchedPeaks[i][j][1]);
        }
      }
      DEBUG_TRACE;
    }

    return specs.size();
  }

  // -------------------------------------------------------------------------
  short SpecSet::SaveSpecSet_pklbin(const char *filename,
                                    const char *binFilename)
  {
    if (binFilename)
    {
      WARN_MSG("Bin filename [" << binFilename << "] passed to SaveSpecSet_pklbin");
      WARN_MSG("Bin file [" << binFilename << "] will NOT be written.");
    }
    return (short)savePklBin(filename);
  }
  // -------------------------------------------------------------------------
  unsigned int SpecSet::LoadSpecSet_pklbin(const char *filename,
                                           bool countSpectraOnly)
  {
    if (countSpectraOnly == true)
    {
      WARN_MSG("countSpectraOnly is no longer supported");
    }
    return loadPklBin(filename);
  }
  // -------------------------------------------------------------------------
  void SpecSet::SaveScanNums(const char *filename)
  {
    vector<unsigned int> scanNums(specs.size());
    for (unsigned int i = 0; i < specs.size(); i++)
      scanNums[i] = specs[i].scan;
    Save_binArray(filename, scanNums);
  }
  // -------------------------------------------------------------------------
  short SaveSpecSet_pklbin(const char *filename, vector<Spectrum *> specs)
  {
    FILE *fp;
    unsigned int numSpecs = specs.size();
    unsigned int i, p;

    fp = fopen(filename, "w");
    if (fp == 0)
    {
      cerr << "ERROR: cannot open " << filename << "\n";
      return -1;
    }

    fwrite(&numSpecs, sizeof(int), 1, fp); // Number of spectra in the file

    unsigned short *numPeaks = (unsigned short *)malloc(sizeof(short)
        * numSpecs);
    unsigned short maxNumPeaks = 0;
    for (i = 0; i < numSpecs; i++)
    {
      numPeaks[i] = specs[i]->size();
      maxNumPeaks = max(maxNumPeaks, numPeaks[i]);
    }
    fwrite(numPeaks, sizeof(short), numSpecs, fp); // Number of peaks per spectrum in the file

    float *peaksBuffer = (float *)malloc(2 * (maxNumPeaks + 1) * sizeof(float));
    unsigned int pbIdx;
    for (i = 0; i < numSpecs; i++)
    {
      peaksBuffer[0] = specs[i]->parentMass;
      peaksBuffer[1] = (float)specs[i]->parentCharge;
      for (pbIdx = 2, p = 0; p < numPeaks[i]; p++)
      {
        peaksBuffer[pbIdx++] = (*specs[i])[p][0];
        peaksBuffer[pbIdx++] = (*specs[i])[p][1];
      }
      fwrite(peaksBuffer, sizeof(float), 2 * (numPeaks[i] + 1), fp); // [parentMass charge] followed by [masses intensities]
    }

    free(peaksBuffer);
    free(numPeaks);
    fclose(fp);
    return 1;
  }
}
