#include "PeptideSpectrumMatchSet.h"

namespace specnets
{
  // -------------------------------------------------------------------------
  psmPtr & PeptideSpectrumMatchSet::operator[](unsigned int i)
  {
    return m_psmSet[i];
  }

  // -------------------------------------------------------------------------
  const psmPtr & PeptideSpectrumMatchSet::operator[](unsigned int i) const
  {
    return m_psmSet[i];
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::removePsmSetItem(psmPtr item)
  {
    for (int i = 0; i < m_psmSet.size(); i++)
    {

      std::tr1::shared_ptr<PeptideSpectrumMatch> p =
          std::tr1::dynamic_pointer_cast<PeptideSpectrumMatch>(m_psmSet[i]); //I AM A PSMNETWORK POINTER

      // 	  		PeptideSpectrumMatch * p = (PeptideSpectrumMatch *)m_psmSet[i].get();

      if (p == item)
      {
        m_psmSet.erase(m_psmSet.begin() + i, m_psmSet.begin() + i + 1);
        return;
      }
    }
    DEBUG_MSG("PeptideSpectrumMatchSet::removePsmSetItem:: pointer not found!");
  }
  // -------------------------------------------------------------------------
  PeptideSpectrumMatchSet & PeptideSpectrumMatchSet::operator=(PeptideSpectrumMatchSet &other)
  {
    m_psmSet.clear();
    for (int i = 0; i < other.m_psmSet.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);
      *currMatch = *(other.m_psmSet[i]);
      m_psmSet.push_back(currMatch);
    }
  }
  // -------------------------------------------------------------------------
  PeptideSpectrumMatchSet & PeptideSpectrumMatchSet::push_back(psmPtr &other)
  {
    psmPtr currMatch(new PeptideSpectrumMatch);
    *currMatch = *other;
    m_psmSet.push_back(currMatch);
  }
  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::size()
  {
    return (unsigned int)m_psmSet.size();
  }
  // -------------------------------------------------------------------------
  unsigned int PeptideSpectrumMatchSet::resize(unsigned int newSize)
  {
    m_psmSet.resize(newSize);
    return m_psmSet.size();
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::parseLine(psmPtr currMatch,
                                          vector<string> &line,
                                          int zeroIndexed,
                                          int scanIndex,
                                          int specIdxIndex,
                                          int spectrumFileIndex,
                                          int annotationIndex,
                                          int chargeIndex,
                                          int scoreIndex,
                                          int pvalueIndex,
                                          int proteinIndex,
                                          bool isInspect)
  {

    bool haveScanNum = false;
    if (scanIndex != -1)
    {
      if (!sscanf(line[scanIndex].c_str(), "%d", &(currMatch->m_scanNum)))
      {
        ERROR_MSG("Unable to get scan number!");
        return false;
      }
      if (currMatch->m_scanNum > 0)
      {
        haveScanNum = true;
      }
      if (haveScanNum && zeroIndexed)
      {
        currMatch->m_scanNum += 1;
      }
    }
    if (!haveScanNum && specIdxIndex != -1)
    {
      if (!sscanf(line[specIdxIndex].c_str(), "%d", &(currMatch->m_scanNum)))
      {
        ERROR_MSG("Unable to get scan number from spectrum index!");
        return false;
      }
      if (zeroIndexed)
      {
        currMatch->m_scanNum += 1;
      }
    }

    if (spectrumFileIndex != -1)
    {
      currMatch->m_spectrumFile = line[spectrumFileIndex];
    }

    if (annotationIndex != -1)
    {
      if (isInspect)
      {
        currMatch->inspectToSpecNets(line[annotationIndex],
                                     currMatch->m_annotation);
      }
      else
      {
        currMatch->m_annotation = line[annotationIndex];
      }
    }

    if (chargeIndex != -1)
    {
      if (!sscanf(line[chargeIndex].c_str(), "%d", &(currMatch->m_charge)))
      {
        ERROR_MSG("Unable to get charge!");
        return false;
      }
    }

    if (scoreIndex != -1)
    {
      if (!sscanf(line[scoreIndex].c_str(), "%f", &(currMatch->m_score)))
      {
        ERROR_MSG("Unable to get score!");
        return false;
      }
    }

    if (pvalueIndex != -1)
    {
      if (!sscanf(line[pvalueIndex].c_str(), "%lf", &(currMatch->m_pValue)))
      {
        ERROR_MSG("Unable to get p-value");
        return false;
      }
    }

    if (proteinIndex != -1)
    {
      currMatch->m_protein = line[proteinIndex];
    }

    if (currMatch->m_protein.substr(0, 3).compare("XXX") == 0)
    {
      currMatch->m_isDecoy = true;
    }
    else
    {
      currMatch->m_isDecoy = false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::getPSMSet(SpecSet * spectra)
  {
    m_psmSet.clear(); //empty any current psm annotations

    for (int i = 0; i < spectra->size(); i++)
    {
      Spectrum spectrum = (*spectra)[i];
      if (spectrum.psmList.size() > 0)
      {
        list<psmPtr>::iterator it;

        for (it = spectrum.psmList.begin(); it != spectrum.psmList.end(); it++)
        {
          psmPtr currMatch(new PeptideSpectrumMatch);
          (*currMatch) = **it;
          currMatch->m_scanNum = spectrum.scan;
          m_psmSet.push_back(currMatch);
        }
      }
    }
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::addSpectra(SpecSet * spectra, bool addToSpectra)
  {
    vector<psmPtr>::iterator it;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      int scanNum = (*it)->m_scanNum;
      Spectrum * spectrum = spectra->getScan(scanNum); //getScan returns null if spectrum is not found

      if (spectrum == NULL)
      {
        WARN_MSG("Unable to find spectrum for scan " << scanNum);
      }
      else
      {
        if (addToSpectra)
        {
          spectrum->psmList.push_back(*it);
        }
      }
      (*it)->m_spectrum = spectrum;
    }
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatchSet::addSpectra(SpecSet * spectra,
                                           string filename,
                                           bool addToSpectra)
  {
    vector<psmPtr>::iterator it;
    for (it = m_psmSet.begin(); it != m_psmSet.end(); it++)
    {
      string psmFile;
      extractFileName((*it)->m_spectrumFile, psmFile, false);

      if (filename.compare(psmFile) == 0)
      {
        int scanNum = (*it)->m_scanNum;

        Spectrum * spectrum = spectra->getScan(scanNum);
        if (spectrum == NULL)
        {
          WARN_MSG("Unable to find spectrum for scan " << scanNum);
        }
        else
        {
          if (addToSpectra)
          {
            spectrum->psmList.push_back(*it);
          }
        }
        (*it)->m_spectrum = spectrum;
      }
    }
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadFromFile(const char * filename)
  {
    // Read the tab delimited file
    vector < vector<string> > lines;
    if (!DelimitedTextReader::loadDelimitedFileNoHeader(filename,
                                                        "\t",
                                                        "",
                                                        lines))
    {
      ERROR_MSG("Unable to open file! " << filename);
      return false;
    }

    for (int i = 0; i < lines.size(); i++)
    {
      if (lines[i].size() == 1)
      {
        break;
      }
      if (lines[i].size() != 15)
      {
        ERROR_MSG("Improper line length [" << lines[i].size() << "]");
        return false;
      }
      psmPtr currMatch(new PeptideSpectrumMatch);

      // Parse results into m_psmSet
      if (!sscanf(lines[i][0].c_str(), "%d", &(currMatch->m_scanNum)))
      {
        ERROR_MSG("Unable to get scan number!");
        return false;
      }

      currMatch->m_spectrumFile = lines[i][1];
      currMatch->m_annotation = lines[i][2];
      currMatch->m_origAnnotation = lines[i][3];
      currMatch->m_protein = lines[i][4];

      if (!sscanf(lines[i][5].c_str(), "%d", &(currMatch->m_dbIndex)))
      {
        ERROR_MSG("Unable to get DB index!");
        return false;
      }
      if (!sscanf(lines[i][6].c_str(), "%d", &(currMatch->m_numMods)))
      {
        ERROR_MSG("Unable to get num mods!");
        return false;
      }
      if (!sscanf(lines[i][7].c_str(), "%d", &(currMatch->m_matchOrientation)))
      {
        ERROR_MSG("Unable to get match orientation!");
        return false;
      }
      if (!sscanf(lines[i][8].c_str(), "%f", &(currMatch->m_startMass)))
      {
        ERROR_MSG("Unable to get start mass!");
        return false;
      }
      if (!sscanf(lines[i][9].c_str(), "%d", &(currMatch->m_charge)))
      {
        ERROR_MSG("Unable to get charge!");
        return false;
      }
      if (!sscanf(lines[i][10].c_str(), "%f", &(currMatch->m_score)))
      {
        ERROR_MSG("Unable to get score!");
        return false;
      }
      if (!sscanf(lines[i][11].c_str(), "%lf", &(currMatch->m_pValue)))
      {
        ERROR_MSG("Unable to get p-value!");
        return false;
      }
      int isDecoy;
      if (!sscanf(lines[i][12].c_str(), "%d", &isDecoy))
      {
        ERROR_MSG("Unable to get score!");
        return false;
      }
      currMatch->m_isDecoy = isDecoy;
      if (!sscanf(lines[i][13].c_str(),
                  "%f",
                  &(currMatch->m_strict_envelope_score)))
      {
        ERROR_MSG("Unable to get strict envelope score!");
        return false;
      }
      if (!sscanf(lines[i][14].c_str(),
                  "%f",
                  &(currMatch->m_unstrict_envelope_score)))
      {
        ERROR_MSG("Unable to get non-strict envelope score!");
        return false;
      }

      m_psmSet.push_back(currMatch);
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::saveToFile(const char * filename,
                                           bool includeHeader)
  {
    ofstream ofs(filename);
    if (!ofs || !ofs.good())
    {
      ERROR_MSG("Unable to read file! " << filename);
      return false;
    }

    if (includeHeader)
    {
      ofs
          << "Scan#\t#SpectrumFile\tAnnotation\tOrigAnnotation\tProtein\tdbIndex\tnumMods\tmatchOrientation\tstartMass\tCharge\tMQScore\tp-value\tisDecoy\tStrictEnvelopeScore\tUnstrictEvelopeScore";
    }

    for (int i = 0; i < size(); i++)
    {
      m_psmSet[i]->saveToFile(ofs);
    }
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadSpecnetsResultsFile(const char * resultsFile,
                                                        bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up inspect headers
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      else if (header[i].compare("#SpectrumFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Annotation") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("MQScore") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("p-value") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
    }

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     -1,
                     spectrumFileIndex,
                     annotationIndex,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     proteinIndex,
                     false))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << endl; */
        m_psmSet.push_back(currMatch);
      }
    }

    return true;
  }
  // -------------------------------------------------------------------------

  bool PeptideSpectrumMatchSet::loadInspectResultsFiles(const char * resultsFileList)
  {
    vector < vector<string> > inspectFileList;
    map<string, unsigned int> inspectFileListHeader;
    vector < string > requiredHeader;
    requiredHeader.push_back("Path");
    vector<int> requiredHeaderIndex;

    if (!DelimitedTextReader::loadDelimitedFile(resultsFileList,
                                                "\t",
                                                "#",
                                                inspectFileListHeader,
                                                inspectFileList,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to load inspect files!" << resultsFileList);
      return false;
    }

    //set up column numbers for optional parameters.
    int isZeroIndexedColumn = -1;

    map<string, unsigned int>::iterator it;

    it = inspectFileListHeader.find("isZeroIndexed");

    if (it != inspectFileListHeader.end())
    {
      isZeroIndexedColumn = it->second;
    }

    for (int i = 0; i < inspectFileList.size(); i++)
    {
      //get path
      const char * path = inspectFileList[i][requiredHeaderIndex[0]].c_str();
      DEBUG_VAR(path);

      //make sure results are loaded with correct indices.
      if (isZeroIndexedColumn != -1)
      {
        //make sure that sscanf doesn't overrun. Set to int then cast to bool.
        int isZeroIndexed = 0;
        if (!sscanf(inspectFileList[i][isZeroIndexedColumn].c_str(),
                    "%d",
                    &isZeroIndexed))
        {
          WARN_MSG("Unable to convert isZeroIndexed to bool!");
        }

        isZeroIndexed = (bool)isZeroIndexed;

        if (!this->loadInspectResultsFile(path, isZeroIndexed))
        {
          ERROR_MSG("Unable to load inspect file! " << inspectFileList[i][requiredHeaderIndex[0]]);
          return false;
        }
      }
      else
      {
        if (!this->loadInspectResultsFile(path))
        {
          ERROR_MSG("Unable to load inspect file! " << inspectFileList[i][requiredHeaderIndex[0]]);
          return false;
        }
      }
    }
    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadInspectResultsFile(const char * resultsFile,
                                                       bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up inspect headers
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Annotation");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      else if (header[i].compare("#SpectrumFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Annotation") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("MQScore") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("p-value") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
    }

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     -1,
                     spectrumFileIndex,
                     annotationIndex,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     proteinIndex,
                     true))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << "\t" << currMatch->m_isDecoy << endl;*/
        m_psmSet.push_back(currMatch);
      }
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSet::loadMSGFDBResultsFile(const char * resultsFile,
                                                      bool mergeConsecutive,
                                                      bool updateScanIdxs,
                                                      bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up inspect headers
    requiredHeader.push_back("SpecIndex");
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Peptide");
    requiredHeader.push_back("Charge");
    requiredHeader.push_back("Protein");

    //need these in order to identify PSMs from the same precursor
    if (mergeConsecutive)
    {
      requiredHeader.push_back("Precursor");
      requiredHeader.push_back("FragMethod");
    }

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int specIdxIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;
    int precusorIndex = -1;
    int fragIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("SpecIndex") == 0)
      {
        specIdxIndex = i;
      }
      else if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      else if (header[i].compare("#SpecFile") == 0)
      {
        spectrumFileIndex = i;
      }
      else if (header[i].compare("Peptide") == 0)
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("SpecProb") == 0)
      {
        scoreIndex = i;
      }
      else if (header[i].compare("P-value") == 0)
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
      else if (header[i].compare("Precursor") == 0)
      {
        precusorIndex = i;
      }
      else if (header[i].compare("FragMethod") == 0)
      {
        fragIndex = i;
      }
    }

    //parse results into m_psmSet
    int numConsecutive = 1;
    bool haveNumConsec = false;
    string prevPrecursorStr("");
    string prevFragStr("");

    for (int i = 0; i < lines.size(); i++)
    {
      if (mergeConsecutive)
      {
        string precursorStr(lines[i][precusorIndex]);
        string fragStr(lines[i][fragIndex]);
        if (prevPrecursorStr == precursorStr && prevFragStr != fragStr)
        {
          if (!haveNumConsec)
          {
            numConsecutive++;
          }
          continue;
        }
        else if (prevPrecursorStr != precursorStr && prevFragStr == fragStr)
        {
          haveNumConsec = true;
        }
        prevFragStr = fragStr;
        prevPrecursorStr = precursorStr;
      }

      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     specIdxIndex,
                     spectrumFileIndex,
                     annotationIndex,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     proteinIndex,
                     true))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /* cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << endl; */
        if (updateScanIdxs && haveNumConsec && numConsecutive > 1)
        {
          int idx = currMatch.get()->m_scanNum;
          idx = (zeroIndexed) ? idx / numConsecutive : ((idx - 1)
              / numConsecutive) + 1;
          currMatch.get()->m_scanNum = idx;
        }
        m_psmSet.push_back(currMatch);
      }
    }

    if (updateScanIdxs && haveNumConsec && numConsecutive > 1)
    {
      int idx = m_psmSet[0].get()->m_scanNum;
      idx = (zeroIndexed) ? idx / numConsecutive : ((idx - 1) / numConsecutive)
          + 1;
      m_psmSet[0].get()->m_scanNum = idx;
    }

    return true;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatchSetSpectralLibraryLoader::loadSpecnetsResultsFile(const char * resultsFile,
                                                                             bool zeroIndexed)
  {
    vector < string > header;
    vector < vector<string> > lines;
    vector < string > requiredHeader;
    vector<int> requiredHeaderIndex;

    //set up inspect headers
    requiredHeader.push_back("Scan#");
    requiredHeader.push_back("Charge");

    if (!DelimitedTextReader::loadDelimitedFile(resultsFile,
                                                "\t",
                                                "",
                                                header,
                                                lines,
                                                requiredHeader,
                                                requiredHeaderIndex))
    {
      ERROR_MSG("Unable to open results! " << resultsFile);
      return false;
    }

    int scanIndex = -1;
    int spectrumFileIndex = -1;
    int annotationIndex = -1;
    int chargeIndex = -1;
    int scoreIndex = -1;
    int pvalueIndex = -1;
    int proteinIndex = -1;
    int strictenvIndex = -1;
    int unstrictenvIndex = -1;

    //map headers
    for (int i = 0; i < header.size(); i++)
    {
      if (header[i].compare("Scan#") == 0)
      {
        scanIndex = i;
      }
      else if ((header[i].compare("#SpectrumFile") == 0)
          or (header[i].compare("#SpecFile") == 0))
      {
        spectrumFileIndex = i;
      }
      else if ((header[i].compare("Annotation") == 0)
          or (header[i].compare("Peptide") == 0))
      {
        annotationIndex = i;
      }
      else if (header[i].compare("Charge") == 0)
      {
        chargeIndex = i;
      }
      else if (header[i].compare("MSGFScore") == 0)
      {
        scoreIndex = i;
      }
      //else if ((header[i].compare("p-value") == 0) or (header[i].compare("P-value") == 0))
      else if ((header[i].compare("SpecProb") == 0))
      {
        pvalueIndex = i;
      }
      else if (header[i].compare("Protein") == 0)
      {
        proteinIndex = i;
      }
      else if (header[i].compare("ScrictEnvelope") == 0)
      {
        strictenvIndex = i;
      }
      else if (header[i].compare("UnstrictEnvelope") == 0)
      {
        unstrictenvIndex = i;
      }
    }

    //parse results into m_psmSet
    for (int i = 0; i < lines.size(); i++)
    {
      psmPtr currMatch(new PeptideSpectrumMatch);

      if (!parseLine(currMatch,
                     lines[i],
                     zeroIndexed,
                     scanIndex,
                     spectrumFileIndex,
                     annotationIndex,
                     chargeIndex,
                     scoreIndex,
                     pvalueIndex,
                     proteinIndex,
                     strictenvIndex,
                     unstrictenvIndex))
      {
        WARN_MSG("Unable to parse line" << i);
      }
      else
      {
        /*cout << currMatch->m_spectrumFile << "\t" << currMatch->m_scanNum << "\t"
         << currMatch->m_annotation << endl; */
        m_psmSet.push_back(currMatch);
      }
    }
  }

  bool PeptideSpectrumMatchSetSpectralLibraryLoader::parseLine(psmPtr currMatch,
                                                               vector<string> &line,
                                                               int zeroIndexed,
                                                               int scanIndex,
                                                               int spectrumFileIndex,
                                                               int annotationIndex,
                                                               int chargeIndex,
                                                               int scoreIndex,
                                                               int pvalueIndex,
                                                               int proteinIndex,
                                                               int strictenvIndex,
                                                               int unstrictenvIndex)
  {

    if (scanIndex != -1)
    {
      if (!sscanf(line[scanIndex].c_str(), "%d", &(currMatch->m_scanNum)))
      {
        ERROR_MSG("Unable to get scan number!");
        return false;
      }
      //cout<<"Scan#: "<<line[scanIndex]<<endl;
    }

    if (spectrumFileIndex != -1)
    {
      currMatch->m_spectrumFile = line[spectrumFileIndex];
      //cout<<"File: "<<line[spectrumFileIndex]<<endl;
    }

    if (annotationIndex != -1)
    {

      currMatch->m_annotation = line[annotationIndex];
      //cout<<"Annotation: "<<line[annotationIndex]<<endl;
    }

    if (chargeIndex != -1)
    {
      if (!sscanf(line[chargeIndex].c_str(), "%d", &(currMatch->m_charge)))
      {
        ERROR_MSG("Unable to get charge!");
        return false;
      }
    }

    if (scoreIndex != -1)
    {
      if (!sscanf(line[scoreIndex].c_str(), "%f", &(currMatch->m_score)))
      {
        ERROR_MSG("Unable to get score!");
        return false;
      }
    }

    if (pvalueIndex != -1)
    {
      float temp;
      if (!sscanf(line[pvalueIndex].c_str(), "%f", &(temp)))
      {
        ERROR_MSG("Unable to get p-value");
        return false;
      }
      currMatch->m_pValue = temp;
    }

    if (proteinIndex != -1)
    {
      currMatch->m_protein = line[proteinIndex];
    }

    if (currMatch->m_protein.substr(0, 3).compare("XXX") == 0)
    {
      currMatch->m_isDecoy = true;
    }
    else
    {
      currMatch->m_isDecoy = false;
    }

    if (strictenvIndex != -1)
    {
      if (!sscanf(line[strictenvIndex].c_str(),
                  "%f",
                  &(currMatch->m_strict_envelope_score)))
      {
        ERROR_MSG("Unable to get strict envelope score");
        return false;
      }
    }

    if (unstrictenvIndex != -1)
    {
      if (!sscanf(line[unstrictenvIndex].c_str(),
                  "%f",
                  &(currMatch->m_unstrict_envelope_score)))
      {
        ERROR_MSG("Unable to get unstrict envelope score");
        return false;
      }
    }

    return true;
  }
}
