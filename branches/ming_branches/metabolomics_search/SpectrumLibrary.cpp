/*
 * SpectrumLibrary.cpp
 *
 *  Created on: Apr 21, 2011
 *      Author: jsnedecor
 */

#include "SpectrumLibrary.h"

namespace specnets
{

  /*
   * Generates key for an associated annotation;
   */
  string getKey(string &annotation, int charge)
  {
    stringstream key;
    key << annotation << '_' << charge;
    return key.str();
  }

  // -------------------------------------------------------------------------
  void SpectrumLibrary::buildKeyMap(void)
  {
    m_map.clear();

    for (int i = 0; i < m_library.size(); i++)
    {
      string key = getKey(m_library[i]->m_annotation, m_library[i]->m_charge);
      m_map[key] = i;
    }
  }

  // -------------------------------------------------------------------------
  psmPtr & SpectrumLibrary::operator[](unsigned int i)
  {
    return m_library[i];
  }

  // -------------------------------------------------------------------------
  const psmPtr & SpectrumLibrary::operator[](unsigned int i) const
  {
    return m_library[i];
  }

  // -------------------------------------------------------------------------
  SpectrumLibrary & SpectrumLibrary::operator=(SpectrumLibrary &other)
  {
    m_library.resize(other.m_library.size());
    m_library = other.m_library;
    buildKeyMap();
  }

  // -------------------------------------------------------------------------
  unsigned int SpectrumLibrary::size()
  {
    return (unsigned int)m_library.size();
  }

  // -------------------------------------------------------------------------
  unsigned int SpectrumLibrary::resize(unsigned int newSize)
  {
    m_library.resize(newSize);
    m_libraryScores.resize(newSize);
    buildKeyMap();
    return m_library.size();
  }

  // -------------------------------------------------------------------------
  // Helper function for addToLibrary
  void SpectrumLibrary::addSpectrumMatch(Spectrum &spectrum,
                                         PeptideSpectrumMatch &psm,
                                         string &key,
                                         float score)
  {
    psmPtr currPtr(new PeptideSpectrumMatch);
    *currPtr = psm;
    m_library.m_psmSet.push_back(currPtr);
    m_libraryScores.push_back(score);
    m_map[key] = m_library.size() - 1;
  }

  // -------------------------------------------------------------------------
  bool SpectrumLibrary::addToLibrary(Spectrum &spectrum,
                                     PeptideSpectrumMatch &psm,
                                     float score,
                                     bool scoreAscending)
  {
    string key = getKey(psm.m_annotation, psm.m_charge);
    std::tr1::unordered_map<string, unsigned int>::const_iterator it;
    it = m_map.find(key);

    if (it != m_map.end())
    {
      psmPtr libMatch = m_library[it->second];

      if (scoreAscending)
      {
        // if new score is better than old score
        if (m_libraryScores[it->second] > score)
        {
          addSpectrumMatch(spectrum, psm, key, score);
          return true;
        }
      }
      else if (m_libraryScores[it->second] < score)
      {
        // if new score is better than old score
        addSpectrumMatch(spectrum, psm, key, score);
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      addSpectrumMatch(spectrum, psm, key, score);
      return true;
    }
  }
  // -------------------------------------------------------------------------
  bool SpectrumLibrary::addToLibrary(Spectrum &spectrum,
                                     PeptideSpectrumMatch &psm,
                                     float score,
                                     string &key,
                                     bool scoreAscending)
  {
    std::tr1::unordered_map<string, unsigned int>::const_iterator it;
    it = m_map.find(key);

    if (it != m_map.end())
    {
      psmPtr libMatch = m_library[it->second];

      if (scoreAscending)
      {
        // if new score is better than old score
        if (m_libraryScores[it->second] > score)
        {
          addSpectrumMatch(spectrum, psm, key, score);
          return true;
        }
      }
      else if (m_libraryScores[it->second] < score)
      {
        // if new score is better than old score
        addSpectrumMatch(spectrum, psm, key, score);
        return true;
      }
      else
      {
        return false;
      }
    }
    else
    {
      addSpectrumMatch(spectrum, psm, key, score);
      return true;
    }
  }
  // -------------------------------------------------------------------------
  bool SpectrumLibrary::getSpectrumLibraryMatch(string &annotation,
                                                int charge,
                                                psmPtr &outputPsm)
  {
    string key = getKey(annotation, charge);
    DEBUG_MSG(key);
    std::tr1::unordered_map<string, unsigned int>::const_iterator it;
    it = m_map.find(key);

    if (it != m_map.end())
    {
      outputPsm = m_library[it->second];
      return true;
    }
    else
    {
      return false;
    }
  }
}

