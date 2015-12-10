#include "FdrPeptide.h"

namespace specnets
{
  //---------------------------------------------------------------------------

  bool compareFDR(psmPtr i, psmPtr j)
  {
    return i->m_score > j->m_score;
  }

  // -------------------------------------------------------------------------
  bool FdrPeptide::calculatePValues(PeptideSpectrumMatchSet &inputPeptides,
                                    PeptideSpectrumMatchSet &outputPeptides)
  {
    outputPeptides = inputPeptides;

    sort(outputPeptides.m_psmSet.begin(),

    outputPeptides.m_psmSet.end(), compareFDR);

    unsigned int correctHits = 0;
    unsigned int incorrectHits = 0;

    for (unsigned int i = 0; i < outputPeptides.size(); i++)
    {
      psmPtr psm = outputPeptides[i];
      if (psm->m_isDecoy)
      {
        incorrectHits++;
      }
      else
      {
        correctHits++;
      }
      if (correctHits > 0)
      {
        outputPeptides[i]->m_pValue = (double)incorrectHits / (double)correctHits;
      }
      else
      {
        outputPeptides[i]->m_pValue = 1;
      }
    }
    if (incorrectHits == 0)
    {
      ERROR_MSG("No decoys set on inputPeptides!");
      return false;
    }
    else
    {
      return true;
    }
  }
  // helper function for calculatePValuesSeparateSearch,
  // adds psm with maximum m_score to map.
  // -------------------------------------------------------------------------
  bool addToMaxMap(PeptideSpectrumMatchSet &inputSet, std::tr1::unordered_map<
      Spectrum *, psmPtr> &maxPsmMap)
  {
    std::tr1::unordered_map<Spectrum *, psmPtr>::const_iterator it;

    for (unsigned int i = 0; i < inputSet.size(); i++)
    {
      psmPtr currMatch = inputSet[i];

      if (currMatch->m_spectrum == NULL)
      {
        WARN_MSG("Input PSM does not have associated spectrum.");
        return false;
      }

      it = maxPsmMap.find(currMatch->m_spectrum);

      if (it != maxPsmMap.end())
      {
        if (maxPsmMap[currMatch->m_spectrum]->m_score < currMatch->m_score)
        {
          maxPsmMap[currMatch->m_spectrum] = currMatch;
        }
      }
      else
      {
        maxPsmMap[currMatch->m_spectrum] = currMatch;
      }
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::mergeTargetDecoy(PeptideSpectrumMatchSet &targetPeptides,
                                    PeptideSpectrumMatchSet &decoyPeptides,
                                    PeptideSpectrumMatchSet &outputPeptides)
  {
    std::tr1::unordered_map<Spectrum *, psmPtr> maxPsmMap; //map between spectrum pointer and current max psm

    unsigned int correctHits = 0;
    unsigned int incorrectHits = 0;

    //cycle through target peptides first

    if (!addToMaxMap(targetPeptides, maxPsmMap))
    {
      ERROR_MSG("Unable to calculate p-values");
      return false;
    }

    if (!addToMaxMap(decoyPeptides, maxPsmMap))
    {
      ERROR_MSG("Unable to calculate p-values");
      return false;
    }

    std::tr1::unordered_map<Spectrum *, psmPtr>::const_iterator it;

    for (it = maxPsmMap.begin(); it != maxPsmMap.end(); it++)
    {
      psmPtr currPtr = it->second;
      outputPeptides.push_back(currPtr);
    }
    return true;
  }
  // -------------------------------------------------------------------------
  bool FdrPeptide::filterByPValue(PeptideSpectrumMatchSet &inputPeptides,
                                  double cutoff)
  {
    int count = 0;
    PeptideSpectrumMatchSet filteredSet;
    for (unsigned int i = 0; i < inputPeptides.size(); i++)
    {
      psmPtr currMatch = inputPeptides[i];

      if (currMatch->m_pValue == -1)
      {
        ERROR_MSG("PValue not set on input!");
        return false;
      }

      if (currMatch->m_pValue <= cutoff)
      {
        filteredSet.m_psmSet.push_back(currMatch);
      }
      else
      {
        break;
      }
    }
    inputPeptides = filteredSet;
    return true;
  }
}
