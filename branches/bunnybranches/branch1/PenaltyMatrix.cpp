#include "PenaltyMatrix.h"

// Module Includes
#include "DelimitedTextReader.h"

// System Includes
#include <stdio.h>
#include <math.h>

#define DEBUG_PENALTY 0

const float UNKNOWN_MULT = 1.0;

using namespace specnets;
using namespace std;

namespace PenaltyMatrix_const {
  std::map<float, float> emptyMap;
  string stringN("N");
  string stringD("D");
  string stringE("E");
  string stringQ("Q");
  string stringK("K");
  string stringI("I");
  string stringL("L");
}
using namespace PenaltyMatrix_const;

// -------------------------------------------------------------------------
PenaltyMatrix::PenaltyMatrix(AAJumps & jumps, float resolution /* = 1.0 */) :
  m_resolution(resolution), m_allSpectraAveragePeakIntensity(1.0), m_unknownPenalty(-1.0)
{
  if (jumps.aaLetters.size() == 0) {
    WARN_MSG("Initializing with empty amino acids");
  }
  // Create AA sequences (first index) in the penalty matrix
  for(int i = 0; i < jumps.aaLetters.size(); i++) {
    char aa = jumps.aaLetters[i];
    if (aa == 0) {
      continue;
    }
    if (DEBUG_PENALTY) DEBUG_VAR(aa);
    string strAA("X");  // Create a dummy string
    strAA[0] = aa;      // Set the string to the AA char
    std::map<float, float> newMap;
    if (DEBUG_PENALTY) DEBUG_VAR(strAA);
    penalties[strAA] = newMap; // Create a new map for this AA

    float mass = jumps.masses[i];

    if (DEBUG_PENALTY) DEBUG_VAR(mass);
    // Sanity check.. shouldn't happen (but seems to sometimes with grid)
    if (mass == 0.0) {
      continue;
    }
    mapCharMods[strAA] = mass; // map of AA to Mass
    mapModChars[mass] = strAA; // map of Mass to AA
  }

  // Fill in masses for I/L and K/Q if not already in map
  if (mapCharMods.find(stringL) == mapCharMods.end() && mapCharMods.find(stringI) != mapCharMods.end()) {
    if (DEBUG_PENALTY) DEBUG_VAR(stringL);
    mapCharMods[stringL] = mapCharMods[stringI];
    if (DEBUG_PENALTY) DEBUG_VAR(stringL);
    if (DEBUG_PENALTY) DEBUG_VAR(mapCharMods[stringL]);
  }
  if (mapCharMods.find(stringI) == mapCharMods.end() && mapCharMods.find(stringL) != mapCharMods.end()) {
    mapCharMods[stringI] = mapCharMods[stringL];
    if (DEBUG_PENALTY) DEBUG_VAR(stringI);
    if (DEBUG_PENALTY) DEBUG_VAR(mapCharMods[stringI]);
  }

  if (mapCharMods.find(stringK) == mapCharMods.end() && mapCharMods.find(stringQ) != mapCharMods.end()) {
    mapCharMods[stringK] = mapCharMods[stringQ];
    if (DEBUG_PENALTY) DEBUG_VAR(stringK);
    if (DEBUG_PENALTY) DEBUG_VAR(mapCharMods[stringK]);
  }
  if (mapCharMods.find(stringQ) == mapCharMods.end() && mapCharMods.find(stringK) != mapCharMods.end()) {
    mapCharMods[stringQ] = mapCharMods[stringK];
    if (DEBUG_PENALTY) DEBUG_VAR(stringQ);
    if (DEBUG_PENALTY) DEBUG_VAR(mapCharMods[stringQ]);
  }

  return;
}

// -------------------------------------------------------------------------
PenaltyMatrix::~PenaltyMatrix(void)
{
}

// -------------------------------------------------------------------------
float PenaltyMatrix::operator()(string & seq, float mass, float averageSpectrumPeakIntensity)
{
  if (penalties.find(seq) == penalties.end()) {
    if (averageSpectrumPeakIntensity != 0.0) {
      return m_unknownPenalty * averageSpectrumPeakIntensity;
    }
    return m_unknownPenalty * m_allSpectraAveragePeakIntensity;
  }

  float massRounded = roundMass(mass);

  if (penalties[seq].find(massRounded) == penalties[seq].end()) {
    return 0.0;
  }

  float basePenalty = penalties[seq][massRounded];
  // Known mods are not adjusted for average peak intensity
  if (m_knownMods[seq].find(massRounded) != m_knownMods[seq].end()) {
    return basePenalty;
  }
  if (averageSpectrumPeakIntensity != 0.0) {
    return basePenalty * averageSpectrumPeakIntensity;
  } else {
    return basePenalty * m_allSpectraAveragePeakIntensity;
  }

  return penalties[seq][massRounded] * averageSpectrumPeakIntensity;
}

// -------------------------------------------------------------------------
bool PenaltyMatrix::isKnown(string & seq, float mass)
{
  if (penalties.find(seq) == penalties.end()) {
    return false;
  }
  float massRounded = roundMass(mass);
  if (penalties[seq].find(massRounded) == penalties[seq].end()) {
    return false;
  }
  return true;
}

// -------------------------------------------------------------------------
void PenaltyMatrix::getPenalties(string & seq, std::map<float, float> & penaltyMap, float averageSpectrumPeakIntensity)
{
  map<string, map<float, float> >::iterator itrMap = penalties.find(seq);

  if (itrMap == penalties.end()) {
    return;
  }

  map<float, float>::iterator itr = itrMap->second.begin();
  map<float, float>::iterator itr_end = itrMap->second.end();
  for (; itr != itr_end; itr++) {

    float basePenalty = itr->second;
    // Known mods are not adjusted for average peak intensity
    if (m_knownMods[itrMap->first].find(itr->first) != m_knownMods[itrMap->first].end()) {
      penaltyMap[itr->first] = basePenalty;
    }
    if (averageSpectrumPeakIntensity != 0.0) {
      penaltyMap[itr->first] = basePenalty * averageSpectrumPeakIntensity;
    } else {
      penaltyMap[itr->first] = basePenalty * m_allSpectraAveragePeakIntensity;
    }
  }

  return;
}

// -------------------------------------------------------------------------
float PenaltyMatrix::getUnknownPenalty(float averageSpectrumPeakIntensity)
{
  if (averageSpectrumPeakIntensity != 0.0) {
    return m_unknownPenalty * averageSpectrumPeakIntensity;
  }
  return m_unknownPenalty * m_allSpectraAveragePeakIntensity;
}

//-----------------------------------------------------------------------------
// Load and process a BLOSUM matrix file
//-----------------------------------------------------------------------------
bool PenaltyMatrix::loadFromBlosum(std::string & filename)
{
  vector<string> aa;

  // Read the tab delimited file
  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(filename.c_str(),
                                                      " \t",
                                                      "",
                                                      lines)) {
    ERROR_MSG("Unable to open file! " << filename);
    return false;
  }

  float maxScore = 0.0;

  int firstLineIdx = -1;

  if (DEBUG_PENALTY) DEBUG_VAR(lines.size());
  for (int i = 0; i < lines.size(); i++) {
    // Discard empty lines
    if (DEBUG_PENALTY) DEBUG_VAR(lines[i].size());
    if (lines[i].size() == 0) {
      continue;
    }
    // Discard comments lines
    if (lines[i][0].compare("#") == 0) {
      continue;
    }

    if (firstLineIdx == -1) {
      firstLineIdx = i;
    }
    if (DEBUG_PENALTY) DEBUG_VAR(firstLineIdx);

    int aaIndex1 = i - firstLineIdx - 1;
    int aaIndex2 = 0;
    if (DEBUG_PENALTY) DEBUG_VAR(aaIndex1);
    if (DEBUG_PENALTY) DEBUG_VAR(aaIndex2);
    for (int j = 0; j < lines[i].size(); j++) {
      // Discard extra spaces and empty strings
      if (lines[i][j].length() == 0 || lines[i][j].compare(" ") == 0) {
        continue;
      }

      if (DEBUG_PENALTY) DEBUG_VAR(lines[i][j]);

      if (aaIndex1 == -1) {
        aa.push_back(lines[i][j]);
      } else {
        if (j == 0) {
          // Do nothing - the first char is the AA
        } else {
          string aa1 = aa[aaIndex1];
          int score;
          sscanf(lines[i][j].c_str(), "%d", &score);
          if (DEBUG_PENALTY) DEBUG_VAR(score);
          float mass1 = getMass(aa[aaIndex1]);
          float mass2 = getMass(aa[aaIndex2]);
          float massDiff = mass2 - mass1;
          // Round to the desired precision
          massDiff = roundMass(massDiff);

          if (massDiff != 0.0) {
            penalties[aa[aaIndex1]][massDiff]  = score;
            // Save the max score so we can subtract it out later so all scores are negative
            if (maxScore < score) {
              maxScore = score;
            }
          }

          aaIndex2++;
        }
      }

    } //for (int j = 0; j < lines[i].size(); j++) {

  } // for (int i = 0; i < lines.size(); i++) {

  if (DEBUG_PENALTY) DEBUG_VAR(maxScore);

  float minPenalty = -1.0;
  map<string, map<float, float> >::iterator itr = penalties.begin();
  map<string, map<float, float> >::iterator itr_end = penalties.end();
  for (; itr != itr_end; itr++) {
    map<float, float>::iterator itr2 = itr->second.begin();
    map<float, float>::iterator itr_end2 = itr->second.end();
    for (; itr2 != itr_end2; itr2++) {
      itr2->second -= maxScore;

      // Make the diagonals (mass diff 0) penalty be 0
      if (itr2->first == 0.0) {
        itr2->second = 0;
      }

      // LARS - DEBUG - EMPTY PENALTIES
      //itr2->second = 0;

      if (itr2->second < minPenalty) {
        minPenalty = itr2->second;
      }
    }
  }

  if (DEBUG_PENALTY) DEBUG_VAR(minPenalty);
  m_unknownPenalty = minPenalty * UNKNOWN_MULT;
  if (DEBUG_PENALTY) DEBUG_VAR(m_unknownPenalty);

#if DEBUG_PENALTY
  for (itr = penalties.begin(); itr != itr_end; itr++) {
    //DEBUG_VAR(itr->first);
    map<float, float>::iterator itr2 = itr->second.begin();
    map<float, float>::iterator itr_end2 = itr->second.end();
    for (; itr2 != itr_end2; itr2++) {
      DEBUG_MSG("[" << itr->first << "][" << itr2->first << "] = " << itr2->second);
    }
    DEBUG_TRACE;
  }
#endif

  return true;
}


//-----------------------------------------------------------------------------
// Create the penalties from modification count data
//-----------------------------------------------------------------------------
bool PenaltyMatrix::createFromModificationFreqs(map<float, float> & modFreq,
                                                float peakEquivalents,
                                                float minFrequency,
                                                float averagePeakIntensity)
{
  if (DEBUG_PENALTY) DEBUG_VAR(peakEquivalents);
  m_allSpectraAveragePeakIntensity = averagePeakIntensity;
  if (DEBUG_PENALTY) DEBUG_VAR(m_allSpectraAveragePeakIntensity);

  float maxPenalty = -1.0;

  // Transmute modification frequencies to penalties
  map<float, float>::iterator itr = modFreq.begin();
  map<float, float>::iterator itr_end = modFreq.end();
  for(; itr != itr_end; itr++) {
    float mass = itr->first;
    float freq = itr->second;

    // Some masses and frequencies may have been zeroed (undesireable)
    if (freq <= minFrequency) {
      continue;
    }

    // Is this an amino acid mass? If so ignore it.. its not a modification
    float roundedMass = roundMass(mass);
    if (mapModChars.find(roundedMass) != mapModChars.end()) {
      if (DEBUG_PENALTY) DEBUG_MSG("Mass is an amino acid [" << mass << "] at resolution [" << m_resolution << "]");
      continue;
    }
    if (DEBUG_PENALTY) DEBUG_MSG("Mass [" << mass << "] has frequency [" << freq << "]");

    map<string, map<float, float> > ::iterator itrP = penalties.begin();
    map<string, map<float, float> > ::iterator itrP_end = penalties.end();
    for(; itrP != itrP_end; itrP++) {

      float modf = modFreq[roundedMass];
      if (modf == 0.0) {
        continue;
      }

      // This amino acid, mass combination is known
      if (m_knownMods[itrP->first].find(roundedMass) != m_knownMods[itrP->first].end()) {
        penalties[itrP->first][roundedMass] = log10(modf);
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has known mod [" << roundedMass
                         << "], penalty = " << penalties[itrP->first][roundedMass]);
      } else {
        penalties[itrP->first][roundedMass] = peakEquivalents * log10(modf);
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has unknown mod [" << roundedMass
                         << "], penalty = " << penalties[itrP->first][roundedMass]);
      }

      float aaMass = mapCharMods[itrP->first];
      // This amino acid, negative mass combination is too small to be real
      if (aaMass - roundedMass < 57.0) {
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] with mod [" << -roundedMass
                                            << "], not allowed (total mass too small)");
        continue;
      // This amino acid, negative mass combination is known
      } else if (m_knownMods[itrP->first].find(-roundedMass) != m_knownMods[itrP->first].end()) {
        penalties[itrP->first][-roundedMass] = log10(modf);
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has known mod [" << -roundedMass
                         << "], penalty = " << penalties[itrP->first][-roundedMass]);
      } else {
        penalties[itrP->first][-roundedMass] = peakEquivalents * log10(modf);
        if (DEBUG_PENALTY) DEBUG_MSG("AA [" << itrP->first << "] has unknown mod [" << -roundedMass
                         << "], penalty = " << penalties[itrP->first][-roundedMass]);
      }

      if (penalties[itrP->first][roundedMass] < maxPenalty) {
        maxPenalty = penalties[itrP->first][roundedMass];
      }
      if (penalties[itrP->first][-roundedMass] < maxPenalty) {
        maxPenalty = penalties[itrP->first][-roundedMass];
      }

      // LARS - DEBUG - EMPTY PENALTIES
      //penalties[itrP->first][mass] = 0;

    } // penalties iterator

  } // modFreq iterator

  m_unknownPenalty = maxPenalty * UNKNOWN_MULT;
  if (DEBUG_PENALTY) DEBUG_VAR(m_unknownPenalty);

#if DEBUG_PENALTY
  for (map<string, map<float, float> >::iterator itr = penalties.begin(); itr != penalties.end(); itr++) {
    map<float, float>::iterator itr2 = itr->second.begin();
    map<float, float>::iterator itr_end2 = itr->second.end();
    for (; itr2 != itr_end2; itr2++) {
      DEBUG_MSG("[" << itr->first << "][" << itr2->first << "] = " << itr2->second);
    }
    DEBUG_TRACE;
  }
#endif

  return true;
}


//-----------------------------------------------------------------------------
float PenaltyMatrix::getMass(string aa)
{
  if (mapCharMods.size() == 0) {
    WARN_MSG("Internal modification map is empty in getMass() method");
    return 0.0;
  }
  if (mapCharMods.find(aa) == mapCharMods.end()) {
    WARN_MSG("No mass for [" << aa << "]");
  }
  if (mapCharMods[aa] == 0.0) {
    WARN_MSG("Mass of [" << aa << "] is 0.0");
  }

  return mapCharMods[aa];
}

//-----------------------------------------------------------------------------
float PenaltyMatrix::roundMass(float mass)
{
  int sign = mass < 0 ? -1 : 1;
  return (int)(fabs(mass) * 0.9995 + 0.5) * sign;
  //return abs((float)((int)(mass / m_resolution)) * m_resolution);
}


//-----------------------------------------------------------------------------
// Load the known modifications
//-----------------------------------------------------------------------------
bool PenaltyMatrix::loadKnownModifications(std::string & filename)
{
  if (DEBUG_PENALTY) DEBUG_TRACE;
  // Read the tab delimited file
  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(filename.c_str(),
                                                      " \t",
                                                      "",
                                                      lines)) {
    ERROR_MSG("Unable to open file! " << filename);
    return false;
  }

  string aa;
  if (DEBUG_PENALTY) DEBUG_VAR(lines.size());
  for (int i = 0; i < lines.size(); i++) {

    // Discard empty lines
    if (lines[i].size() == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }
    // Discard comments lines
    if (lines[i][0].compare("#") == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }

    float mass = 0.0;
    float roundedMass = 0.0;
    if (DEBUG_PENALTY) DEBUG_VAR(lines[i].size());
    for (int j = 0; j < lines[i].size(); j++) {

      if (j == 0) {
        // First item is the amino acid
        aa = lines[i][j];
      } else {
        // Rest of the items are masses of the modifications
        sscanf(lines[i][j].c_str(), "%f", &mass);
        if (DEBUG_PENALTY) DEBUG_VAR(mass);
        roundedMass = roundMass(mass);
        if (DEBUG_PENALTY) DEBUG_VAR(roundedMass);
        m_knownMods[aa].insert(roundedMass);
        if (DEBUG_PENALTY) DEBUG_MSG("Known modification of amino acid [" << aa << "] [" << roundedMass << "]");
      }

    } //for (int j = 0; j < lines[i].size(); j++) {

  } // for (int i = 0; i < lines.size(); i++) {

  return true;
}

//-----------------------------------------------------------------------------
// Load a generic matrix
//-----------------------------------------------------------------------------
bool PenaltyMatrix::load(std::string & filename)
{
  if (DEBUG_PENALTY) DEBUG_TRACE;
  // Read the tab delimited file
  vector<vector<string> > lines;
  if (!DelimitedTextReader::loadDelimitedFileNoHeader(filename.c_str(),
                                                      " \t",
                                                      "",
                                                      lines)) {
    ERROR_MSG("Unable to open file! " << filename);
    return false;
  }

  float minPenalty = 0.0;

  if (DEBUG_PENALTY) DEBUG_VAR(lines.size());
  for (int i = 0; i < lines.size(); i++) {

    // Discard empty lines
    if (lines[i].size() == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }
    // Discard comments lines
    if (lines[i][0].compare("#") == 0) {
      if (DEBUG_PENALTY) DEBUG_TRACE;
      continue;
    }

    string aa;
    float mass = 0.0;
    float roundedMass = 0.0;
    float penalty = 0.0;
    if (DEBUG_PENALTY) DEBUG_VAR(lines[i].size());
    for (int j = 0; j < lines[i].size(); j++) {

      if (j == 0) {
        // First item is the amino acid
        aa = lines[i][j];
      } else if (j == 1) {
        // Second item is the mass of the modification
        sscanf(lines[i][j].c_str(), "%f", &mass);
        if (DEBUG_PENALTY) DEBUG_VAR(mass);
        roundedMass = roundMass(mass);
        if (DEBUG_PENALTY) DEBUG_VAR(roundedMass);
      } else if (j == 2) {
        // Third item is the penalty
        sscanf(lines[i][j].c_str(), "%f", &penalty);
        if (DEBUG_PENALTY) DEBUG_VAR(penalty);

        penalties[aa][roundedMass] = penalty;

        if (penalty < minPenalty) {
          minPenalty = penalty;
        }
      }

    } //for (int j = 0; j < lines[i].size(); j++) {

  } // for (int i = 0; i < lines.size(); i++) {

  m_unknownPenalty = minPenalty * UNKNOWN_MULT;
  DEBUG_VAR(m_unknownPenalty);

  return true;
}

//-----------------------------------------------------------------------------
// Get the known modifications
//-----------------------------------------------------------------------------
const map<string, set<float> > & PenaltyMatrix::getKnownMods(void)
{
  return m_knownMods;
}

//-----------------------------------------------------------------------------
// Save the penalty matrix
//-----------------------------------------------------------------------------
bool PenaltyMatrix::saveMatrix(std::string & filename)
{
  ofstream ofs(filename.c_str(), ios_base::out | ios_base::binary);
  if (!ofs) {
    return false;
  }

  std::map<std::string, std::map<float, float> >::iterator itr = penalties.begin();
  std::map<std::string, std::map<float, float> >::iterator itrEnd = penalties.end();
  for ( ; itr != itrEnd; itr++) {
    std::map<float, float>::iterator itrMap = itr->second.begin();
    std::map<float, float>::iterator itrMapEnd = itr->second.end();
    for ( ; itrMap != itrMapEnd; itrMap++) {
      ofs << itr->first << "\t" << itrMap->first <<  "\t" << itrMap->second << endl;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Save the known mods
//-----------------------------------------------------------------------------
bool PenaltyMatrix::saveKnownMods(std::string & filename)
{
  ofstream ofs(filename.c_str(), ios_base::out | ios_base::binary);
  if (!ofs) {
    return false;
  }

  map<string, set<float> >::iterator itr = m_knownMods.begin();
  map<string, set<float> >::iterator itrEnd = m_knownMods.end();
  for ( ; itr != itrEnd; itr++) {
    ofs << itr ->first;
    set<float>::iterator itrSet = itr->second.begin();
    set<float>::iterator itrSetEnd = itr->second.end();
    for ( ; itrSet != itrSetEnd; itrSet++) {
      ofs << "\t" << *itrSet;
    }
    ofs << endl;
  }
  ofs.close();

  return true;
}

//-----------------------------------------------------------------------------
// Save the known mods
//-----------------------------------------------------------------------------
bool PenaltyMatrix::saveAminoAcids(std::string & filename)
{
  ofstream ofs(filename.c_str(), ios_base::out | ios_base::binary);
  if (!ofs) {
    return false;
  }

  ofs << mapCharMods.size() << endl;
  std::map<std::string, float>::iterator itr =  mapCharMods.begin();
  std::map<std::string, float>::iterator itrEnd = mapCharMods.end();
  for ( ; itr != itrEnd; itr++) {
    ofs << itr ->first << "=" << itr->second << endl;
  }
  ofs.close();

  return true;
}

