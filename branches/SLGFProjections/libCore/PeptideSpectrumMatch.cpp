#include "Logger.h"
#include "PeptideSpectrumMatch.h"

#define DEBUG_GMPFA 0
#define DEBUG_GAFMP 0


namespace specnets
{
  // -------------------------------------------------------------------------
  PeptideSpectrumMatch::PeptideSpectrumMatch()
  {
    m_spectrumFile = "";
    m_scanNum = -1;
    m_origAnnotation = "";
    m_annotation = "";
    m_protein = "";
    m_charge = 0;
    m_score = -1;
    m_pValue = -1;
    m_isDecoy = false;
    m_spectrum = (Spectrum *)0x0;
    m_peakAnnotations.resize(0);
    m_ionTypes.resize(0);
    m_strict_envelope_score = 0.f;
    m_unstrict_envelope_score = 0.f;
    m_protein = "";
    m_dbIndex = -1;
    m_numMods = -1;
    m_matchOrientation = 0;
    m_startMass = -1.0;
    m_notes = "";
    m_ionmode = "";
    //m_submission_user = "";
    //m_submission_id = "";
    //m_submission_date = "";
    //m_organism = "";
    //m_molecule_mass = -1.0;
    //m_compound_name = "";
 }
  // -------------------------------------------------------------------------
  PeptideSpectrumMatch::~PeptideSpectrumMatch()
  {
    //EMPTY
  }
  // -------------------------------------------------------------------------
  PeptideSpectrumMatch::PeptideSpectrumMatch(const PeptideSpectrumMatch &other)
  {
    internalCopy(other);
  }

  // -------------------------------------------------------------------------
  PeptideSpectrumMatch &PeptideSpectrumMatch::operator=(const PeptideSpectrumMatch &other)
  {
    internalCopy(other);
    return (*this);
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::internalCopy(const PeptideSpectrumMatch &other)
  {
    m_spectrumFile = other.m_spectrumFile;
    m_scanNum = other.m_scanNum;
    m_annotation = other.m_annotation;
    m_origAnnotation = other.m_origAnnotation;
    m_peakAnnotations.resize(other.m_peakAnnotations.size());
    m_ionTypes.resize(other.m_ionTypes.size());

    map<const ftIonFragment*, const ftIonFragment*> pointerMapping; //mapping between old pointer (pointer_mapping.first) and new pointer (pointer_mapping.second

    for (unsigned int i = 0; i < other.m_peakAnnotations.size(); i++)
    {
      //set old annotation fragment pointer to new annotation fragment pointer.
      const ftIonFragment* oldFragPtr = other.m_peakAnnotations[i].first;
      m_peakAnnotations[i].first = pointerMapping[oldFragPtr];
      m_peakAnnotations[i].second = other.m_peakAnnotations[i].second;
    }

    m_matchedPeaks = other.m_matchedPeaks;

    for (unsigned int i = 0; i < other.m_ionTypes.size(); i++)
    {
      m_ionTypes[i] = other.m_ionTypes[i];
      pointerMapping[&(other.m_ionTypes[i])] = &(m_ionTypes[i]);
    }

    m_protein = other.m_protein;
    m_charge = other.m_charge;
    m_score = other.m_score;
    m_pValue = other.m_pValue;
    m_isDecoy = other.m_isDecoy;
    m_spectrum = other.m_spectrum;
    m_strict_envelope_score = other.m_strict_envelope_score;
    m_unstrict_envelope_score = other.m_unstrict_envelope_score;

    m_protein = other.m_protein;
    m_dbIndex = other.m_dbIndex;
    m_numMods = other.m_numMods;
    m_matchOrientation = other.m_matchOrientation;
    m_startMass = other.m_startMass;
    m_submission_metadata = other.m_submission_metadata;
    m_organism = other.m_organism;
    m_compound_name = other.m_compound_name;
    m_smiles = other.m_smiles;
    m_InChI = other.m_InChI;
    m_InChI_Aux = other.m_InChI_Aux;
    m_notes = other.m_notes;
    m_ionmode = other.m_ionmode;
    
    m_ion_extraction = other.m_ion_extraction;
    SLGF_distribution = other.SLGF_distribution;
    m_mz = other.m_mz;
    
    return;
  }

  // -------------------------------------------------------------------------
  /*
   * Helper function for annotate
   */
  bool compareProbs(const ftIonFragment &i, const ftIonFragment &j)
  {
    return (i.prob > j.prob);
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::setAnnotationToMatches(vector<int> &matches,
                                                    vector<pair<
                                                        const ftIonFragment*,
                                                        short> > &annotation,
                                                    int ionIdx,
                                                    const ftIonFragment* currIonFrag)
  {
    for (int i = 0; i < matches.size(); i++)
    {
      short annot_index = ionIdx + 1;
      const ftIonFragment* last_ion = annotation[matches[i]].first;
      if (last_ion && last_ion->prob >= currIonFrag->prob)
      {
        // if the last annotation had greater probabiliy, ignore current annotation
      }
      else
      {
        annotation[matches[i]].first = currIonFrag;
        annotation[matches[i]].second = annot_index;
      }
    }
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::setAnnotationToMatchesNoDups(vector<int> &matches,
                                                          vector<
                                                              pair<
                                                                  const ftIonFragment*,
                                                                  short> > &annotation,
                                                          int ionIdx,
                                                          const ftIonFragment* currIonFrag)
  {
    float maxIntensity = 0.0;
    int maxIntensityIndex = -1;

    for (int i = 0; i < matches.size(); i++)
    {
      float intensity = (*m_spectrum)[matches[i]][1];
      const ftIonFragment* last_ion = annotation[matches[i]].first;
      // Update max intensity
      if (last_ion && last_ion->prob >= currIonFrag->prob)
      {
        // if the last annotation had greater probabiliy, ignore current annotation
      }
      else
      {
        if (intensity > maxIntensity)
        {
          maxIntensity = intensity;
          maxIntensityIndex = matches[i];
        }
      }
    }

    short annot_index = ionIdx + 1;
    if (maxIntensityIndex > -1)
    {
      annotation[maxIntensityIndex].first = currIonFrag;
      annotation[maxIntensityIndex].second = annot_index;
    }

  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::setDbMatch(const string & protein,
                                        const int dbIndex,
                                        const float startMass)
  {
    m_protein = protein;
    m_dbIndex = dbIndex;
    m_startMass = startMass;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::getDbMatch(string & protein,
                                        int & dbIndex,
                                        float & startMass) const
  {
    protein = m_protein;
    dbIndex = m_dbIndex;
    startMass = m_startMass;
  }

  // -------------------------------------------------------------------------

  /**
   * helper function for annotate
   * sets m_ionTypes based on inclusion list
   * @param ionNamesInclude -  comma delimited string indicating which ions with which names should be annotated. If string is simply "all"
   * then just include all ions in MS2Model. ex. "y,b,y++,b++"
   * @param inputIonTypes - definition of ion types: mass offsets, ion probabilities, prefix/suffix.
   * is copied into m_ionTypes
   * @param ionTypes - output vector for fragments
   */
  bool copyIonTypes(const string &_ionNamesInclude,
                    const vector<ftIonFragment> &inputIonTypes,
                    vector<ftIonFragment> &outputIonTypes)
  {
    string ionNamesInclude(_ionNamesInclude);
    std::transform(ionNamesInclude.begin(),
                   ionNamesInclude.end(),
                   ionNamesInclude.begin(),
                   ::tolower);

    if (ionNamesInclude.compare("all") == 0)
    { //if we're including all ion types
      outputIonTypes.resize(inputIonTypes.size());
      copy(inputIonTypes.begin(), inputIonTypes.end(), outputIonTypes.begin());
      return inputIonTypes.size() == outputIonTypes.size(); // check to make sure our vector sizes match
    }
    else
    {
      vector < string > ionNames;

      splitText(ionNamesInclude.c_str(), ionNames, (const char*)",");

      for (int i = 0; i < inputIonTypes.size(); i++)
      {
        string lowerCaseName(inputIonTypes[i].name);
        std::transform(lowerCaseName.begin(),
                       lowerCaseName.end(),
                       lowerCaseName.begin(),
                       ::tolower);
        for (int j = 0; j < ionNames.size(); j++)
        {
          if (ionNames[j].compare(lowerCaseName) == 0)
          { //if current ion matches include list.
            ftIonFragment ion = inputIonTypes[i];
            outputIonTypes.push_back(ion); //copy ftIonFragment object.
            break;
          }
        }
      }
      return ionNames.size() == outputIonTypes.size();
    }
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::annotate(const string &peptide,
                                      const string &ionNamesInclude,
                                      const MS2ScoringModel &inputIonTypes,
                                      const float prmOffset,
                                      const float srmOffset,
                                      const float peakTol,
                                      const bool removeDuplicates,
                                      const bool ignoreParentCharge,
                                      const bool retainOldAnnotations)
  {
    AAJumps jumps(1);
    return annotate(peptide,
                    ionNamesInclude,
                    inputIonTypes,
                    prmOffset,
                    srmOffset,
                    peakTol,
                    jumps,
                    removeDuplicates,
                    ignoreParentCharge,
                    retainOldAnnotations);
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::annotate(const string &peptide,
                                      const string &ionNamesInclude,
                                      const MS2ScoringModel &inputIonTypes,
                                      const float prmOffset,
                                      const float srmOffset,
                                      const float peakTol,
                                      const AAJumps &jumps,
                                      const bool removeDuplicates,
                                      const bool ignoreParentCharge,
                                      const bool retainOldAnnotations)
  {
    //check to make sure spectrum is defined:
    if (m_spectrum == 0x0)
    {
      WARN_MSG("Spectrum not defined for PSM!");
      return false;
    }
    m_spectrum->rememberTolerances();
    m_spectrum->setPeakTolerance(peakTol);
    bool res = annotate(peptide,
                        ionNamesInclude,
                        inputIonTypes,
                        prmOffset,
                        srmOffset,
                        jumps,
                        removeDuplicates,
                        ignoreParentCharge,
                        retainOldAnnotations);
    m_spectrum->revertTolerances();

    return res;
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::annotate(const string &peptide,
                                      const string &ionNamesInclude,
                                      const MS2ScoringModel &inputIonTypes,
                                      const float prmOffset,
                                      const float srmOffset,
                                      const AAJumps &jumps,
                                      const bool removeDuplicates,
                                      const bool ignoreParentCharge,
                                      const bool retainOldAnnotations)
  {
    //check to make sure spectrum is defined:
    if (m_spectrum == 0x0)
    {
      WARN_MSG("Spectrum not defined for PSM!");
      return false;
    }
    m_annotation = peptide;

    if (!retainOldAnnotations || m_ionTypes.size() == 0)
    {
      //clear existing annotations
      m_peakAnnotations.clear();
      m_ionTypes.clear();

      //copy ftIonFragments into PSM ionTypes
      if (!copyIonTypes(ionNamesInclude, inputIonTypes.probs, m_ionTypes))
      {
        cerr << "Unable to copy ionNames from MS2ScoringModel! ionNamesInclude "
            << ionNamesInclude << endl;
        return false;
      }

      //resize m_peakAnnotations to contain annotations for each peak
      m_peakAnnotations.resize(m_spectrum->size());

      //initialize annotate vector to be null values
      for (int i = 0; i < m_peakAnnotations.size(); i++)
      {
        m_peakAnnotations[i].first = (ftIonFragment*)NULL;
        m_peakAnnotations[i].second = 0;
      }

      sort(m_ionTypes.begin(), m_ionTypes.end(), compareProbs); //sort fragment types in descending order so that low probability
      //annotations will be overwritten.
    }
    //generate srm and prm masses, no offsets.
    vector<float> prm_masses;
    vector<float> srm_masses;

    //generate total peptide mass (summed value of aa masses)
    float peptide_mass;

    jumps.getPRMandSRMMasses(peptide, prm_masses, srm_masses, peptide_mass);

    short peptide_length = prm_masses.size(); //length of peptide

    short ionIdx;
    int last_srm_index = -1; //keep track of last srm index so other srm searches start from there
    int last_prm_index = -1; //keep track of last prm index so other prm searches start from there

    for (ionIdx = 0; ionIdx < peptide_length - 1; ionIdx++)
    {

      vector<ftIonFragment>::const_iterator currIonFrag;

      for (currIonFrag = m_ionTypes.begin(); currIonFrag != m_ionTypes.end(); currIonFrag++)
      {
        vector<int> matches;
        if (currIonFrag->isIF) //we're looking at internal fragment
        {
          if (currIonFrag->isNTerm)
          {
            float curr_mass = (prm_masses[ionIdx] + currIonFrag->massOffset
                + prmOffset) / currIonFrag->charge;

#ifdef DEBUG
            std::cout << "mass " << curr_mass << endl;
            std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif

            if (m_charge > 0 && !ignoreParentCharge)
            {
              if (m_charge >= currIonFrag->charge)
                m_spectrum->findMatches(curr_mass,
                                        -1.0,
                                        matches,
                                        -1);
            }
            else
            {
              m_spectrum->findMatches(curr_mass,
                                          -1.0,
                                      matches,
                                      -1);
            }

            if (matches.size() > 0)
              last_prm_index = matches[matches.size() - 1];

          }
          else
          {
            float curr_mass = (srm_masses[ionIdx] + currIonFrag->massOffset
                + srmOffset) / currIonFrag->charge;
#ifdef DEBUG
            std::cout << "mass " << curr_mass << endl;
            std::cout << "name " << ionIdx << " " << currIonFrag->name << endl;
#endif

            if (m_charge > 0 && !ignoreParentCharge)
            {
              if (m_charge >= currIonFrag->charge)
                m_spectrum->findMatches(curr_mass,
                                                -1.0,
                                        matches,
                                        -1);
            }
            else
            {
              m_spectrum->findMatches(curr_mass,
                                      -1.0,
                                      matches,
                                      -1);
            }

            if (matches.size() > 0)
              last_srm_index = matches[matches.size() - 1];
          }
        }
        if (removeDuplicates)
          setAnnotationToMatchesNoDups(matches,
                                       m_peakAnnotations,
                                       ionIdx,
                                       &(*currIonFrag));
        else
          setAnnotationToMatches(matches,
                                 m_peakAnnotations,
                                 ionIdx,
                                 &(*currIonFrag));
      }
    }

    int last_pm_index = -1; //keep track of last pm index so other searches start from there.

    vector<ftIonFragment>::const_iterator currIonFrag;

    // looking at peptide mass shift (i.e., total peptide, total peptide - water, etc.)
    for (currIonFrag = m_ionTypes.begin(); currIonFrag != m_ionTypes.end(); currIonFrag++)
    {
      vector<int> matches;

      if (!currIonFrag->isIF)
      { // we're looking at peptide mass shift
        float curr_mass = (peptide_mass + currIonFrag->massOffset + prmOffset
            + srmOffset) / currIonFrag->charge;

        if (m_charge >= currIonFrag->charge)
          m_spectrum->findMatches(curr_mass, -1.0, matches, -1);

        if (matches.size() > 0)
          last_pm_index = matches[matches.size() - 1];
      }
      if (removeDuplicates)
        setAnnotationToMatchesNoDups(matches,
                                     m_peakAnnotations,
                                     ionIdx,
                                     &(*currIonFrag));
      else
        setAnnotationToMatches(matches,
                               m_peakAnnotations,
                               ionIdx,
                               &(*currIonFrag));
    }
    return true;
  }
  // -------------------------------------------------------------------------
  /*pair<int, float> PeptideSpectrumMatch::countAnnotatedPeaks(string& _ionNamesInclude,
   map<string, pair<int, float> >* outputIonCounts) {
   if (outputIonCounts != 0) {
   outputIonCounts->clear();
   }

   int numMatched = 0;
   float matchedIntensity = 0;

   string ionNamesInclude(_ionNamesInclude);
   std::transform(ionNamesInclude.begin(), ionNamesInclude.end(),
   ionNamesInclude.begin(), ::tolower);

   vector < string > ionNames;
   splitText(ionNamesInclude.c_str(), ionNames, (const char*)",");
   set < string > ionNamesSet;
   for (int i = 0; i < ionNames.size(); i++) {
   ionNamesSet.insert(ionNames[i]);
   }

   for (int i = 0; i < m_spectrum->size(); i++) {
   if (m_peakAnnotations[i].first == (ftIonFragment*)NULL) {
   continue;
   }
   string ionName(m_peakAnnotations[i].first->name);
   string ionNameLower(ionName);
   std::transform(ionNameLower.begin(), ionNameLower.end(),
   ionNameLower.begin(), ::tolower);

   if (ionNamesSet.count(ionNameLower) > 0) {
   numMatched++;
   matchedIntensity += (*m_spectrum)[i][1];

   if (outputIonCounts) {
   if (outputIonCounts->count(ionName) == 0) {
   pair<int, float> pairCount(1, (*m_spectrum)[i][1]);
   (*outputIonCounts)[ionName] = pairCount;
   } else {
   (*outputIonCounts)[ionName].first += 1;
   (*outputIonCounts)[ionName].second += (*m_spectrum)[i][1];
   }
   }
   }
   }

   return pair<int, float>(numMatched, matchedIntensity);
   }*/
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::isModified(void) const
  {
    size_t found = m_annotation.find_first_of("([");

    return found != string::npos;
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::setChargeByAnnotation(void)
  {
    AAJumps jumps(1);
    return setChargeByAnnotation(jumps);
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::setChargeByAnnotation(const AAJumps &jumps)
  {
    if (m_spectrum == 0x0)
    {
      WARN_MSG("Spectrum not defined for PSM!");
      return false;
    }

    if (m_spectrum->parentMZ == 0)
    {
      WARN_MSG("Parent mass to charge not defined for spectrum!");
      return false;
    }

    if (m_annotation.empty())
    {
      WARN_MSG("Annotation not defined for this peptide!");
      return false;
    }

    double charge = floor((jumps.getPeptideMass(m_annotation)
        / m_spectrum->parentMZ) + .5);

    stringstream ss;
    ss << charge;
    ss >> m_charge;
    return true;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::getMatchedPeaksFromAnnotation(Spectrum & dbSpec, 
                                                           AAJumps & aaJumps,
                                                           bool      useOriginal)
  {
    if (!m_spectrum || m_spectrum->size() == 0) {
      return;
    }

    vector<float> massValues;
    if (DEBUG_GMPFA) DEBUG_VAR(useOriginal);
    string annotation;
    if (useOriginal) {
      annotation = m_origAnnotation;
    } else {
      annotation = m_annotation;
    }
    if (annotation.empty()) {
      return;
    }
    if (DEBUG_GMPFA) DEBUG_VAR(annotation);
    vector<string> tokens;
    aaJumps.getPRMMasses(annotation, massValues, 0.0, &tokens, true);

    if (DEBUG_GMPFA) DEBUG_TRACE;
    vector<int> tokenLengths;
    for (int x=0; x < tokens.size(); x++) {
      if (DEBUG_GMPFA) DEBUG_VAR(tokens[x]);
      if (tokens[x][0] == '[') {
        tokenLengths.push_back(-1);
      } else {
        string stripped = aaJumps.stripMods(tokens[x]);
        if (DEBUG_GMPFA) DEBUG_VAR(stripped);
        tokenLengths.push_back(stripped.length());
      }
      if (DEBUG_GMPFA) DEBUG_VAR(tokenLengths[tokenLengths.size() - 1]);
    }
    bool startWithGap = tokenLengths[0] == -1;
    bool endWithGap = tokenLengths[tokenLengths.size() - 1] == -1;

    if (DEBUG_GMPFA) DEBUG_TRACE;
    // Find DB start index
    if (DEBUG_GMPFA) DEBUG_VAR(m_startMass);
    int dbStartIndex = 0;
    for (int idb = 1; idb < dbSpec.size(); idb++) {
      if (fabs(dbSpec[idb][0] - m_startMass) < 1.0) {
        dbStartIndex = idb;
        break;
      }
    }
    if (DEBUG_GMPFA) DEBUG_VAR(dbStartIndex);

    if (DEBUG_GMPFA) DEBUG_TRACE;
    if (DEBUG_GMPFA) DEBUG_VAR(massValues[0]);
    for (int i = 1; i < massValues.size(); i++) {
      if (DEBUG_GMPFA) DEBUG_VAR(massValues[i]);
      float massDiff = massValues[i] - massValues[i-1];
      if (DEBUG_GMPFA) DEBUG_VAR(massDiff);
    }
    if (DEBUG_GMPFA) DEBUG_VAR(massValues.size());


    if (DEBUG_GMPFA) DEBUG_TRACE;
    if (DEBUG_GMPFA) DEBUG_VAR((*m_spectrum)[0][0]);
    for (int i = 1; i < m_spectrum->size(); i++) {
      if (DEBUG_GMPFA) DEBUG_VAR((*m_spectrum)[i][0]);
      float massDiff = (*m_spectrum)[i][0] - (*m_spectrum)[i-1][0];
      if (DEBUG_GMPFA) DEBUG_VAR(massDiff);
    }
    if (DEBUG_GMPFA) DEBUG_VAR(m_spectrum->size());

    if (DEBUG_GMPFA) DEBUG_TRACE;
    m_matchedPeaks.clear();
    float peakTolerance = 1.0;
    vector<int> annoPeaks;
    vector<int> specPeaks;
    float firstSpecMass = startWithGap ? 0.0 : (*m_spectrum)[0][0];
    int is = startWithGap ? 1 : 0;
    int ia = startWithGap ? 1 : 0;
    while (ia < massValues.size() && is < m_spectrum->size()) {
      float diffMass;
      do {
        if (DEBUG_GMPFA) DEBUG_VAR(is);
        if (DEBUG_GMPFA) DEBUG_VAR(ia);
        if (DEBUG_GMPFA) DEBUG_VAR(massValues[ia]);
        float specMass = (*m_spectrum)[is][0] - firstSpecMass;
        if (DEBUG_GMPFA) DEBUG_VAR(specMass);
        diffMass = fabs(massValues[ia] - specMass);
        if (DEBUG_GMPFA) DEBUG_VAR(diffMass);
        if (diffMass < peakTolerance) {
          break;
        }
        if (specMass < massValues[ia]) {
          is++;
        } else {
          ia++;
        }

      } while (diffMass > peakTolerance && ia < massValues.size() && is < m_spectrum->size());

      if (ia < massValues.size() && is < m_spectrum->size()) {
        annoPeaks.push_back(ia);
        specPeaks.push_back(is);
        if (DEBUG_GMPFA) DEBUG_VAR(is);
        if (DEBUG_GMPFA) DEBUG_VAR(ia);
        is++;
        ia++;
      }
   } 

   if (DEBUG_GMPFA) DEBUG_TRACE;

   int id = dbStartIndex;
   int it = startWithGap ? 1 : 0;
   if (DEBUG_GMPFA) DEBUG_MSG("MATCHING PEAKS ARE:  " << specPeaks[0] << "  " << id);
   m_matchedPeaks.push_back(TwoValues<int>(specPeaks[0], id));

   int lengthMod = (endWithGap ? -1 : 0);
   for (int ip = 1; ip < annoPeaks.size() + lengthMod; ip++) {
     if (DEBUG_GMPFA) DEBUG_VAR(ip);
     if (DEBUG_GMPFA) DEBUG_VAR(annoPeaks[ip]);
     for (int x = 0; x < annoPeaks[ip] - annoPeaks[ip - 1]; x++) {
       id += tokenLengths[it];
       it++; 
     }
     if (DEBUG_GMPFA) DEBUG_MSG("MATCHING PEAKS ARE:  " << specPeaks[ip] << "  " << id);
     m_matchedPeaks.push_back(TwoValues<int>(specPeaks[ip], id));
   }

    return;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::getAnnotationFromMatchedPeaks(Spectrum & dbSpec, const string & dbSeqStr, string & annotation)
  {
    annotation.clear();  // Make sure string is clear since we will be appending to it
    // Annotate gap at beginning
    if (m_matchedPeaks.size()== 0)
    {
      WARN_MSG("No matched peaks!");
      annotation = m_annotation;
      return;
    }

    if (m_matchedPeaks[0][0] != 0 && m_spectrum) {
      char modString[32];
      sprintf(modString, "%.3f", (*m_spectrum)[m_matchedPeaks[0][0]][0]);
      annotation = "[";
      annotation += modString;
      annotation += "]";
    }

    int firstDbIndex = m_matchedPeaks[0][1];
    if (DEBUG_GAFMP) DEBUG_VAR(m_matchedPeaks.size());
    for (int k = 0; k < m_matchedPeaks.size() - 1; k++) {
      if (DEBUG_GAFMP) DEBUG_VAR(k);
      int thisDbIndex = m_matchedPeaks[k][1];
      int nextDbIndex = m_matchedPeaks[k+1][1];
      int peakWidth = nextDbIndex - thisDbIndex;
      int thisSpecIndex = m_matchedPeaks[k][0];
      int nextSpecIndex = m_matchedPeaks[k+1][0];
      if (DEBUG_GAFMP) DEBUG_VAR(thisDbIndex);
      if (DEBUG_GAFMP) DEBUG_VAR(nextDbIndex);
      if (DEBUG_GAFMP) DEBUG_VAR(peakWidth);
      if (DEBUG_GAFMP) DEBUG_VAR(thisSpecIndex);
      if (DEBUG_GAFMP) DEBUG_VAR(nextSpecIndex);

      float dbMassDiff = dbSpec[nextDbIndex][0] - dbSpec[thisDbIndex][0];
      float specMassDiff = (*m_spectrum)[nextSpecIndex][0] - (*m_spectrum)[thisSpecIndex][0];
      if (DEBUG_GAFMP) DEBUG_VAR(dbMassDiff);
      if (DEBUG_GAFMP) DEBUG_VAR(specMassDiff);

      float gapMassDiff = specMassDiff - dbMassDiff;

      if (fabs(gapMassDiff) < 0.5) {
        for (int a = 0; a < peakWidth; a++) {
          annotation += dbSeqStr[a + thisDbIndex];
        }
      } else {
        annotation += "(";
        for (int a = 0; a < peakWidth; a++) {
          annotation += dbSeqStr[a + thisDbIndex];
        }
        annotation += ",";
        char modString[32];
        sprintf(modString, "%.3f", gapMassDiff);
        annotation += modString;
        annotation += ")";
      }
      if (DEBUG_GAFMP) DEBUG_VAR(annotation);
    }
    if (DEBUG_GAFMP) DEBUG_VAR(annotation);

    // Annotate gap at end
    int lastMatchedPeak = m_matchedPeaks[m_matchedPeaks.size() - 1][0];
    if (m_spectrum && (lastMatchedPeak != m_spectrum->size() - 1)) {
      char modString[32];
      sprintf(modString, "%.3f", (*m_spectrum)[m_spectrum->size() - 1][0] - (*m_spectrum)[lastMatchedPeak][0]);
      annotation += "[";
      annotation += modString;
      annotation += "]";
    }
    if (DEBUG_GAFMP) DEBUG_VAR(annotation);

    return;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::getUnmodifiedPeptide(const string &inputPeptide,
                                                  string &outputPeptide)
  {
    outputPeptide = inputPeptide;
    //erase inspect syle mods
    size_t start = outputPeptide.find_first_of("(");

    while (start != string::npos)
    {
      outputPeptide.erase(start, 1);
      size_t modStart = outputPeptide.find_first_of(",", start);
      size_t modEnd = outputPeptide.find_first_of(")", modStart);
      outputPeptide.erase(modStart, modEnd - modStart + 1);

      start = outputPeptide.find_first_of("(", modStart);
    }
    //erase specnets style mods
    start = outputPeptide.find_first_of("[");

    while (start != string::npos)
    {
      outputPeptide.erase(start, 1);
      size_t modStart = start;
      size_t modEnd = outputPeptide.find_first_of("]", modStart);
      outputPeptide.erase(modStart, modEnd - modStart + 1);

      start = outputPeptide.find_first_of("[", modStart);
    }
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::generateTheoreticalMasses(const string &peptide,
                                                       const string &ionNamesInclude,
                                                       const MS2ScoringModel &inputIonTypes,
                                                       const float prmOffset,
                                                       const float srmOffset,
                                                       vector<string> &ionNames,
                                                       vector<float> &theoreticalMasses)
  {
    AAJumps jumps(1);
    return generateTheoreticalMasses(peptide,
                                     ionNamesInclude,
                                     inputIonTypes,
                                     prmOffset,
                                     srmOffset,
                                     jumps,
                                     ionNames,
                                     theoreticalMasses);

  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::generateTheoreticalMasses(const string &peptide,
                                                       const string &ionNamesInclude,
                                                       const MS2ScoringModel &inputIonTypes,
                                                       const float prmOffset,
                                                       const float srmOffset,
                                                       const AAJumps &jumps,
                                                       vector<string> &ionNames,
                                                       vector<float> &theoreticalMasses)
  {
    //generate srm and prm masses, no offsets.
    vector<float> prmMasses;
    vector<float> srmMasses;

    //generate total peptide mass (summed value of aa masses)
    float peptideMass;

    jumps.getPRMandSRMMasses(peptide, prmMasses, srmMasses, peptideMass);

    short peptideLength = prmMasses.size(); //length of peptide

    vector<ftIonFragment> ionTypes;
    //copy fragments from m_model.

    //copy ftIonFragments into PSM ionTypes
    if (!copyIonTypes(ionNamesInclude, inputIonTypes.probs, ionTypes))
    {
      cerr << "Unable to copy ionNames from MS2ScoringModel! ionNamesInclude "
          << ionNamesInclude << endl;
      return false;
    }
    sort(ionTypes.begin(), ionTypes.end(), compareProbs); //sort fragment types in ascending order so that low probability
    //annotations will be overwritten.

    short ionIdx;

    for (ionIdx = 0; ionIdx < peptideLength - 1; ionIdx++)
    {
      vector<ftIonFragment>::const_iterator currIonFrag;

      for (currIonFrag = ionTypes.begin(); currIonFrag != ionTypes.end(); currIonFrag++)
      {
        if (currIonFrag->isIF) //we're looking at internal fragment
        {
          float currMass;
          if (currIonFrag->isNTerm)
          {
            currMass
                = (prmMasses[ionIdx] + currIonFrag->massOffset + prmOffset)
                    / currIonFrag->charge;
          }
          else
          {
            currMass
                = (srmMasses[ionIdx] + currIonFrag->massOffset + srmOffset)
                    / currIonFrag->charge;
          }
          theoreticalMasses.push_back(currMass);
          stringstream ss;
          ss << currIonFrag->name << ionIdx + 1;
          string key = ss.str();
          ionNames.push_back(key);
        }
      }
    }

    vector<ftIonFragment>::const_iterator currIonFrag;

    // looking at peptide mass shift (i.e., total peptide, total peptide - water, etc.)
    for (currIonFrag = ionTypes.begin(); currIonFrag != ionTypes.end(); currIonFrag++)
    {
      if (!currIonFrag->isIF)
      { // we're looking at peptide mass shift
        float currMass = (peptideMass + currIonFrag->massOffset + prmOffset
            + srmOffset) / currIonFrag->charge;
        theoreticalMasses.push_back(currMass);
        stringstream ss;
        ss << currIonFrag->name << ionIdx + 1;
        string key = ss.str();
        ionNames.push_back(key);
      }
    }
    return true;
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::getUnmodifiedPeptide(string &outputPeptide) const
  {
    outputPeptide = m_annotation;
    size_t start = outputPeptide.find_first_of("(");

    while (start != string::npos)
    {
      outputPeptide.erase(start, 1);
      size_t modStart = outputPeptide.find_first_of(",", start);
      size_t modEnd = outputPeptide.find_first_of(")", modStart);
      outputPeptide.erase(modStart, modEnd - modStart + 1);

      start = outputPeptide.find_first_of("(", modStart);
    }
    //erase specnets style mods
    start = outputPeptide.find_first_of("[");

    while (start != string::npos)
    {
      outputPeptide.erase(start, 1);
      size_t modStart = start;
      size_t modEnd = outputPeptide.find_first_of("]", modStart);
      outputPeptide.erase(modStart, modEnd - modStart + 1);

      start = outputPeptide.find_first_of("[", modStart);
    }
  }
  // -------------------------------------------------------------------------
  int PeptideSpectrumMatch::countGaps(void)
  {
    if (m_spectrum == (Spectrum *)NULL)
    {
      WARN_MSG("Warning: Spectrum is not defined!");
      return -1.0;
    }

    Spectrum * spectrum = m_spectrum;

    AAJumps jumps(1);

    size_t start = m_annotation.find_first_of("[");

    int gapCount = 0;

    while (start != string::npos)
    {
      if (start == 1)
      {
        //not necessarily a gap, we need to check to make sure we're not looking at an n-term mod
        size_t end = m_annotation.find_first_of("]",start);
        float modification;
        stringstream ss;
        ss << m_annotation.substr(start + 1, end - start - 1);
        ss >> modification;
        if (modification > -100 && modification < 100)
        {
          gapCount++;
        }
      }
      else
      {
        gapCount++;
        start = m_annotation.find_first_of("[", start+1);
      }
    }
    return gapCount;
  }

  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::getModifications(vector<float> &modifications) const
  {
    size_t start = m_annotation.find_first_of("([");

    while (start != string::npos)
    {
      size_t mod;
      size_t end;

      if (m_annotation.substr(start, 1).compare("(") == 0)
      {
        mod = m_annotation.find_first_of(",", start);
        end = m_annotation.find_first_of(")", mod);
        float modification;
        stringstream ss;
        ss << m_annotation.substr(mod + 1, end - mod - 1);
        ss >> modification;
        if (ss.fail())
        {
          ss.clear();
          WARN_MSG("Unable to convert " << m_annotation.substr(mod+1, end-mod-1) << " to float!");
          return false;
        }
        else
        {
          modifications.push_back(modification);
        }
      }
      else
      {
        mod = start;
        end = m_annotation.find_first_of("]", mod);
        float modification;
        stringstream ss;
        ss << m_annotation.substr(mod + 1, end - mod - 1);
        ss >> modification;
        if (ss.fail())
        {
          ss.clear();
          WARN_MSG("Unable to convert " << m_annotation.substr(mod+1, end-mod-1) << " to float!");
          return false;
        }
        else
        {
          modifications.push_back(modification);
        }
      }
      start = m_annotation.find_first_of("([", end);
    }
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::insertModifications(const string &unmodifiedPeptide,
                                                 const vector<float> &modifications,
                                                 const vector<unsigned int> &positions,
                                                 string &outputPeptide)
  {
    stringstream ss;
    unsigned int peptidePosition = 0; //current peptide position.

    if (modifications.size() != positions.size())
    {
      ERROR_MSG("Unequal modifications and positions array sizes!");
      return false;
    }

    for (unsigned int i = 0; i < positions.size(); i++)
    {
      if (positions[i] >= unmodifiedPeptide.length())
      {
        ERROR_MSG("Modification position outside peptide!" << positions[i]);
        return false;
      }

      if (positions[i] < peptidePosition)
      {
        ERROR_MSG("Positions vector not sorted!");
        return false;
      }

      bool positionFound = false;
      while (!positionFound)
      {
        if (peptidePosition == positions[i])
        {
          ss << '(' << unmodifiedPeptide[peptidePosition] << ','
              << modifications[i] << ')';
          positionFound = true;
        }
        else
        {
          ss << unmodifiedPeptide[peptidePosition];
        }
        peptidePosition++;
      }
    }
    if (peptidePosition < unmodifiedPeptide.length())
    {
      ss << unmodifiedPeptide.substr(peptidePosition);
    }
    outputPeptide = ss.str();
    return true;
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::insertModifications(const vector<float> &modifications,
                                                 const vector<unsigned int> &positions,
                                                 string &outputPeptide) const
  {
    return insertModifications(m_annotation,
                               modifications,
                               positions,
                               outputPeptide);
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::mapIons(vector<vector<int> > &outputMatches,
                                     std::tr1::unordered_map<string, int> &ionMap) const
  {
    outputMatches.resize(0);

    for (int i = 0; i < m_peakAnnotations.size(); i++)
    {
      if (m_peakAnnotations[i].first != NULL)
      {
        stringstream ss;
        ss << m_peakAnnotations[i].first->name << m_peakAnnotations[i].second;
        string key = ss.str();

        if (ionMap.find(key) != ionMap.end())
        {
          outputMatches[ionMap[key]].push_back(i);
        }
        else
        {
          vector<int> newVector;
          newVector.push_back(i);
          outputMatches.push_back(newVector);
          ionMap[key] = outputMatches.size() - 1;
        }
      }
    }
  }

  // -------------------------------------------------------------------------

  void PeptideSpectrumMatch::inspectToSpecNets(const string &inspectOriginal,
                                               string &specnets)
  {

    //specnets = "*.";

    //check for trailing characters
        
        if (inspectOriginal.length() == 0) return;

    string inspect = inspectOriginal;
    size_t test = inspect.find_first_of(".");
    char firstChar = inspect.at(0);

    if (test == string::npos || (test != 1 && !(firstChar >= '0' && firstChar <= '9')))
    {
      stringstream ss;
      ss << "*." << inspect << ".*";
      inspect = ss.str();
    }

    int inspectLength = inspect.length();
    bool modification = false;
    bool startMass = false;

    const char * ptrInspect = inspect.c_str();
    char prevChar = 0x0;

    int i;
    //note: inspect_length - 3 means we ignore the trailing characters around annotation
    for (i = 2; i <= inspectLength - 2; i++)
    {
      char currChar = ptrInspect[i];

      //handle modifications
      if ('p' == currChar)
      {
        //phosphorylation
        specnets += "(";
        specnets += prevChar;
        specnets += ",80)";

        // look for starting + ou -
      }
      else if ((i == 2) && ('+' == currChar || '-' == currChar))
      {
        startMass = true;
        if ('-' == currChar)
          specnets += "[-";
        else
          specnets += "[";

        // look for a mass in the sequence (modification)
      }
      else if ('+' == currChar || '-' == currChar)
      {
        modification = true;

        if ('-' == currChar)
        {
          specnets += "(";
          specnets += prevChar;
          specnets += ",-";
        }
        else
        {
          specnets += "(";
          specnets += prevChar;
          specnets += ",";
        }

      }
      else if (modification)
      {
        if (0 == currChar) // check for null
        {
          // do nothing
        }
        else if ('9' >= currChar || '.' == currChar)
        { //this is a number (or decimal)
          specnets += currChar;
        }
        else
        {
          specnets += ')';
          modification = false;
        }
      }
      else if (startMass)
      {
        if (0 == currChar) // check for null
        {
          // do nothing
        }
        else if ('9' >= currChar || '.' == currChar)
        { //this is a number (or decimal)
          specnets += currChar;
        }
        else
        {
          specnets += ']';
          startMass = false;
        }
      }
      //handle normal AAs
      else if ('a' > currChar) //this is a capital letter
      {
        if (prevChar && 'a' > prevChar)
        {
          specnets += prevChar;
        }
      }
      prevChar = currChar;
    }

    if ('9' <= prevChar && '.' != prevChar)
    {
      specnets += prevChar;
    }

    if (modification)
    {
      specnets += ')';
    }

    if (startMass)
    {
      specnets += ']';
    }

    //specnets += ".*";
  }

  // -------------------------------------------------------------------------

  bool PeptideSpectrumMatch::loadFromFile(std::ifstream & ifs)
  {
    ifs >> m_spectrumFile;
    ifs >> m_scanNum;
    ifs >> m_annotation;
    ifs >> m_origAnnotation;
    ifs >> m_protein;
    ifs >> m_dbIndex;
    ifs >> m_numMods;
    ifs >> m_matchOrientation;
    ifs >> m_startMass;
    ifs >> m_charge;
    ifs >> m_score;
    ifs >> m_pValue;
    ifs >> m_isDecoy;
    ifs >> m_strict_envelope_score;
    ifs >> m_unstrict_envelope_score;
    
    
    
    std::string compound_name;
    ifs >> compound_name;
    m_compound_name.push_back(compound_name);

    return true;
  }

  // -------------------------------------------------------------------------

  bool PeptideSpectrumMatch::saveToFile(std::ofstream & ofs)
  {
    ofs << m_scanNum << "\t";
    if (m_spectrumFile.empty())
    {
      ofs << " " << "\t";
    }
    else
    {
      ofs << m_spectrumFile << "\t";
    }
    if (m_annotation.empty())
    {
      ofs << " " << "\t";
    }
    else
    {
      ofs << m_annotation << "\t";
    }
    if (m_origAnnotation.empty())
    {
      ofs << " " << "\t";
    }
    else
    {
      ofs << m_origAnnotation << "\t";
    }
    if (m_protein.empty())
    {
      ofs << " " << "\t";
    }
    else
    {
      ofs << m_protein << "\t";
    }
    ofs << m_dbIndex << "\t";
    ofs << m_numMods << "\t";
    ofs << m_matchOrientation << "\t";
    ofs << m_startMass << "\t";
    ofs << m_charge << "\t";
    ofs << m_score << "\t";
    ofs << m_pValue << "\t";
    ofs << m_isDecoy << "\t";
    ofs << m_strict_envelope_score << "\t";
    ofs << m_unstrict_envelope_score << "\t";
    
    if(m_compound_name.size() == 0){
        ofs << "N/A" ;
    }
    else{
        if(m_compound_name[m_compound_name.size()-1].length() == 0){
            ofs << "N/A" ;
        }
        else{
            ofs << m_compound_name[m_compound_name.size()-1];
        }
    }
    
    ofs << "\t";
    
    if(m_organism.size() == 0){
        ofs << "N/A" ;
    }
    else{
        if(m_organism[m_organism.size()-1].length() == 0){
            ofs << "N/A" ;
        }
        else{
            ofs << m_organism[m_organism.size()-1];
        }
    }
    
    ofs << "\t";
    
    if (m_spectrumFile.empty())
    {
      ofs << " ";
    }
    else
    {
        ofs << m_spectrumFile << m_scanNum;
    }
    
    ofs << endl;

    return true;
  }

}
