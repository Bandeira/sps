#include "Logger.h"
#include "PeptideSpectrumMatch.h"

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
    m_charge = -1;
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
    
    m_submission_user = "";
    m_submission_id = "";
    m_submission_date = "";
    m_organism = "";
    m_molecule_mass = -1.0;
    m_compound_name = "";
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
    
    m_submission_user = other.m_submission_user;
    m_submission_id = other.m_submission_id;
    m_submission_date = other.m_submission_date;
    m_organism = other.m_organism;
    m_molecule_mass = other.m_molecule_mass;
    m_compound_name = other.m_compound_name;

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
  void PeptideSpectrumMatch::setDbMatch(string & protein,
                                        int dbIndex,
                                        float startMass)
  {
    m_protein = protein;
    m_dbIndex = dbIndex;
    m_startMass = startMass;
  }

  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::getDbMatch(string & protein,
                                        int & dbIndex,
                                        float & startMass)
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
  bool copyIonTypes(string &_ionNamesInclude,
                    vector<ftIonFragment> &inputIonTypes,
                    vector<ftIonFragment> &outputIonTypes)
  {
	string ionNamesInclude(_ionNamesInclude);
	std::transform(ionNamesInclude.begin(), ionNamesInclude.end(), ionNamesInclude.begin(),
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
    	std::transform(lowerCaseName.begin(), lowerCaseName.end(), lowerCaseName.begin(), ::tolower);
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
  bool PeptideSpectrumMatch::annotate(string &peptide,
                                      string &ionNamesInclude,
                                      MS2ScoringModel &inputIonTypes,
                                      float prmOffset,
                                      float srmOffset,
                                      float peakTol,
                                      bool removeDuplicates,
                                      bool ignoreParentCharge)
  {
    AAJumps jumps(1);
    return annotate(peptide,
                    ionNamesInclude,
                    inputIonTypes,
                    prmOffset,
                    srmOffset,
                    peakTol,
                    jumps,
                    removeDuplicates);
  }

  bool PeptideSpectrumMatch::annotate(string &peptide,
                                      string &ionNamesInclude,
                                      MS2ScoringModel &inputIonTypes,
                                      float prmOffset,
                                      float srmOffset,
                                      float peakTol,
                                      AAJumps &jumps,
                                      bool removeDuplicates,
                                      bool ignoreParentCharge) {
	//check to make sure spectrum is defined:
	if (m_spectrum == 0x0) {
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
			            ignoreParentCharge);
	m_spectrum->revertTolerances();

	return res;
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::annotate(string &peptide,
                                      string &ionNamesInclude,
                                      MS2ScoringModel &inputIonTypes,
                                      float prmOffset,
                                      float srmOffset,
                                      AAJumps &jumps,
                                      bool removeDuplicates,
                                      bool ignoreParentCharge)
  {
    //check to make sure spectrum is defined:
    if (m_spectrum == 0x0)
    {
      WARN_MSG("Spectrum not defined for PSM!");
      return false;
    }
    //clear existing annotations
    m_peakAnnotations.clear();
    m_ionTypes.clear();

    m_annotation = peptide;

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

    //generate srm and prm masses, no offsets.
    vector<float> prm_masses;
    vector<float> srm_masses;

    //generate total peptide mass (summed value of aa masses)
    float peptide_mass;

    jumps.getPRMandSRMMasses(peptide, prm_masses, srm_masses, peptide_mass);

    short peptide_length = prm_masses.size(); //length of peptide

    sort(m_ionTypes.begin(), m_ionTypes.end(), compareProbs); //sort fragment types in descending order so that low probability
    //annotations will be overwritten.

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
  bool PeptideSpectrumMatch::isModified(void)
  {
    size_t found = m_annotation.find_first_of("(");

    return found != string::npos;
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::getUnmodifiedPeptide(string &inputPeptide,
                                                  string &outputPeptide)
  {
    outputPeptide = inputPeptide;
    size_t start = outputPeptide.find_first_of("(");

    while (start != string::npos)
    {
      outputPeptide.erase(start, 1);
      size_t modStart = outputPeptide.find_first_of(",", start);
      size_t modEnd = outputPeptide.find_first_of(")", modStart);
      outputPeptide.erase(modStart, modEnd - modStart + 1);

      start = outputPeptide.find_first_of("(", modStart);
    }
  }
  // -------------------------------------------------------------------------
    bool PeptideSpectrumMatch::generateTheoreticalMasses(string &peptide,
                                                         string &ionNamesInclude,
                                                         MS2ScoringModel &inputIonTypes,
                                                         float prmOffset,
                                                         float srmOffset,
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
  bool PeptideSpectrumMatch::generateTheoreticalMasses(string &peptide,
                                                       string &ionNamesInclude,
                                                       MS2ScoringModel &inputIonTypes,
                                                       float prmOffset,
                                                       float srmOffset,
                                                       AAJumps &jumps,
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
            currMass = (prmMasses[ionIdx] + currIonFrag->massOffset
                + prmOffset) / currIonFrag->charge;
          }
          else
          {
            currMass = (srmMasses[ionIdx] + currIonFrag->massOffset
                + srmOffset) / currIonFrag->charge;
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
  void PeptideSpectrumMatch::getUnmodifiedPeptide(string &outputPeptide)
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
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::getModifications(vector<float> &modifications)
  {
    size_t start = m_annotation.find_first_of("(");

    while (start != string::npos)
    {
      size_t mod = m_annotation.find_first_of(",", start);
      size_t end = m_annotation.find_first_of(")", mod);
      float modification;
      if (sscanf(m_annotation.substr(mod + 1, end - mod - 1).c_str(),
                 "%f",
                 &modification))
      {
        modifications.push_back(modification);
      }
      else
      {
        WARN_MSG("Unable to convert " << m_annotation.substr(mod+1, end-mod-1) << " to float!");
        return false;
      }
      start = m_annotation.find_first_of("(", end);
    }
  }
  // -------------------------------------------------------------------------
  bool PeptideSpectrumMatch::insertModifications(string &unmodifiedPeptide,
                                                 vector<float> &modifications,
                                                 vector<unsigned int> &positions,
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
  bool PeptideSpectrumMatch::insertModifications(vector<float> &modifications,
                                                 vector<unsigned int> &positions,
                                                 string &outputPeptide)
  {
    return insertModifications(m_annotation,
                               modifications,
                               positions,
                               outputPeptide);
  }
  // -------------------------------------------------------------------------
  void PeptideSpectrumMatch::mapIons(vector<vector<int> > &outputMatches,
                                     std::tr1::unordered_map<string, int> &ionMap)
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

  void PeptideSpectrumMatch::inspectToSpecNets(string &inspect,
                                               string &specnets)
  {
    int inspectLength = inspect.length();
    bool modification = false;
    bool startMass = false;
    //specnets = "*.";

    const char * ptrInspect = inspect.c_str();
    char prevChar = 0x0;

    int i;
    //note: inspect_length - 3 means we ignore the trailing characters around annotation
    for (i = 2; i <= inspectLength - 3; i++)
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
      } else if ( (i == 2) && ('+' == currChar || '-' == currChar)) {
        startMass = true;
        if ('-' == currChar)
          specnets += "[-";
        else
          specnets += "[";
        
        // look for a mass in the sequence (modification)
      } else if ('+' == currChar || '-' == currChar)  {
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
      
      } else if (modification) {
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
      } else if (startMass)  {
        if (0 == currChar) // check for null
        {
          // do nothing
        }
        else if ('9' >= currChar || '.' == currChar)
        { //this is a number (or decimal)
          specnets += currChar;
        } else {
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
    ofs << m_unstrict_envelope_score << endl;

    return true;
  }

}
