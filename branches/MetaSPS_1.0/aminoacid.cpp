#include "aminoacid.h"

// Blocked cysteines
const float AAmasses[] = { 71.0371137870, 156.1011110260, 115.0269430310,
                           114.0429274460, 160.0306482000, 129.0425930950,
                           128.0585775100, 57.0214637230, 137.0589118610,
                           113.0840639790, 113.0840639790, 128.0949630160,
                           131.0404846050, 147.0684139150, 97.0527638510,
                           87.0320284090, 101.0476784730, 186.0793129520,
                           163.0633285370, 99.0684139150 };
// Blocked cysteines - NIPIA
//const float AAmasses[]  = {71.0371137870,156.1011110260,115.0269430310,114.0429274460,202.077598392,129.0425930950,128.0585775100,57.0214637230,137.0589118610,113.0840639790,113.0840639790,128.0949630160,131.0404846050,147.0684139150,97.0527638510,87.0320284090,101.0476784730,186.0793129520,163.0633285370,99.0684139150};
// Non-blocked cysteines
//const float AAmasses[]  = {71.0371137870,156.1011110260,115.0269430310,114.0429274460,103.009184477,129.0425930950,128.0585775100,57.0214637230,137.0589118610,113.0840639790,113.0840639790,128.0949630160,131.0404846050,147.0684139150,97.0527638510,87.0320284090,101.0476784730,186.0793129520,163.0633285370,99.0684139150};
const char AAletters[] = { 'A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I',
                           'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };
const char* AAnames[] = { "Alanine", "Arginine", "Aspartic acid", "Asparagine",
                          "Cysteine", "Glutamic acid", "Glutamine", "Glycine",
                          "Histidine", "Isoleucine", "Leucine", "Lysine",
                          "Methionine", "Phenylalanine", "Proline", "Serine",
                          "Threonine", "Tryptophan", "Tyrosine", "Valine" };
// Blocked cysteines + B = SILAC labeled Lysine (K+8)
//const unsigned int AAcount = 21;
//const float AAmasses[]  = {71.0371137870,156.1011110260,115.0269430310,114.0429274460,160.0306482000,129.0425930950,128.0585775100,57.0214637230,137.0589118610,113.0840639790,113.0840639790,128.0949630160,131.0404846050,147.0684139150,97.0527638510,87.0320284090,101.0476784730,186.0793129520,163.0633285370,99.0684139150,136.109162016};
//const char  AAletters[] = {'A','R','D','N','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','B'};
//const char* AAnames[]   = {"Alanine","Arginine","Aspartic acid","Asparagine","Cysteine","Glutamic acid","Glutamine","Glycine","Histidine","Isoleucine","Leucine","Lysine","Methionine","Phenylalanine","Proline","Serine","Threonine","Tryptophan","Tyrosine","Valine","silacK"};

vector<float> AAJumps::glbMasses;
vector<char> AAJumps::glbLetters;
const short AAJumps::NO_MODS = 0, AAJumps::USE_MODS = 1;
const double AAJumps::massHion = 1.0072763;
const double AAJumps::minAAmass = 57.0214637230;
const double AAJumps::massMH = 18.010564686 + 1.0072763;
const double AAJumps::massH2O = 18.010564686;
const double AAJumps::massNH3 = 17.026549101;
const double AAJumps::massCO = 27.994914622;

void getMasses(char *sequence, vector<float> &masses)
{
  masses.resize(strlen(sequence));
  int m;
  for (unsigned int i = 0; i < masses.size(); i++)
  {
    m = 0;
    while (m < 20 and AAletters[m] != sequence[i])
      m++;
    if (m >= 20)
      masses[i] = 0;
    else
      masses[i] = AAmasses[m];
  }
}

void getPepSeq(const char* sequence, char* destination)
{
  int m;
  int dest_index = 0;
  for (unsigned int i = 0; sequence[i] != 0; i++)
  {
    m = 0;
    while (m < 20 and AAletters[m] != sequence[i])
      m++;
    if (m < 20)
    {
      destination[dest_index] = sequence[i];
      dest_index++;
    }
  }
  destination[dest_index] = 0;
}

void getMassesCummulative(string sequence, vector<float> &masses, float offset) {
	const char * sequence_ptr = sequence.c_str();

	    int size = strlen(sequence_ptr);

	    masses.resize(size + 2);

	    float f;
	    int float_i;
	    char float_seq[25];
	    int spec_index = 0;
	    float total = offset;
	    int start = 2;

	    masses[spec_index] = total;
	      spec_index++;


	    if (sequence_ptr[1] != '.')
	    {
	      start = 0;
	    }

	    //note: we include the total peptide mass in this vector, which is a bogus
	    //"prefix mass", but is useful for generating suffix ions. Never use the
	    //last ion for any prefix ion calculations!
	    for (unsigned int i = start; i < size - start; i++)
	    {
	      // find global matching aa for curr sequence aa
	      int max_aa = AAcount;

	      // Square brackets: an unknown mass
	      if ('[' == sequence_ptr[i])
	      {
	        float_i = i + 1;
	        for (; sequence_ptr[i] != ']'; i++)
	          ; //iterate to closing bracket
	        strncpy(float_seq, &sequence_ptr[float_i], i - float_i);
	        float_seq[i - float_i] = '\0';

	        f = getFloat(float_seq);
	        //if (spec_index > 0) {
	        // if we have a modification at the very beginning of the peptide, we don't want it to appear as separate mass
	        //masses[spec_index - 1] = masses[spec_index - 1] + f;
	        //}
	        // Square brackets define a fixed mass for which we don't know the sequence
	        masses[spec_index] = f + total; // added
	        spec_index++; // added
	        total += f;

	        // round brackets: a modification
	      }
	      else if (sequence_ptr[i] == '(')
	      {
	        float_i = i + 1;
	        for (; sequence_ptr[i] != ')'; i++)
	          ; //iterate to closing bracket
	        strncpy(float_seq, &sequence_ptr[float_i], i - float_i);
	        float_seq[i - float_i] = '\0';

	        // find comma
	        int j;
	        for (j = 0; (j < i - float_i) && (float_seq[j] != ','); j++)
	          ;

	        // Parse aa's in sequence
	        float massAux = 0.0;
	        // Cycle through all amino acids before the comma
	        for (int k = 0; k < j; k++)
	        {
	          // skip spaces in the middle of the sequence
	          if (float_seq[k] == ' ' || float_seq[k] == '\t')
	            continue;
	          // initial aa index in aa table
	          int aa_index = 0;
	          // search for the aa in aa table
	          while (aa_index < max_aa and AAletters[aa_index] != float_seq[k])
	            aa_index++;
	          // if not found, issue message and skip
	          if (aa_index >= max_aa)
	          {
	            cerr << "Warning1: unknown aa! " << float_seq[k] << endl;
	            // otherwise, add it's mass
	          }
	          else
	          {
	            massAux += AAmasses[aa_index];
	          }
	        }
	        // Now process the modification value - number after the comma
	        // initialize modification mass
	        f = 0.0;
	        // if there is a comma, we have a modification
	        if (j < i - float_i)
	        {
	          // define space for value string
	          char float_seq2[25];
	          // get the string
	          strncpy(float_seq2, &float_seq[j + 1], i - float_i - j);
	          // add string terminator
	          float_seq2[i - float_i - j] = '\0';
	          // get the value
	          f = getFloat(float_seq2);
	        }
	        // add aa sequence mass value plus modification
	        masses[spec_index] = massAux + f + total; // added
	        total += massAux + f;
	        spec_index++;

	        //cout << "i2 " << i << endl;

	        // An aminoacid.
	      }
	      else
	      {

	        //cout << "i3 " << i << endl;

	        int aa_index = 0;
	        while (aa_index < max_aa and AAletters[aa_index] != sequence[i])
	          aa_index++;
	        if (aa_index >= max_aa)
	        {
	          cerr << "Warning2: unknown aa! " << sequence[i] << " in sequence " << sequence << endl;
	        }
	        else
	        {
	          masses[spec_index] = AAmasses[aa_index] + total;
	          total += AAmasses[aa_index];
	          spec_index++;
	        }
	      }
	    }
	    masses.resize(spec_index);
}


void getMassesCummulativeNoShift(char *sequence, vector<float> &masses) {
    masses.resize(strlen(sequence)+1);
    int m;
    int count = 1;
    float total = 0.0;
    masses[0] = 0.0;
    for(unsigned int i=0; sequence[i] != 0; i++) {
        m=0;
        while (m<20 and AAletters[m] != sequence[i]) m++;
        if(m<20) {
            masses[count] += AAmasses[m] + total;
            total += AAmasses[m];
            count ++;
        }
    }
    masses.resize(count);
}

void getMasses(vector<char> &sequence, vector<float> &masses) {
    masses.resize(sequence.size());
    int m;
    for(unsigned int i=0; i<masses.size(); i++) {
        m=0;
        while (m<20 and AAletters[m]!=sequence[i]) m++;
        if(m>=20) masses[i]=0; else masses[i]=AAmasses[m];
    }
}

//
// ************************************************************************************************
//  AAJumps::AAJumps(short maxJumpSize, float resolution, short useMods)
// ************************************************************************************************
//
AAJumps::AAJumps() {
	init(0);
}
AAJumps::AAJumps(short maxJumpSize,
                 float resolution,
                 float peakTol,
                 short useMods,
                 bool uniqueMasses)
{
  init(maxJumpSize, resolution, peakTol, useMods, uniqueMasses);
}

void AAJumps::init(short maxJumpSize,
                 float resolution,
                 float peakTol,
                 short useMods,
                 bool uniqueMasses)
{
  unsigned int aaIndex;
  if (glbMasses.empty())
  {
    glbMasses.resize(AAcount);
    glbLetters.resize(AAcount);
    for (aaIndex = 0; aaIndex < AAcount; aaIndex++)
      glbMasses[aaIndex] = AAmasses[aaIndex];
    for (aaIndex = 0; aaIndex < AAcount; aaIndex++)
      glbLetters[aaIndex] = AAletters[aaIndex];
  }
  getjumps(maxJumpSize, resolution, peakTol, useMods, uniqueMasses);
}

//
// ************************************************************************************************
//  void AAJumps::getjumps(short maxJumpSize, float resolution, short useMods)
// ************************************************************************************************
//
void AAJumps::getjumps(short maxJumpSize,
                       float resolution,
                       float peakTol,
                       short useMods,
                       bool uniqueMasses)
{
  unsigned int step = 0, aaIndex = 0, numMasses;

  if (maxJumpSize == 0)
  {
    masses.resize(1);
    masses[0] = 0;
    aaLetters.resize(1);
    aaLetters[0] = '0';
    return;
  }
  if (maxJumpSize < 0)
  {
    masses.resize(0);
    aaLetters.resize(1);
    return;
  }

  if (useMods == USE_MODS)
  {
    cerr << "(ERROR) AAjumps::getjumps: USE_MODS not implemented yet\n";
    masses.resize(0);
    exit(-1);
  }
  else
    modsUsed = useMods;

  if (refMasses.empty())
  {
    refMasses.resize(glbMasses.size());
    refLetters.resize(glbMasses.size());
    for (aaIndex = 0; aaIndex < glbMasses.size(); aaIndex++)
      refMasses[aaIndex] = glbMasses[aaIndex];
    for (aaIndex = 0; aaIndex < glbMasses.size(); aaIndex++)
      refLetters[aaIndex] = glbLetters[aaIndex];
  }

  int numMassesPerStep = refMasses.size(), totalNumMasses = numMassesPerStep;
  for (step = 2; step <= maxJumpSize; step++)
  {
    numMassesPerStep *= refMasses.size();
    totalNumMasses += numMassesPerStep;
  }
  masses.resize(totalNumMasses);
  //cout << "getjumps: masses.size() = " << masses.size() << ", " << endl;

  // This takes care of step 1
  for (aaIndex = 0; aaIndex < refMasses.size(); aaIndex++)
    masses[aaIndex] = refMasses[aaIndex];

  int prevStepStart = 0, // Start of previous step's jumps
      prevStepEnd = refMasses.size() - 1, // End of previous step's jumps
      prevIter, // Iterator for previous step's jumps
      curStepStart = refMasses.size(), // Start of current step's jumps
      curStepIter = curStepStart; // Iterates through current step's jumps

  for (step = 2; step <= maxJumpSize; step++)
  {
    //cerr << "Step = " << step << ": ["<<prevStepStart<<","<<prevStepEnd<<"] -> " << curStepStart << endl;
    for (prevIter = prevStepStart; prevIter <= prevStepEnd; prevIter++)
      for (aaIndex = 0; aaIndex < refMasses.size(); aaIndex++)
        masses[curStepIter++] = masses[prevIter] + refMasses[aaIndex];
    prevStepStart = curStepStart;
    prevStepEnd = curStepIter - 1;
    curStepStart = curStepIter;
  }

  if (maxJumpSize == 1)
  {
    if (uniqueMasses)
    {
      vector<unsigned int> idx;
      Utils::unique(masses, resolution, &idx);
      aaLetters.resize(idx.size());
      for (unsigned int aaIdx = 0; aaIdx < idx.size(); aaIdx++)
        aaLetters[aaIdx] = refLetters[idx[aaIdx]];
    }
    else
    {
      aaLetters.resize(refLetters.size());
      for (unsigned int aaIdx = 0; aaIdx < refLetters.size(); aaIdx++)
        aaLetters[aaIdx] = refLetters[aaIdx];
    }
  }
  else
  {
    aaLetters.resize(0);
    Utils::unique(masses, resolution);
  }

  //cout << "getjumps: masses.size() = " << masses.size() << endl;
  if (peakTol < 0)
    index.resize(0);
  else
    computeIndex(peakTol, resolution);
}

//
// ************************************************************************************************
//  void AAJumps::alljumps(short maxJumpMass, float resolution, short useMods)
// ************************************************************************************************
//
void AAJumps::alljumps(float maxJumpMass,
                       float resolution,
                       float peakTol,
                       short useMods)
{
  unsigned int vecSize = (unsigned int)ceil(maxJumpMass / resolution) + 1,
      jumpIdx, destIdx;
  unsigned int aaIdx;
  vector<bool> *jumpOk = new vector<bool> (vecSize);

  for (jumpIdx = 0; jumpIdx < vecSize; jumpIdx++)
    (*jumpOk)[jumpIdx] = false;
  getjumps(4, resolution, useMods); // Compute exact jump masses up to 4 amino acid jumps
  for (jumpIdx = 0; jumpIdx < masses.size(); jumpIdx++) // and initialize valid jumps
    (*jumpOk)[(int)round(masses[jumpIdx] / resolution)] = true;

  // Add new valid jumps from every valid jump
  unsigned int countOk = 0;
  for (jumpIdx = 0; jumpIdx < vecSize; jumpIdx++)
    if ((*jumpOk)[jumpIdx])
    {
      countOk++;
      for (aaIdx = 0; aaIdx < refMasses.size(); aaIdx++)
      {
        destIdx = jumpIdx + (unsigned int)round(refMasses[aaIdx] / resolution);
        if (destIdx < vecSize)
          (*jumpOk)[destIdx] = true;
      }
    }

  masses.resize(countOk);
  destIdx = 0;
  for (jumpIdx = 0; jumpIdx < vecSize; jumpIdx++)
    if ((*jumpOk)[jumpIdx])
      masses[destIdx++] = jumpIdx * resolution;
  if (peakTol < 0)
    index.resize(0);
  else
    computeIndex(peakTol, resolution, maxJumpMass);
  delete jumpOk;
}

//
// ************************************************************************************************
//  void AAJumps::addJumps(vector<float> &newJumps, vector<char> *newNames = 0)
//
//       Adds new jumps to masses (set of current jumps)
// ************************************************************************************************
//
void AAJumps::addJumps(vector<float> &newJumps, vector<char> *newNames)
{
  unsigned int numJumps = masses.size();
  masses.resize(numJumps + newJumps.size());
  aaLetters.resize(masses.size());
  index.resize(0);
  for (unsigned int i = 0; i < newJumps.size(); i++)
  {
    masses[numJumps + i] = newJumps[i];
    if (newNames != (vector<char> *)0 and i < newNames->size())
      aaLetters[numJumps + i] = (*newNames)[i];
    else
      aaLetters[numJumps + i] = 'X';
  }

  // Sort by increasing mass
  vector<pair<float, unsigned int> > massesIdx(masses.size());
  vector<float> oldMasses(masses.size());
  vector<char> oldLetters(masses.size());
  for (unsigned int i = 0; i < masses.size(); i++)
  {
    oldMasses[i] = masses[i];
    oldLetters[i] = aaLetters[i];
    massesIdx[i].first = masses[i];
    massesIdx[i].second = i;
  }
  sort(massesIdx.begin(), massesIdx.end());

  for (unsigned int i = 0; i < masses.size(); i++)
  {
    masses[i] = oldMasses[massesIdx[i].second];
    aaLetters[i] = oldLetters[massesIdx[i].second];
  }
}

//
// ************************************************************************************************
//  void AAJumps::loadJumps(char *filename)
//
//       Adds new jumps to masses (set of current jumps)
// ************************************************************************************************
//
bool AAJumps::loadJumps(const char *filename, bool setGlobal)
{
  InputParams aaMassParams;
  if (not aaMassParams.readParams(filename))
    return false;
  //cerr<<" -- Loaded "<<filename<<"\n"; cerr.flush();

  refMasses.resize(aaMassParams.size());
  refLetters.resize(aaMassParams.size());
  const char *paramName;
  for (unsigned int i = 0; i < aaMassParams.size(); i++)
  {
    refMasses[i] = aaMassParams.getValueDouble(i);
    paramName = aaMassParams.getParamName(i);
    if (strlen(paramName) > 0)
      refLetters[i] = paramName[0];
    else
      refLetters[i] = '-';
    //cerr<<" == Got -"<<refLetters[i]<<"- with mass "<<refMasses[i]<<"\n"; cerr.flush();
    //printf(" == Got - %c - with mass %.9f\n",refLetters[i],refMasses[i]);
  }
  //cerr<<" -- Done reading file\n"; cerr.flush();
  if (setGlobal)
  {
    glbMasses = refMasses;
    glbLetters = refLetters;
  }
  addJumps(refMasses, &refLetters);
  //cerr<<" -- Returning\n"; cerr.flush();
  return true;
}

//
// ************************************************************************************************
//  void AAJumps::forceDoubleSided()
//
//       Makes masses = [-masses; masses];
// ************************************************************************************************
//
void AAJumps::forceDoubleSided()
{
  vector<float> tmp(masses);

  masses.resize(2 * masses.size());
  index.resize(0);
  for (unsigned int i = 0; i < tmp.size(); i++)
  {
    masses[i] = -tmp[i];
    masses[i + tmp.size()] = tmp[i];
  }
}

//
// ************************************************************************************************
//  void AAJumps::forceJump(float mass)
//
//      Adds a given mass to the set of valid jumps
// ************************************************************************************************
//
void AAJumps::forceJumps(vector<float> newMasses)
{
  masses.reserve(masses.size() + newMasses.size());
  index.resize(0);
  for (unsigned int i = 0; i < newMasses.size(); i++)
    masses.push_back(newMasses[i]);
  sort(masses.begin(), masses.end());
}

//
// ************************************************************************************************
//  void AAJumps::forceTolerance(float tolerance, float resolution)
//
//      Every mass m in masses is replaced by a set of masses m+[-tolerance:resolution:tolerance]
// ************************************************************************************************
//
void AAJumps::forceTolerance(float tolerance, float resolution)
{
  float iter;

  int szRange = 0;
  for (iter = -tolerance; iter <= tolerance + 0.00001; iter += resolution)
    szRange++;
  if (szRange == 0)
    return;

  vector<float> tmp(masses);
  masses.resize(masses.size() * szRange);
  index.resize(0);
  for (unsigned int i = 0, curMass = 0; i < tmp.size(); i++)
    for (iter = -tolerance; iter <= tolerance + 0.00001; iter += resolution)
      masses[curMass++] = tmp[i] + iter;
}

//
// ************************************************************************************************
//  void AAJumps::removeHigherJumps(float highestMass)
//
//      Removes all jumps of mass higher than largestJump
// ************************************************************************************************
//
void AAJumps::removeHigherJumps(float largestJump)
{
  unsigned int idx;
  for (idx = 0; idx < masses.size() and masses[idx] <= largestJump; idx++)
    ;
  masses.resize(idx);
}

//
// ************************************************************************************************
//  bool AAJumps::isValid(float mass, float tolerance)
//
//      Test whether a given mass is a valid jump in the current set
// ************************************************************************************************
//
bool AAJumps::isValid(float mass, float tolerance)
{
  if (masses.size() == 0)
    return false;
  unsigned int curIdx = (masses.size()) / 2, high = masses.size(), low = 0;

  while (low < high)
  {
    if (fabs(masses[curIdx] - mass) <= tolerance + 0.000001)
      return true;
    if (masses[curIdx] < mass)
      low = curIdx + 1;
    else
      high = curIdx; // High is first non-checked position, low is lowest eligible position
    curIdx = (int)(high + low) / 2;
    //		szIdxJump = (int)abs(prev-curIdx)/2;   prev = curIdx;
    //		if(masses[curIdx]<mass) curIdx+=szIdxJump; else curIdx-=szIdxJump;
  }
  if (curIdx < masses.size() and fabs(masses[curIdx] - mass) <= tolerance
      + 0.0001)
    return true;
  //cerr<<"Closest was "<<masses[curIdx]<<endl;
  return false;
}

//
// ************************************************************************************************
//  unsigned int AAJumps::find(float mass, float tolerance, TwoValues<unsigned int> &idxBounds)
//
//      Find the indices of all jumps within tolerance of mass. On exit, idxBounds[0] is the
//    index of the first matching mass and idxBounds[1] is the index of the last match plus one.
// ************************************************************************************************
//
unsigned int AAJumps::find(float mass,
                           float tolerance,
                           TwoValues<unsigned int> &idxBounds)
{
  idxBounds.set(0, 0);
  if (masses.size() == 0)
    return 0;

  unsigned int curIdx = (masses.size()) / 2, high = masses.size(), low = 0;
  bool found = false;
  while (low < high and not found)
  {
    if (fabs(masses[curIdx] - mass) <= tolerance + 0.000001)
    {
      found = true;
      break;
    }
    if (masses[curIdx] < mass)
      low = curIdx + 1;
    else
      high = curIdx; // High is first non-checked position, low is lowest eligible position
    curIdx = (int)(high + low) / 2;
  }
  if (not found)
    return 0;

  for (; curIdx > 0 and fabs(masses[curIdx] - mass) <= tolerance + 0.000001; curIdx--)
    ; // Find first match
  if (fabs(masses[curIdx] - mass) > tolerance + 0.000001)
    curIdx++;
  idxBounds[0] = curIdx;
  for (; curIdx < masses.size() and masses[curIdx] - mass <= tolerance
      + 0.000001; curIdx++)
    ; // Find first index after last match
  idxBounds[1] = curIdx;

  return idxBounds[1] - idxBounds[0];
}

void AAJumps::computeIndex(float peakTol, float resolution, float maxMass)
{
  if (maxMass < 0)
    maxMass = masses[masses.size() - 1];
  maxMass += peakTol + resolution;
  unsigned int jStart, jEnd; // start/end indices of the jumps within tolerance of the current mass
  index.resize((unsigned int)ceil(maxMass / resolution));

  // Initialize jStart/jEnd
  for (jStart = 0; jStart < masses.size() and masses[jStart] <= 0; jStart++)
    ;
  if (jStart == masses.size())
  {
    index.resize(0);
    return;
  }
  jEnd = jStart;

  // Fill-out index
  float curMass;
  for (unsigned int idxMass = 0; idxMass < index.size(); idxMass++)
  {
    curMass = ((float)idxMass) * resolution;
    while (jStart < masses.size() and masses[jStart] < curMass - peakTol)
      jStart++;
    jEnd = max(jEnd, jStart);
    while (jEnd < masses.size() and masses[jEnd] <= curMass + peakTol)
      jEnd++;
    index[idxMass].set(jStart, jEnd);
  }
}

/** AAJumps::getPRMMasses
 * Calculates PRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
 * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
 * @param masses: 	Vector of masses to be modified.
 * @param offset: offset for prefix masses
 */
void AAJumps::getPRMMasses(string sequence, vector<float> &masses, float offset)
{
  const char * sequence_ptr = sequence.c_str();

  int size = strlen(sequence_ptr);

  masses.resize(size + 2);

  float f;
  int float_i;
  char float_seq[25];
  int spec_index = 0;
  float total = offset;
  int start = 2;

  if (sequence_ptr[1] != '.') {
    start = 0;
  }

  //note: we include the total peptide mass in this vector, which is a bogus
  //"prefix mass", but is useful for generating suffix ions. Never use the
  //last ion for any prefix ion calculations!
  for (unsigned int i = start; i < size - start; i++)
  {
    // find global matching aa for curr sequence aa
    int max_aa = glbLetters.size();

    // Square brackets: an unkwon mass
    if ('[' == sequence_ptr[i]) {
      float_i = i + 1;
      for (; sequence_ptr[i] != ']'; i++)
        ; //iterate to closing bracket
      strncpy(float_seq, &sequence_ptr[float_i], i - float_i);
      float_seq[i - float_i] = '\0';

      f = getFloat(float_seq);
      //if (spec_index > 0) {
        // if we have a modification at the very beginning of the peptide, we don't want it to appear as separate mass
        //masses[spec_index - 1] = masses[spec_index - 1] + f;
      //}
      // Square brackets define a fixed mass for which we don't know the sequence
      masses[spec_index] = f + total;   // added
      spec_index++;                     // added
      total += f;

    // round brackets: a modification
    } else if(sequence_ptr[i] == '(') {
      float_i = i + 1;
      for (; sequence_ptr[i] != ')'; i++)
        ; //iterate to closing bracket
      strncpy(float_seq, &sequence_ptr[float_i], i - float_i);
      float_seq[i - float_i] = '\0';

      // find comma
      int j;
      for(j = 0 ; (j < i - float_i) && (float_seq[j] != ',') ; j++)
        ;

      // Parse aa's in sequence
      float massAux = 0.0;
      // Cycle thru all aminoacids before the comma
      for(int k = 0 ; k < j ; k++) {
        // skip spaces in the middle of the sequence
        if(float_seq[k] == ' ' || float_seq[k] == '\t') continue;
        // initial aa index in aa table
        int aa_index = 0;
        // search for the aa in aa table
        while (aa_index < 20 and glbLetters[aa_index] != float_seq[k])
          aa_index++;
        // if not found, issue message and skip
        if (aa_index >= 20) {
          cerr << "Warning: unknown aa! " << float_seq[k] << endl;
        // otherwise, add it's mass
        } else {
          massAux += glbMasses[aa_index];
        }
      }
      // Now process the modification value - number after the comma
      // initialize modificatino mass
      f = 0.0;
      // if there is a comma, we have a modification
      if( j < i - float_i) {
        // define space for value string
        char float_seq2[25];
        // get the string
        strncpy(float_seq2, &float_seq[j+1], i - float_i - j);
        // add string terminator
        float_seq2[i - float_i - j] = '\0';
        // get the value
        f = getFloat(float_seq2);
      }
      // add aa sequence mass value plus modification
      masses[spec_index] = massAux + f + total;   // added
      total += massAux + f;
      spec_index++;

    // An aminoacid.
    } else {
      int aa_index = 0;
      while (aa_index < 20 and glbLetters[aa_index] != sequence[i])
        aa_index++;
      if (aa_index >= 20)
      {
        cerr << "Warning: unknown aa! " << sequence[i] << endl;
      }
      else
      {
        masses[spec_index] = glbMasses[aa_index] + total;
        total += glbMasses[aa_index];
        spec_index++;
      }
    }
  }
  masses.resize(spec_index);
}

/** AAJumps::getSRMMasses
 * Calculates SRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
 * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
 * @param masses: 	Vector of masses to be modified.
 * @param offset: offset for suffix masses
 */
void AAJumps::getSRMMasses(string sequence, vector<float> &masses, float offset)
{
  vector<float> prm_masses; //this is a bit of a cheat, but I'm lazy. Calculate SRM from PRM masses

  getPRMMasses(sequence, prm_masses);

  int prm_length = prm_masses.size();

  masses.resize(prm_length - 1);
  int i;
  int mass_i = 0;

  for (i = prm_length - 2; i >= 0; i--)
  {
    masses[mass_i] = prm_masses[prm_length - 1] - prm_masses[i] + offset; //total mass of peptide minus masses starting from end.
    mass_i++;
  }
}

/** AAJumps::getPRMandSRMMasses
 * Calculates SRM masses for input annotation in SpecNets format *.SEQ[-17]UENCE.*
 * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
 * @param prm_masses: 	Vector of prefix masses to be modified.
 * @param srm_masses: 	Vector of suffix masses to be modified.
 * @param peptide_mass: Summed mass of all amino acids. NOT THE SAME AS PRECURSOR MASS.
 *
 */
void AAJumps::getPRMandSRMMasses(string sequence,
                                 vector<float> &prm_masses,
                                 vector<float> &srm_masses,
                                 float &peptide_mass)
{
  getPRMMasses(sequence, prm_masses);

  int prm_length = prm_masses.size();

  peptide_mass = prm_masses[prm_length - 1];

  srm_masses.resize(prm_length - 1);
  int i;
  int mass_i = 0;

  for (i = prm_length - 2; i >= 0; i--)
  {
    srm_masses[mass_i] = prm_masses[prm_length - 1] - prm_masses[i]; //total mass of peptide minus masses starting from end.
    mass_i++;
  }
}

/** AAJumps::getPeptideMass
 * Calculates sum of all aa masses for input annotation in SpecNets format *.SEQ[-17]UENCE.* NOT PARENT MASS
 * @param sequence: Annotation of sequence in SpecNets format (i.e. brackets around modifications)
 */
double AAJumps::getPeptideMass(string &sequence)
{
  const char * sequence_ptr = sequence.c_str();

  int size = strlen(sequence_ptr);

  float f;
  int float_i;
  char float_seq[25];
  int spec_index = 0;
  double total = 0;

  for (unsigned int i = 2; i < size - 2; i++)
  {
    // find global matching aa for curr sequence aa
    int max_aa = glbLetters.size();

    if ('[' == sequence_ptr[i])
    {
      float_i = i + 1;
      for (; sequence_ptr[i] != ']'; i++)
        ; //iterate to closing bracket
      strncpy(float_seq, &sequence_ptr[float_i], i - float_i);
      float_seq[i - float_i] = '\0';

      f = getFloat(float_seq);
      total += f;
    }
    else
    {
      int aa_index = 0;
      while (aa_index < 20 and glbLetters[aa_index] != sequence[i])
        aa_index++;
      if (aa_index >= 20)
      {
        cerr << "Warning: unknown aa! " << sequence[i] << endl;
      }
      else
      {
        total += glbMasses[aa_index];
        spec_index++;
      }
    }
  }
  return total;
}

int AAJumps::getPeptideLength(string &sequence)
{
  const char * sequence_ptr = sequence.c_str();

  int size = strlen(sequence_ptr);
  int spec_index = 0;

  for (unsigned int i = 2; i < size - 2; i++)
  {
    if ('[' == sequence_ptr[i])
    {
      for (; sequence_ptr[i] != ']'; i++)
        ;
    }
    else
    {
      int aa_index = 0;
      while (aa_index < 20 and glbLetters[aa_index] != sequence[i])
        aa_index++;
      if (aa_index >= 20)
      {
        cerr << "Warning: unknown aa! " << sequence[i] << endl;
      }
      else
      {
        spec_index++;
      }
    }
  }
  return spec_index;
}
