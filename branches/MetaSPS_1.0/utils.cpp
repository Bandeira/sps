#include "utils.h"

//include proper form of mkdir function
#ifdef win32
#include <direct.h>
#endif

const char* CSV_SEP = ",";
const float PPM_FACTOR = 1.0 / 1000000.0;

////////////////////////////////////////////////////////////////////////////////
// wildcarToRegex(string <pattern to translate> 
//
// Translates a regular expression to a wildcard expression
////////////////////////////////////////////////////////////////////////////////
string wildcardToRegex(string pattern)
{
  // Make sure dots are not regex patterns
  replaceAll(pattern, ".", "\\.");
  // Replace the * to the regex .*
  replaceAll(pattern, "*", ".*");
  // The ? should be replace by a .
  replaceAll(pattern, "?", ".");
  // Returns the changed string
  return pattern;
}

void inspectToSpecNets(const string &inspectOriginal,
                                               string &specnets)
  {

    //specnets = "*.";
	
	if (inspectOriginal.length() == 0) {
	    return;
	}
	
    //check for trailing characters

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
  
////////////////////////////////////////////////////////////////////////////////
// relpaceAll(string <context>, string <from>, string <to>
//
// String pattern find and replace. Used by wildcardToRegex()
////////////////////////////////////////////////////////////////////////////////
string& replaceAll(string& context, const string& from, const string& to)
{
  // Where the parsing pointer is
  size_t lookHere = 0;
  // Where a string (from) was found
  size_t foundHere;
  // Cycle until the end of the string. FoundHere hold where the pattern was found
  while ((foundHere = context.find(from, lookHere)) < context.size())
  {
    // First part of the string copyed to an auxiliary string
    string aux = context.substr(0, foundHere);
    // Followed by the replacement pattern
    aux += to;
    // and finally the string last part
    aux += context.substr(foundHere + from.size());
    // assign it to context var for further processing
    context = aux;
    // the position where our pattern was found plus it's size is the next starting point
    lookHere = foundHere + to.size();
  }
  // Also returns the changed string for usability purposes
  return context;
}

string makeBracketsMods(string& peptide)
{

  if (peptide.length() == 0)
  {
    return string("");
  }

  if (peptide[0] == '[')
  {
    int length = 0;
    while (peptide[length] != ']')
    {
      ++length;
    }
    ++length;
    string modStr = peptide.substr(0, length);
    peptide.erase(0, length);
    peptide.insert(1, modStr);
  }

  list<int> insertpos;
  int addIdx = 0;
  for (int i = 0; i < peptide.length(); i++)
  {
    if (peptide[i] == '[')
    {
      insertpos.push_back(i - 1 + addIdx);
      addIdx++;
    }
  }

  string modPeptide(peptide);
  for (list<int>::iterator posIt = insertpos.begin(); posIt != insertpos.end(); posIt++)
  {
    modPeptide.insert(*posIt, "(");
  }
  modPeptide = replaceAll(modPeptide, "[", ",");
  modPeptide = replaceAll(modPeptide, "]", ")");
  return string(modPeptide);
}

////////////////////////////////////////////////////////////////////////////////
// regularExpressionMatch(string <regular expression>, string <string to match>
//
// Matches a string against a regular expression. Returns 0 on a succesful match.
////////////////////////////////////////////////////////////////////////////////
int regularExpressionMatch(string reg1, string str)
{
  // Get var space
  regex_t reg;

  // compile the RE given by reg1 parameter
  // Compile regex. If error, not 0 is returned
  if (regcomp(&reg, reg1.c_str(), REG_EXTENDED | REG_NOSUB | REG_ICASE) != 0)
  {
    return -1;
  }

  // Does the match. Returns 0 if it is a match
  return ((regexec(&reg, str.c_str(), 0, (regmatch_t *)NULL, 0)));
}
////////////////////////////////////////////////////////////////////////////////

bool mkdir_if_not_exist(const char* path)
{

#ifdef win32

  if (mkdir(path) == -1)
  {
    return false;
  }
  return true;
#else
  if (mkdir(path, 0777) == -1)
  {
    return false;
  }
  return true;

#endif
}

void parseArguments(char *argv[], int argc, set<string>& flags, map<string,
    string>& parsedArgs)
{
  parsedArgs.clear();
  int i = 0;
  while (i < argc)
  {
    string arg(argv[i]);
    if (flags.count(arg) > 0)
    {
      parsedArgs[arg] = "";
      if (i + 1 < argc)
      {
        string val(argv[i + 1]);
        if (flags.count(val) == 0)
        {
          parsedArgs[arg] = val;
          i += 2;
          continue;
        }
      }
    }
    i++;
  }
}

void getHistogramInfo(vector<float>& data, vector<TwoValues<float> >& binFreq)
{

  float last_bin = 0;
  for (int i = 0; i < binFreq.size(); i++)
  {
    float bin = binFreq[i][0];
    binFreq[i][1] = 0.0;
    for (int j = 0; j < data.size(); j++)
    {
      if (data[j] >= last_bin && (i == binFreq.size() - 1 || data[j] < bin))
        binFreq[i][1] += 1.0;
    }
    last_bin = bin;
  }
}

/**
 * Writes a vector of strings to a file in binary format
 * @param fp file ptr
 * @param strings vector of strings
 * @return true if written successfully, false otherwise
 */
bool writeStringsToBinaryStream(FILE* fp, vector<string>& strings)
{
  if (fp == 0)
  {
    return false;
  }

  unsigned int count;
  unsigned int numStrings = strings.size();

  count = fwrite(&numStrings, sizeof(unsigned int), 1, fp);
  if (count == 0)
  {
    return false;
  }
  unsigned int strLen;

  for (unsigned int i = 0; i < numStrings; i++)
  {
    strLen = strlen(strings[i].c_str());

    count = fwrite(&strLen, sizeof(unsigned int), 1, fp);
    if (strLen == 0)
    {
      continue;
    }
    else if (count == 0)
    {
      return false;
    }

    count = fwrite(strings[i].c_str(), sizeof(char), strLen, fp);
    if (count == 0)
    {
      return false;
    }
  }
  return true;
}

/**
 * Reads a vector of strings from a file in binary format
 * @param fp file ptr
 * @param strings vector of strings
 * @return true if read successfully, false otherwise
 */
bool readStringsFromBinaryStream(FILE* fp, vector<string>& strings)
{
  if (fp == 0)
  {
    return false;
  }

  unsigned int count;
  unsigned int numStrings = 0;

  count = fread(&numStrings, sizeof(unsigned int), 1, fp);
  if (count == 0)
  {
    return false;
  }
  strings.resize(numStrings);
  unsigned int strLen;

  for (unsigned int i = 0; i < numStrings; i++)
  {

    count = fread(&strLen, sizeof(unsigned int), 1, fp);
    if (count == 0)
    {
      return false;
    }

    if (strLen == 0)
    {
      strings[i] = "";
      continue;
    }

    char* labelArr = (char*)malloc(strLen + 1);

    count = fread(labelArr, sizeof(char), strLen, fp);
    if (count == 0)
    {
      free(labelArr);
      return false;
    }
    labelArr[strLen] = 0;

    strings[i] = labelArr;

    free(labelArr);
  }
  return true;
}

float getResolution(float peakTol)
{
  float divisor = peakTol;
  while (divisor < 1)
    divisor *= 10.0;
  return peakTol / divisor;
}

float getFloat(const char* str)
{
  return (float)atof(str);
}

int floatToInt(float x)
{
  if (x > 0)
  {
    return (int)(round(x) + 0.01);
  }
  else if (x == 0)
  {
    return 0;
  }
  else
  {
    return (int)(round(x) - 0.01);
  }
}

int getInt(const char* str)
{
  istringstream ins(str);
  int f;
  ins >> f;
  return f;
}

/**
 * Depreciated. See MZRange::EqualWithinRange
 */
bool isEqual(float f1, float f2, float range)
{
  return (f1 >= f2 - range) && (f1 <= f2 + range);
}

string parseInt(int x, int equalizeLength)
{
  ostringstream out;
  out << x;
  string formatted = out.str();
  //cout << "formatted: " << formatted << endl; cout.flush();
  if (equalizeLength <= formatted.length())
    return string(formatted);

  string pad = "0";
  string temp = "";
  for (int diff = equalizeLength - formatted.length(); diff != 0; diff--)
    temp += pad;

  temp += formatted;
  return string(temp);
}

string parseFloat(float x, int prec)
{
  ostringstream out;
  out << fixed << setprecision(prec) << x;
  return out.str();
}

bool splitText(const char* text, list<string>& vec, const char* delim)
{
  string str_res;
  vec.clear();
  unsigned int textLength = strlen(text);
  char* textCpy = (char*)malloc(textLength + 1);
  if (strncpy(textCpy, text, strlen(text)) == NULL)
  {
    cerr << "ERROR SPLITTING TEXT\n";
    return false;
  }
  textCpy[textLength] = '\0';
  const char* result = strtok(textCpy, delim);
  while (result != NULL)
  {
    str_res = result;
    vec.push_back(str_res);
    result = strtok(NULL, delim);
  }
  free(textCpy);
  return true;
}

bool splitText(const char* text, vector<string>& vec, const char* delim)
{
  string str_res;
  vec.resize(0);
  unsigned int textLength = strlen(text);
  char* textCpy = (char*)malloc(textLength + 1);
  if (strncpy(textCpy, text, strlen(text)) == NULL)
  {
    cerr << "ERROR SPLITTING TEXT\n";
    return false;
  }
  textCpy[textLength] = '\0';
  const char* result = strtok(textCpy, delim);
  while (result != NULL)
  {
    str_res = result;
    vec.push_back(str_res);
    result = strtok(NULL, delim);
  }
  free(textCpy);
  return true;
}

bool getdir(string dir, list<string> &files)
{
  files.clear();
  DIR *dp;
  struct dirent *dirp;
  if ((dp = opendir(dir.c_str())) == NULL)
    return false;

  while ((dirp = readdir(dp)) != NULL)
  {
    files.push_back(string(dirp->d_name));
  }
  closedir(dp);
  files.sort();
  return true;
}

ProgressDisplay::ProgressDisplay(ostream& output, long e, const char* msg)
{
  end = e;
  output << msg;
  output.flush();
  len = 0;
}

void ProgressDisplay::showProgress(ostream& output, long pos)
{
  char back[12];
  const char* perc =
      parseFloat(100.0 * ((float)pos) / ((float)(end)), 3).c_str();
  for (int x = 0; x < 11; x++)
    back[x] = '\b';
  back[len] = 0;
  output << back << perc;
  output.flush();
  len = strlen(perc);
}

void ProgressDisplay::clear(ostream& output)
{
  output << endl;
}

struct value_index
{ // Used in Utils::unique() to recover the indices of the sorted values
  float value;
  unsigned int index;
};
bool cmp_value_index(value_index vi1, value_index vi2)
{
  return vi1.value < vi2.value;
}

vector<float> &Utils::unique(vector<float> &v,
                             float r,
                             vector<unsigned int> *idx)
{
  if (v.size() == 0)
    return v;

  if (idx)
  {
    vector<value_index> tmpV(v.size());
    idx->resize(v.size());
    for (unsigned int pivot = 0; pivot < v.size(); pivot++)
    {
      tmpV[pivot].value = v[pivot];
      tmpV[pivot].index = pivot;
    }
    sort(tmpV.begin(), tmpV.end(), cmp_value_index);
    for (unsigned int pivot = 0; pivot < v.size(); pivot++)
    {
      v[pivot] = tmpV[pivot].value;
      (*idx)[pivot] = tmpV[pivot].index;
    }
  }
  else
    sort(v.begin(), v.end());
  int i;
  vector<float> newV(v);
  for (i = 0; i < newV.size(); i++)
    newV[i] = round(newV[i] / r);

  // Count number of different rounded elements
  int numElems = 1, lastElem = (int)newV[0];
  for (i = 1; i < newV.size(); i++)
    if ((int)newV[i] > lastElem)
    {
      numElems++;
      lastElem = (int)newV[i];
    }

  // Store rounded elements
  v.resize(numElems);
  if (idx)
    idx->resize(numElems);
  numElems = 0;
  v[0] = newV[0];
  for (i = 1; i < newV.size(); i++)
    if ((int)newV[i] > (int)v[numElems])
    {
      numElems++;
      v[numElems] = newV[i];
      if (idx)
        (*idx)[numElems] = (*idx)[i]; // Note that i >= numElems
    }

  // Set resolution to what was defined by parameter r
  for (i = 0; i < v.size(); i++)
    v[i] = r * v[i];

  return v;
}

int Utils::intersect(vector<float> &v1,
                     vector<float> &v2,
                     vector<float> &putHere,
                     float tolerance)
{
  return 0;
}

int Utils::save_tcb(char *filename, vector<TwoValues<float> > data)
{
  FILE *fp;
  unsigned int numEntries = data.size();
  float curData[2];
  unsigned int i, p;

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries, sizeof(int), 1, fp); // Number of entries in the file

  for (i = 0; i < numEntries; i++)
  {
    curData[0] = data[i][0];
    curData[1] = data[i][1];
    fwrite(curData, sizeof(float), 2, fp);
  }

  fclose(fp);
  return 1;
}

int Utils::save_tcb_list(char *filename, list<TwoValues<float> > data)
{
  FILE *fp;
  unsigned int numEntries = data.size();
  float curData[2];
  unsigned int i, p;

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries, sizeof(int), 1, fp); // Number of entries in the file

  for (list<TwoValues<float> >::iterator iter = data.begin(); iter
      != data.end(); iter++)
  {
    curData[0] = (*iter)[0];
    curData[1] = (*iter)[1];
    fwrite(curData, sizeof(float), 2, fp);
  }

  fclose(fp);
  return 1;
}

// Converts a set of values to an indicator boolean vector with true whenever a value is present in values
unsigned int Utils::list_to_vector_bool(vector<float> &values,
                                        vector<bool> &output,
                                        float tolerance,
                                        float multFactor)
{
  unsigned int offset = (int)ceil(-(values[0] - tolerance) * multFactor);
  int tolRange = (int)round(tolerance * multFactor);

  output.resize(offset + (int)ceil((values[values.size() - 1] + tolerance)
      * multFactor) + 1); // resize and initialize
  for (unsigned int i = 0; i < output.size(); i++)
    output[i] = false;

  for (unsigned int vIdx = 0; vIdx < values.size(); vIdx++)
  { // Set value positions to true
    int intValue = offset + (int)round(values[vIdx] * multFactor);
    for (int tolOffset = -tolRange; tolOffset <= tolRange; tolOffset++)
      output[intValue + tolOffset] = true;
  }
  return offset;
}

// Computes a table of binomial probabilities for n=1..maxN, k=1..n; probs[i][j]=p(k=j|n=j+1)
void Utils::binomial(unsigned short maxN,
                     float p,
                     vector<vector<float> > &probs)
{
  unsigned short n, k, idxN;

  probs.resize(maxN);
  if (maxN == 0)
    return;
  for (idxN = 0; idxN < maxN; idxN++)
    probs[idxN].resize(maxN + 1);

  // Compute the n-choose-k values
  for (idxN = 0; idxN < maxN; idxN++)
    probs[idxN][0] = 1;
  for (k = 1; k <= maxN; k++)
    for (idxN = 0; idxN < maxN; idxN++)
    {
      n = idxN + 1; // index=idxN corresponds to n=idxN+1 trials in the binomial distribution
      if (k > n)
      {
        probs[idxN][k] = 0;
        continue;
      }
      if (k == n)
        probs[idxN][k] = 1;
      else
        probs[idxN][k] = probs[idxN - 1][k] * n / (n - k);
    }

  // Multiply by the p^k * (1-p)^(n-k)
  vector<TwoValues<float> > tmpProbs(maxN + 1); // pos[k] contains ( p^k , (1-p)^k )
  tmpProbs[0].set(1, 1);
  tmpProbs[1].set(p, 1 - p);
  for (k = 2; k <= maxN; k++)
    tmpProbs[k].set(p * tmpProbs[k - 1][0], (1 - p) * tmpProbs[k - 1][1]);
  for (k = 0; k <= maxN; k++)
    for (idxN = 0; idxN < maxN; idxN++)
      if (k <= idxN + 1)
        probs[idxN][k] *= tmpProbs[k][0] * tmpProbs[idxN + 1 - k][1];
}

double Utils::gaussiancdf(double x, double mean, double stddev)
{

  if (fabs(stddev) < 0.000001)
    return -1;
  x -= mean;
  x /= stddev;

  //  Start code obtained from http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
  //	double N(const double x) {
  const double b1 = 0.319381530;
  const double b2 = -0.356563782;
  const double b3 = 1.781477937;
  const double b4 = -1.821255978;
  const double b5 = 1.330274429;
  const double p = 0.2316419;
  const double c = 0.39894228;

  if (x >= 0.0)
  {
    double t = 1.0 / (1.0 + p * x);
    return (1.0 - c * exp(-x * x / 2.0) * t * (t * (t
        * (t * (t * b5 + b4) + b3) + b2) + b1));
  }
  else
  {
    double t = 1.0 / (1.0 - p * x);
    return (c * exp(-x * x / 2.0) * t * (t
        * (t * (t * (t * b5 + b4) + b3) + b2) + b1));
  }
  //	}
  //  End code obtained from http://www.sitmo.com/doc/Calculating_the_Cumulative_Normal_Distribution
}

// Computes scores[i][j]=log(probsSignal[i][j]/probsNoise[i][j]) for all probsNoise[i][j]>0
void Utils::logscores(vector<vector<float> > &probsSignal,
                      vector<vector<float> > &probsNoise,
                      vector<vector<float> > &scores)
{
  unsigned int i, j;
  if (probsSignal.size() != probsNoise.size())
  {
    cerr
        << "ERROR in Utils::logscores(): probs vector dimensions do not match!\n";
    return;
  }
  scores.resize(probsSignal.size());
  for (i = 0; i < probsSignal.size(); i++)
    if (probsSignal[i].size() != probsNoise[i].size())
    {
      cerr
          << "ERROR in Utils::logscores(): probs vector dimensions do not match!\n";
      scores.resize(0);
      return;
    }
    else
    {
      scores[i].resize(probsSignal[i].size());
      for (j = 0; j < probsSignal[i].size(); j++)
        if (probsNoise[i][j] > 0)
          scores[i][j] = log(probsSignal[i][j] / probsNoise[i][j]);
        else
          scores[i][j] = 0;
    }
}

/*
 *  LoadGaussianParams - Load and merges a set of files with means and variances
 *     of different samples from the same random variables.
 *
 *  meansIndexFN - Text file where each line contains the name of a file with values for the means
 *  varsIndexFN  - Text file where each line contains the name of a file with values for the variances
 *  params       - Output Gaussian parameters (means and standard devs) are returned in this variable
 */
bool LoadGaussianParams(const char *meansIndexFN,
                        const char *varsIndexFN,
                        vector<TwoValues<float> > &params)
{
  BufferedLineReader blrMeans, blrVars;
  unsigned int pivot, fIdx, szParams, numFiles;
  params.resize(0);
  if (blrMeans.Load(meansIndexFN) <= 0)
  {
    cerr << "ERROR reading " << meansIndexFN << "!\n";
    return false;
  }
  numFiles = blrMeans.size();
  if (blrVars.Load(varsIndexFN) <= 0)
  {
    cerr << "ERROR reading " << varsIndexFN << "!\n";
    return false;
  }
  if (numFiles != blrVars.size())
  {
    cerr << "ERROR: different number of files for means and variances ( "
        << numFiles << " / " << blrVars.size() << " )\n";
    return false;
  }

  // Read first means/vars
  vector<vector<float> > bufferM, bufferV;
  vector<float> counts;
  if (Load_binArray(blrMeans.getline(0), bufferM) <= 0)
  {
    cerr << "ERROR reading " << blrMeans.getline(0) << "!\n";
    return false;
  }
  if (bufferM[0].size() != 2)
  {
    cerr << "ERROR: incorrect file format for " << blrMeans.getline(0) << "!\n";
    return false;
  }
  if (Load_binArray(blrVars.getline(0), bufferV) <= 0)
  {
    cerr << "ERROR reading " << blrVars.getline(0) << "!\n";
    return false;
  }
  if (bufferV.size() != bufferM.size() or bufferV[0].size() != 1)
  {
    cerr << "ERROR: incorrect file format for " << blrVars.getline(0) << "!\n";
    return false;
  }
  szParams = bufferM.size();
  params.resize(szParams);
  counts.resize(szParams);
  for (pivot = 0; pivot < szParams; pivot++)
  {
    params[pivot][0] = bufferM[pivot][0];
    params[pivot][1] = bufferV[pivot][0];
    counts[pivot] = bufferM[pivot][1];
  }

  // Read remaining means/vars and update existing values
  float ratio1, ratio2;
  for (fIdx = 1; fIdx < numFiles; fIdx++)
  {
    if (strlen(blrMeans.getline(fIdx)) == 0)
      continue;
    if (Load_binArray(blrMeans.getline(fIdx), bufferM) <= 0)
    {
      cerr << "ERROR reading " << blrMeans.getline(fIdx) << "!\n";
      params.resize(0);
      return false;
    }
    if (bufferM.size() != szParams or bufferM[0].size() != 2)
    {
      cerr << "ERROR: incorrect file format for " << blrMeans.getline(fIdx)
          << "!\n";
      params.resize(0);
      return false;
    }
    if (Load_binArray(blrVars.getline(fIdx), bufferV) <= 0)
    {
      cerr << "ERROR reading " << blrVars.getline(fIdx) << "!\n";
      params.resize(0);
      return false;
    }
    if (bufferV.size() != szParams or bufferV[0].size() != 1)
    {
      cerr << "ERROR: incorrect file format for " << blrVars.getline(fIdx)
          << "!\n";
      params.resize(0);
      return false;
    }

    for (pivot = 0; pivot < szParams; pivot++)
    {
      if (bufferM[pivot][1] > 0.1)
      { // Avoid computations when there are no new observations
        ratio1 = counts[pivot] / (counts[pivot] + bufferM[pivot][1]);
        ratio2 = 1 - ratio1;
        params[pivot][0] = params[pivot][0] * ratio1 + bufferM[pivot][0]
            * ratio2;
        params[pivot][1] = params[pivot][1] * ratio1 + bufferV[pivot][0]
            * ratio2;
        counts[pivot] += bufferM[pivot][1];
      }
    }
  }

  // Compute final Gaussian parameters
  for (pivot = 0; pivot < szParams; pivot++)
    params[pivot][1] = sqrt(params[pivot][1] - params[pivot][0]
        * params[pivot][0]);
  return true;
}

short BufferedLineReader::Load(const char *filename)
{
  reset();
  FILE *input = fopen(filename, "r");
  if (!input or ferror(input))
  {
    cerr << "ERROR opening " << filename << "!\n";
    return -1;
  }
  if (fseek(input, 0, SEEK_END))
  {
    cerr << "ERROR finding the end of the file (fseek)!\n";
    return -1;
  }
  unsigned int numBytes = ftell(input);
  if (numBytes == 0)
  {
    cerr << "ERROR file size = " << numBytes << "?\n";
    return -1;
  }
  if (fseek(input, 0, SEEK_SET))
  {
    cerr << "ERROR finding the start of the file (fseek)!\n";
    return -1;
  }
  numBytes++;

  lines = (char *)malloc(numBytes);
  lines[numBytes - 1] = '\0';
  unsigned int readBytes = fread(lines, 1, numBytes - 1, input);
  if (readBytes != numBytes - 1)
  {
    cerr << "ERROR reading " << filename << " - only got " << readBytes
        << " bytes instead of " << numBytes - 1 << "?\n";
    return -1;
  }

  unsigned int numLines = 0, pivot;
  for (pivot = 0; pivot < numBytes; pivot++)
  {
    if (lines[pivot] == '\n')
      numLines++;
    else if (lines[pivot] == '\r')
      lines[pivot] = '\0';
  }
  linesIdx.resize(numLines + 1); // +1 because last line may not have a '\n'
  linesIdx[0] = 0;
  numLines = 1;
  for (pivot = 0; pivot < numBytes; pivot++)
  {
    if (lines[pivot] == '\n')
    {
      linesIdx[numLines++] = pivot + 1;
      lines[pivot] = '\0';
    }
  }
  fclose(input);
  return 1;
}

// Hack to force the instantiation of a few variants of the function template
void binArray_instantiator()
{
  vector<vector<float> > dataF;
  Save_binArray("", dataF);
  Load_binArray("", dataF);
  Load_binArray_multiple("", dataF);
  vector<vector<int> > dataI;
  Save_binArray("", dataI);
  Load_binArray("", dataI);
  Load_binArray_multiple("", dataI);
  vector<vector<unsigned int> > dataUI;
  Save_binArray("", dataUI);
  Load_binArray("", dataUI);
  Load_binArray_multiple("", dataUI);
  vector<vector<unsigned short> > dataUS;
  Save_binArray("", dataUS);
  Load_binArray("", dataUS);
  Load_binArray_multiple("", dataUS);
  vector<vector<double> > dataD;
  Save_binArray("", dataD);
  Load_binArray("", dataD);
  Load_binArray_multiple("", dataD);
  vector<vector<short> > dataS;
  Save_binArray("", dataS);
  Load_binArray("", dataS);
  Load_binArray_multiple("", dataS);

  vector<list<float> > dataFL;
  Save_binListArray<float, list<float> , list<float>::iterator> ("", dataFL);
  Save_binListArray<float, vector<float> , vector<float>::iterator> ("", dataF);
  Load_binListArray<float, list<float> , list<float>::iterator> ("", dataFL);
  vector<list<int> > dataIL;
  Save_binListArray<int, list<int> , list<int>::iterator> ("", dataIL);
  Save_binListArray<int, vector<int> , vector<int>::iterator> ("", dataI);
  Load_binListArray<int, list<int> , list<int>::iterator> ("", dataIL);
  vector<list<double> > dataDL;
  Save_binListArray<double, list<double> , list<double>::iterator> ("", dataDL);
  Save_binListArray<double, vector<double> , vector<double>::iterator> ("",
                                                                        dataD);
  Load_binListArray<double, list<double> , list<double>::iterator> ("", dataDL);
  vector<list<short> > dataSL;
  Save_binListArray<short, list<short> , list<short>::iterator> ("", dataSL);
  Save_binListArray<short, vector<short> , vector<short>::iterator> ("", dataS);
  Load_binListArray<short, list<short> , list<short>::iterator> ("", dataSL);

  vector<float> dataF1d;
  Save_binArray<float> ("", dataF1d);
  Load_binArray<float> ("", dataF1d);
  vector<int> dataI1d;
  Save_binArray<int> ("", dataI1d);
  Load_binArray<int> ("", dataI1d);
  vector<unsigned int> dataUI1d;
  Save_binArray<unsigned int> ("", dataUI1d);
  Load_binArray<unsigned int> ("", dataUI1d);
  vector<double> dataD1d;
  Save_binArray<double> ("", dataD1d);
  Load_binArray<double> ("", dataD1d);
  vector<short> dataS1d;
  Save_binArray<short> ("", dataS1d);
  Load_binArray<short> ("", dataS1d);

  vector<TwoValues<short> > dataStv;
  Save_binArray("", dataStv);
  vector<TwoValues<int> > dataItv;
  Save_binArray("", dataItv);
  vector<TwoValues<unsigned int> > dataUItv;
  Save_binArray("", dataUItv);
  vector<TwoValues<float> > dataFtv;
  Save_binArray("", dataFtv);
  vector<TwoValues<double> > dataDtv;
  Save_binArray("", dataDtv);
  list<TwoValues<short> > dataSLtv;
  Save_binArray("", dataSLtv);
  list<TwoValues<float> > dataFLtv;
  Save_binArray("", dataFLtv);
  list<TwoValues<double> > dataDLtv;
  Save_binArray("", dataDLtv);

  list<vector<int> > dataIV;
  Save_binListArray<int, vector<int> , vector<int>::iterator> ("", dataIV);
  Save_binArray<int> ("", dataIV);
  list<vector<unsigned int> > dataUIV;
  Save_binListArray<unsigned int, vector<unsigned int> ,
      vector<unsigned int>::iterator> ("", dataUIV);
  Save_binArray<unsigned int> ("", dataUIV);
  list<vector<float> > dataFV;
  Save_binListArray<float, vector<float> , vector<float>::iterator> ("", dataFV);
  Save_binArray<float> ("", dataFV);
  vector<vector<unsigned int> > dataUIVV;
  Save_binListArray<unsigned int, vector<unsigned int> ,
      vector<unsigned int>::iterator> ("", dataUIVV);
  Save_binArray<unsigned int> ("", dataUIVV);
  Load_binListArray<unsigned int, vector<unsigned int> ,
      vector<unsigned int>::iterator> ("", dataUIVV);
  vector<list<unsigned int> > dataUIVL;
  Save_binListArray<unsigned int, list<unsigned int> ,
      list<unsigned int>::iterator> ("", dataUIVL);
  Load_binListArray<unsigned int, list<unsigned int> ,
      list<unsigned int>::iterator> ("", dataUIVL);
}

bool loadCsvToTxt(const char* infilecsv,
                  const char* outfiletxt,
                  float threshold)
{
  BufferedLineReader buf;
  if (!buf.Load(infilecsv))
    return false;
  char* line = buf.getline(buf.size() - 1);
  int end = buf.size() - 1;
  if (strlen(line) == 0)
  {
    line = buf.getline(buf.size() - 2);
    end--;
  }
  int size = getInt(line);
  int lastIndex = -1;
  string s1 = "N[0.984]";
  string s2 = "M[16]";
  string s3 = "Q[1]";
  string res;
  vector<string> seqs(size + 1);
  for (int i = 1; i <= end; i++)
  {
    line = buf.getline(i);
    char* result = strtok(line, ",");
    int index = getInt(result);

    if (index == lastIndex)
      continue;
    else
      lastIndex = index;

    for (int c = 0; c < 1; c++)
      result = strtok(NULL, ",");
    res = result;

    result = strtok(NULL, ",");
    //cout << index << "," << getFloat(result) << "\n";
    if (getFloat(result) > threshold)
      continue;
    //cout << "outputting: " << res << "\n";
    int idx;
    while (string::npos != (idx = res.find("n")))
    {
      res = res.erase(idx, 1);
      res = res.insert(idx, s1, 0, s1.length());
    }
    while (string::npos != (idx = res.find("m")))
    {
      res = res.erase(idx, 1);
      res = res.insert(idx, s2, 0, s2.length());
    }
    while (string::npos != (idx = res.find("q")))
    {
      res = res.erase(idx, 1);
      res = res.insert(idx, s3, 0, s3.length());
    }
    seqs[index] = res;
  }
  FILE* output = fopen(outfiletxt, "w");
  for (int i = 0; i < seqs.size(); i++)
  {
    if (seqs[i].length() > 0)
      fprintf(output, "%s", seqs[i].c_str());
    fprintf(output, "\n");
  }
  fprintf(output, "\n");
  fclose(output);

}

/*
 template<class T,class T2> int getIndexData(vector<vector<T> > &data, unsigned int **index, T **dataArray) {
 unsigned int numEntries1 = data.size(), numEntries2=0;
 (*index) = (unsigned int*)malloc(sizeof(unsigned int)*numEntries1);
 for(unsigned int i=0;i<numEntries1;i++) { (*index)[i]=data[i].size(); numEntries2+=(*index)[i]; }
 (*dataArray) = (T*)malloc(sizeof(T)*numEntries2);

 int dataArrayIdx=0;
 for(unsigned int i=0;i<numEntries1;i++)
 for(unsigned int j=0;j<data[i].size();j++)
 (*dataArray)[dataArrayIdx++] = data[i][j];
 return numEntries2;
 }

 template<class T,class T2> int getIndexData(vector<list<T> > &data, unsigned int **index, T **dataArray) {
 unsigned int numEntries1 = data.size(), numEntries2=0;
 (*index) = (unsigned int*)malloc(sizeof(unsigned int)*numEntries1);
 for(unsigned int i=0;i<numEntries1;i++) { (*index)[i]=data[i].size(); numEntries2+=(*index)[i]; }
 (*dataArray) = (T*)malloc(sizeof(T)*numEntries2);

 int dataArrayIdx=0;   T2 iter;  //list<T>::iterator iter;
 for(unsigned int i=0;i<numEntries1;i++)
 for(iter=data[i].begin(); iter!=data[i].end(); iter++)
 (*dataArray)[dataArrayIdx++] = *iter;
 return numEntries2;
 }

 template<class T,class T2,class T3> int getIndexData2(vector<T2> &data, unsigned int **index, T **dataArray) {
 unsigned int numEntries1 = data.size(), numEntries2=0;
 (*index) = (unsigned int*)malloc(sizeof(unsigned int)*numEntries1);
 for(unsigned int i=0;i<numEntries1;i++) { (*index)[i]=data[i].size(); numEntries2+=(*index)[i]; }
 (*dataArray) = (T*)malloc(sizeof(T)*numEntries2);

 int dataArrayIdx=0;   T3 iter;  //list<T>::iterator iter;
 for(unsigned int i=0;i<numEntries1;i++)
 for(iter=data[i].begin(); iter!=data[i].end(); iter++)
 (*dataArray)[dataArrayIdx++] = *iter;
 return numEntries2;
 }

 //
 //  Example usage: T1 = int, T2 = list<int>, T3 = list<int>::iterator
 //
 template<class T1,class T2,class T3> int Load_binListArray(char *filename, vector<T2> &data) {
 FILE *fp;   unsigned int numLists=0, numElems=0;
 unsigned int *index;
 unsigned int i,streamIdx; int n;

 fp = fopen(filename,"r");
 if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return -1; }

 // Read number of lists
 n=fread(&numLists,sizeof(unsigned int),1,fp); if(numLists<=0 | n!=1) { cerr << "Invalid number of lists ("<<numLists<<")\n" ;return -1; }

 // Read index
 index = (unsigned int*)malloc(sizeof(unsigned int)*numLists);
 n=fread(index,sizeof(unsigned int),numLists,fp);  if(n!=numLists)  { cerr << "Invalid index!\n" ; free(index); return -1; }
 numElems=0; for(i=0;i<numLists;i++) numElems+=index[i];

 // Read data
 T1 *streamData = (T1 *)malloc(sizeof(T1)*numElems);
 n=fread(streamData,sizeof(T1),numElems,fp);    fclose(fp);
 if(n!=(int)numElems)  { cerr << "Error reading binary list array: not enough elements in the file.\n"; free(streamData); free(index); return -1; }

 // Convert back to C++ structure
 data.resize(numLists);   T3 iter;    streamIdx=0;
 for(i=0; i<numLists; i++){ data[i].resize(index[i]);
 for(iter=data[i].begin(); iter!=data[i].end(); iter++) (*iter) = streamData[streamIdx++];
 }

 free(index);    free(streamData);
 return numLists;
 }

 //
 //  Example usage: T1 = int, T2 = list<int>, T3 = list<int>::iterator
 //
 /* SAME THING AS ABOVE, YOU CANNOT DEFINE TEMPLATE FUNCTIONS IN .CPP FILES,
 * THE LINKER WILL TRY TO KILL YOU.
 template<class T1,class T2,class T3> int Save_binListArray(char *filename, vector<T2> &data) {
 unsigned int numEntries1 = data.size(), numEntries2=0;
 unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)*numEntries1);
 for(unsigned int i=0;i<numEntries1;i++) { index[i]=data[i].size(); numEntries2+=index[i]; }
 T1 *streamData = (T1*)malloc(sizeof(T1)*numEntries2);

 int dataArrayIdx=0;   T3 iter;
 for(unsigned int i=0;i<numEntries1;i++)
 for(iter=data[i].begin(); iter!=data[i].end(); iter++)
 streamData[dataArrayIdx++] = *iter;

 FILE *fp = fopen(filename,"w");
 if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return -1; }

 fwrite(&numEntries1,sizeof(int),1,fp);         // Number of lists in data
 fwrite(index,sizeof(int),numEntries1,fp);      // List sizes
 fwrite(streamData,sizeof(T1),numEntries2,fp);  // List elements

 fclose(fp);
 free(index);   free(streamData);
 return 1;
 }

 //
 //  Example usage: T1 = int, T2 = list<int>, T3 = list<int>::iterator
 //
 template<class T1,class T2,class T3> int Save_binListArray(char *filename, list<T2> &data) {
 unsigned int numEntries1 = data.size(), numEntries2=0;
 unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)*numEntries1);
 typename list<T2>::iterator dataIter;
 //	for(unsigned int i=0;i<numEntries1;i++) { index[i]=data[i].size(); numEntries2+=index[i]; }
 unsigned int indexIdx=0;
 for(dataIter = data.begin(); dataIter != data.end(); dataIter++, indexIdx++) { index[indexIdx]=dataIter->size(); numEntries2+=index[indexIdx]; }
 T1 *streamData = (T1*)malloc(sizeof(T1)*numEntries2);

 int dataArrayIdx=0;   T3 iter;
 for(dataIter = data.begin(); dataIter != data.end(); dataIter++)
 for(iter=dataIter->begin(); iter!=dataIter->end(); iter++)
 streamData[dataArrayIdx++] = *iter;

 FILE *fp = fopen(filename,"w");
 if (fp==0) { cerr << "ERROR: cannot open " << filename << "\n";   return -1; }

 fwrite(&numEntries1,sizeof(int),1,fp);         // Number of lists in data
 fwrite(index,sizeof(int),numEntries1,fp);      // List sizes
 fwrite(streamData,sizeof(T1),numEntries2,fp);  // List elements

 fclose(fp);
 free(index);   free(streamData);
 return 1;
 } */
