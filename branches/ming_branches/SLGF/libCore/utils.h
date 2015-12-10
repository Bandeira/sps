#ifndef UTILS_H
#define UTILS_H 

#include "twovalues.h"
#include "Logger.h"

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <set>
#include <map>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <errno.h>
#include <time.h>

using namespace std;

//delimiter for .csv files
extern const char* CSV_SEP;

vector<string> directoryContents(const string & dir,
                                 const string & prefix,
                                 const string &sufix,
                                 bool sort);

/**
 * Compares two strings by size, then alphabetically.
 *@Param: first string
 *@Param: second string
 *@return: Boolean result
 */
bool compare_nocase(string first, string second);


string intToString(int i);
/**
 * string right trim. Aplies to string objects.
 */
#define rtrim(_STR) _STR.erase(_STR.find_last_not_of(string(" \t\f\v\n\r"))+1)
#define ltrim(_STR) _STR.erase(0, _STR.find_first_not_of(string(" \t\f\v\n\r"))+1)


/**
 * Adds an environment variable at the beggining of a command. Used on a system call when in need to previously define the environment variable. The original value is not destroyd; the needed value is added to the existing one.
 *@param: Variable containing the command
 *@param: Variable to be added
 *@param: Variable value
 *@retur : None
 */

void addEnvironmentVariable(string &command,
                            const char *variable,
                            const string &variableValue);

/**
 * Writes a file index into a file
 *@param: file name
 *@param: strings list (file index)
 *@retur : return status: true = ok ; false = error
 */
bool writeFileIndex(const string &fileName, vector<string> &index);
/**
 * Reads a text file content and returns it's contents in lines
 *@param: file name
 *@param: returned set of sub strings
 *@return: return status: true = ok ; false = error
 */
bool readFilesFromFile(const string& str,
                       vector<string>& subStrings,
                       bool fromCurrentDir = false);
/**
 * Reads a file into a buffer
 *@param: file name
 *@param: file size
 *@return: File buffer; null if error
 */
char *readFile(const string &iFileName, unsigned &length);
/**
 * Splits a string into a vector of strings, not discarding empty string
 *@param: string to split
 *@param: returned set of sub strings
 *@param: delimiters; defaultings to space
 */
void stringSplit2(string str, vector<string> &results, const string &delim =
    " ");

/**
 * Splits a string into a vector of strings, discarding empty strings
 *@param: string to split
 *@param: returned set of sub strings
 *@param: delimiters; defaultings to space
 */
void stringSplit(const string& str,
                 vector<string>& subStrings,
                 const string& delimiters = " ");

/**
 * Gets the filename from a full path
 *@param outputFilename - Input path string
 *@param includeExt - include the file extension
 *@return: string containing file name
 */
void extractFileName(const string & path,
                     string &outputFilename,
                     bool includeExt = true);

/**
 * Copies a file to another location
 * @param src Source filename
 * @param dest destination filename
 * @return Returns true, if the file is copied successfully.  Returns
 * false if the src file does not exists
 */
bool copyFile(string src, string dest);

/**
 * Concatenates 2 files and writes them to a new file
 * @param f1 Filename 1
 * @param f2 Filename 2
 * @param dest The file to create from the concatenation of f1 and f2
 * @return Returns true if concatenation succeeds, false otherwise
 */
bool concatenateFiles(string f1, string f2, string dest);

/**
 * Gets file extension from full path
 * @param filename - Input filename
 * @param fileExt - output file extension
 *
 */
void extractFileExt(const string & filename, string & fileExt);

/**
 * Replace occurences in a string
 *@param string context (the string to search. This string WILL be changed)
 *@param string from (what to change)
 *@param string to (replace pattern)
 *@return string: a copy of context string
 */
string& replaceAll(string& context, const string& from, const string& to);

string makeBracketsMods(string& peptide);

/**
 * Makes a directory/folder if it can be made (path to folder is valid and it does  not exist)
 *@param path folder path to create
 *@return true is folder was created, false if not
 */
bool mkdir_if_not_exist(const char* path);

/**
 * Recursively generate directory/folder and all subdirectories if they can be made.
 *@param path folder path to create
 *@return true is folder was created or already exists, false if not
 */
bool mkdirRecurse(const string &path);

/**
 * Maps command-line arguments to their appropriate flags
 *@param argv array of command-line arguments
 *@param argc number of arguments in argv
 *@param flags argument flags to look for
 *@param parsedArgs where to map found flags to their value
 *@return
 */
void parseArguments(char *argv[], int argc, set<string>& flags, map<string,
    string>& parsedArgs);

void
getHistogramInfo(vector<float>& data, vector<TwoValues<float> >& binFreq);

float getResolution(float peakTol);

/**
 * Parses a c string to a float
 *@param str char*
 *@return float value
 */
float getFloat(const char* str);

/**
 * Rounds a decimal number for a # of decimal places
 *@param double value to round
 *@param int number of decimal places
 *@return double rounded value
 */
double roundDouble(double x, int places);

/**
 * Rounds and concatenates a float to an int
 *@param float x
 *@return int value
 */
int floatToInt(float x);
/**
 * Parses a c string to an int
 *@param str char*
 *@return int value
 */
int getInt(const char* str);

/**
 * Depreciated. See MZRange::EqualWithinRange
 */
bool isEqual(float f1, float f2, float range);

/**
 * Parses an int to a string
 *@param x integer
 *@param equalizeLength pad zeros to beginning of string
 *  if its length is less than this
 *@return string
 */
string parseInt(int x, int equalizeLength = 0);

/**
 * Parses a float, which is rounded, to a string
 *@param x float
 *@param prec number of decimal points to round to
 *@return string
 */
string parseFloat(float x, int prec);

/**
 * Splits text by delim to a vector of string
 *@param text char*
 *@param vec output string vector
 *@param delim char*
 */
bool splitText(const char* text, list<string>& vec, const char* delim);

/**
 * Splits string by delim to a vector of strings
 *@param text string
 *@param vec output string vector
 *@param delim string
 */
bool splitText(const char* text, vector<string>& vec, const char* delim);

/**
 * Converts text containing numbers delimeted by delim to a vector of floats
 *@param text char*
 *@param vec output string list
 *@param delim char*
 */
//bool readTextToVector(const char* text, list<float>& vec, const char* delim);

/**
 * From a string of positive integers and integer ranges (ex. "-3,5,7-9"), retrieves all
 *   numbers specified (ex. a set of 0,1,2,3,5,7,8,9).
 *@param text input string of integers and/or integer ranges separated by commas
 *@param nums output set of all numbers specifed by text
 *@return true if parsing suceeded, false if not
 */
bool getRanges(const char* text, set<short>& nums);

bool getdir(string dir, list<string> &files);

class ProgressDisplay
{
public:
  long end;
  int len;
  ProgressDisplay(ostream& output, long e, const char* msg = "Progress: % ");
  void showProgress(ostream& output, long pos);
  void clear(ostream& output);
};

/**
 * TODO: add description
 */
class Utils
{
public:

  /**
   * TODO: add description
   */
  static const double FLOAT_ERR = 0.0001;

  /**
   * TODO: add description
   */
  static vector<float> &unique(vector<float> &v, float resolution, vector<
      unsigned int> *idx = 0);

  /**
   * Intersects two sorted sets of masses.
   *
   *@param v1
   *@param v2
   *@param putHere
   *@param tolerance
   *@return
   */
  static int intersect(vector<float> &v1,
                       vector<float> &v2,
                       vector<float> &putHere,
                       float tolerance);

  /**
   * Saves a set of pairs of floats to binary file.
   * Format: 1 int with number of pairs followed by numEntries pairs of floats.
   *
   *@param filename
   *@param data
   *@return
   */
  static int save_tcb(char *filename, vector<TwoValues<float> > data);

  /**
   * TODO: add description
   *
   *@param filename
   *@param data
   *@return
   */
  static int save_tcb_list(char *filename, list<TwoValues<float> > data);

  /**
   * Converts a set of values to an indicator boolean sps::vector with true
   * whenever a value is present in values.
   *
   *@param values
   *@param output
   *@param tolerance
   *@param multFactor
   *@return
   */
  static unsigned int list_to_vector_bool(vector<float> &values,
                                          vector<bool> &output,
                                          float tolerance,
                                          float multFactor);

  /**
   * Computes a table of binomial probabilities for n=1..maxN, k=1..n;
   * probs[i][j]=p(k=j|n=j+1).
   *
   *@param maxN
   *@param p
   *@param probs
   */
  static void binomial(unsigned short maxN,
                       float p,
                       vector<vector<float> > &probs);

  /**
   * Computes the Gaussian P(X<=x).
   *
   *@param x
   *@param mean
   *@param stddev
   *@return
   */
  static double gaussiancdf(double x, double mean = 0.0, double stddev = 1.0);

  /**
   * Computes scores[i][j]=log(probsSignal[i][j]/probsNoise[i][j]) for all
   * probsNoise[i][j]>0.
   *
   *@param probsSignal
   *@param probsNoise
   *@param scores
   */
  static void logscores(vector<vector<float> > &probsSignal, vector<vector<
      float> > &probsNoise, vector<vector<float> > &scores);

};

/**
 * Filter values from the input argument values - keep only those within
 * tolerance of some entry in validValues. NOTE: Assumes that both values and
 * validValues are sorted in increasing order.
 *
 *@param values
 *@param validValues
 *@param tolerance
 *@return
 */
template<class T> void filterValues(vector<T> &values,
                                    vector<T> &validValues,
                                    T &tolerance,
                                    vector<T> &finalValues)
{
  finalValues.resize(values.size());
  unsigned idxVal = 0, idxValid, idxF = 0;
  while (idxVal < values.size() and idxValid < validValues.size())
  {
    if (fabs(values[idxVal] - validValues[idxValid]) <= tolerance + 0.0001)
    {
      finalValues[idxF++] = values[idxVal++];
      continue;
    }
    if (values[idxVal] < validValues[idxValid])
      idxVal++;
    else
      idxValid++;
  }
  finalValues.resize(idxF);
}

/**
 * TODO: add description
 *
 *@param meansIndexFN
 *@param varsIndexFN
 *@param params
 *@return
 */
bool LoadGaussianParams(const char *meansIndexFN,
                        const char *varsIndexFN,
                        vector<TwoValues<float> > &params);

bool loadCsvToTxt(const char* infilecsv,
                  const char* outfiletxt,
                  float threshold = 0.1);
/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename,
                                    vector<vector<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols;

  if (numLines == 0)
    return -1;
  else
    numCols = data[0].size();

  T *streamData = (T *)malloc(sizeof(T) * numLines * numCols);
  unsigned int i, j;

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for (i = 0; i < numLines; i++)
  {
    if (data[i].size() != numCols)
    {
      cerr << "ERROR in Save_binArray: Not enough data at position " << i
          << "!\n";
      exit(-1);
    }
    for (j = 0; j < numCols; j++)
      streamData[i * numCols + j] = data[i][j];
  }

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename, vector<T> &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols = 1;
  T *streamData = (T *)malloc(sizeof(T) * numLines);
  unsigned int i;

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for (i = 0; i < numLines; i++)
    streamData[i] = data[i];

  fwrite(streamData, sizeof(T), numLines, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename,
                                    vector<TwoValues<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols = 2;
  T *streamData = (T *)malloc(sizeof(T) * numCols * numLines);
  unsigned int i;
  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for (i = 0; i < numLines; i++)
  {
    streamData[i << 1] = data[i][0];
    streamData[(i << 1) + 1] = data[i][1];
  }

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename,
                                    list<TwoValues<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols = 2;
  T *streamData = (T *)malloc(sizeof(T) * numCols * numLines);

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  typename list<TwoValues<T> >::iterator iter;
  unsigned int i = 0;
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    streamData[i++] = (*iter)[0];
    streamData[i++] = (*iter)[1];
  }

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * Writes two binary arrays as if they were 1.  So if matrix 1 is aXb and matrix 2 is a'xb, then the
 * final array will be (a+a')xb
 */
template<class T> int Save_doubleBinArray(const char * filename, vector<vector<T> > & data1, vector<vector<T> > & data2)
{
  FILE * fp;  

  unsigned int numLines = data1.size() + data2.size();
  if (numLines == 0)
    return -1;
  
  typename vector<vector<T> >::iterator iter1 = data1.begin();
  typename vector<vector<T> >::iterator iter2 = data2.begin();

  //num cols must be the same!!
  if(iter1->size() != iter2->size())
    {
      return -1;
    }
  
  unsigned int numCols = iter1->size();
  T *streamData = (T *)malloc(sizeof(T)*numLines*numCols);
  unsigned int i,j;

  //cout << "We have 2 arrays (" << data1.size() << "x" << numCols << ") and (" << data2.size() << "x" << numCols << ")" << endl;

  fp = fopen(filename,"w");

  fwrite(&numLines,sizeof(int),1,fp);
  fwrite(&numCols,sizeof(int),1,fp);
  for(i=0; iter1!=data1.end(); iter1++,i++)
    for(j=0;j<numCols;j++) 
      streamData[i*numCols+j]=(*iter1)[j];
  //cout << "i = " << i << endl;
  for(; iter2!=data2.end(); iter2++,i++)
    for(j=0;j<numCols;j++) 
      streamData[i*numCols+j]=(*iter2)[j];

  fwrite(streamData,sizeof(T),numLines*numCols,fp);
  free(streamData);
  
  fclose(fp);
  return 1;	
}


/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Save_binArray(const char *filename,
                                    list<vector<T> > &data)
{
  FILE *fp;
  unsigned int numLines = data.size(), numCols;
  typename list<vector<T> >::iterator iter = data.begin();
  if (numLines == 0)
    return -1;
  else
    numCols = iter->size();
  T *streamData = (T *)malloc(sizeof(T) * numLines * numCols);
  unsigned int i, j;

  fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numLines, sizeof(int), 1, fp); // Number of lines in data
  fwrite(&numCols, sizeof(int), 1, fp); // Number of cols in data

  for (i = 0; iter != data.end(); iter++, i++)
    for (j = 0; j < numCols; j++)
      streamData[i * numCols + j] = (*iter)[j];

  fwrite(streamData, sizeof(T), numLines * numCols, fp);
  free(streamData);

  fclose(fp);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */

template<class T, template < typename U, typename = std::allocator<U> > class C> int
Load_binArray(const char *filename, C< C<T> > &data)
{
  FILE *fp; int numLines=0, numCols=0;
  unsigned int i,j; int n;

  fp = fopen(filename,"r");
  if (fp==0)
  { std::cerr << "ERROR: cannot open " << filename << "\n"; return -1;}

  n=fread(&numLines,sizeof(int),1,fp); if(numLines<0 | n!=1)
  { std::cerr << "Invalid number of lines ("<<numLines<<") in "<<filename<<"\n";return -1;}
  n=fread(&numCols,sizeof(int),1,fp); if(numCols<0 | n!=1)
  { std::cerr << "Invalid number of columns ("<<numCols<<") in "<<filename<<"\n";return -1;}
  if(numLines==0 or numCols==0)
  { data.resize(0); return 0;}

  T *streamData = (T *)malloc(sizeof(T)*numLines*numCols);
  if(!streamData)
  { cerr<<"(Load_binArray) Error allocating memory to read "<<numLines<<"/"<<numCols<<" lines/cols of "<<sizeof(T)<<"-byte elements!\n"; return -1;}
  n=fread(streamData,sizeof(T),numLines*numCols,fp);
  fclose(fp);
  if(n!=(numLines*numCols))
  { std::cerr << "Error reading binary array: not enough elements in the file.\n"; free(streamData); return -1;}

  data.resize(numLines);
  if(data.size()!=numLines)
  { cerr<<"(Load_binArray) Error allocating memory to store "<<numLines<<"/"<<numCols<<" lines/cols of "<<sizeof(T)<<"-byte elements!\n"; free(streamData); return(-1);}
  unsigned int streamIdx=0;
  for(i=0;i<(unsigned int)numLines;i++)
  { data[i].resize(numCols);
    for(j=0;j<(unsigned int)numCols;j++) data[i][j] = streamData[streamIdx++];
  }

  free(streamData);
  return 1;
}



/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T, template < typename U, typename = std::allocator<U> > class C> int
Load_binArray(const char *filename, C<T> &data)
{
  FILE *fp; int numLines=0, numCols=0;
  unsigned int i,j; int n;

  fp = fopen(filename,"r");
  if (fp==0)
  { std::cerr << "ERROR: cannot open " << filename << "\n"; return -1;}

  n=fread(&numLines,sizeof(int),1,fp); if(numLines<0 | n!=1)
  { std::cerr << "Invalid number of lines ("<<numLines<<") in "<<filename<<"\n";return -1;}
  n=fread(&numCols,sizeof(int),1,fp); if(numCols!=1 | n!=1)
  { cerr << "Invalid number of columns ("<<numCols<<") in "<<filename<<"\n";return -1;}
  if(numLines==0)
  { data.resize(0); return 0;}

  T *streamData = (T *)malloc(sizeof(T)*numLines);
  if(!streamData)
  { cerr<<"(Load_binArray) Error allocating memory to read "<<numLines<<" lines of "<<sizeof(T)<<"-byte elements!\n"; return -1;}
  n=fread(streamData,sizeof(T),numLines,fp);
  fclose(fp);
  if(n!=(numLines))
  { std::cerr << "Error reading binary array: not enough elements in the file.\n"; free(streamData); return -1;}

  data.resize(numLines);
  if(data.size()!=numLines)
  { cerr<<"(Load_binArray) Error allocating memory to store "<<numLines<<" lines of "<<sizeof(T)<<"-byte elements!\n"; free(streamData); return -1;}
  for(i=0;i<(unsigned int)numLines;i++) data[i]=streamData[i];

  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> int Load_binArray(const char *filename,
                                    vector<TwoValues<T> > &data)
{

  FILE *fp;
  int numLines = 0, numCols = 0, n = 0;

  fp = fopen(filename, "r");
  if (fp == 0)
  {
    std::cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  n = fread(&numLines, sizeof(int), 1, fp);
  if (numLines < 1)
  {
    std::cerr << "Invalid number of lines (" << numLines << ") in " << filename
        << "\n";
    return -1;
  }
  n = fread(&numCols, sizeof(int), 1, fp);
  if (numCols != 2)
  {
    cerr << "Invalid number of columns (" << numCols << ") in " << filename
        << "\n";
    return -1;
  }
  if (numLines == 0)
  {
    data.resize(0);
    return 0;
  }

  data.resize(numLines);
  if (data.size() != numLines)
  {
    cerr << "(Load_binArray) Error allocating memory to store " << numLines
        << " lines of " << sizeof(T) << "-byte elements!\n";
    return -1;
  }

  char * streamData = (char*)malloc(numLines * numCols * sizeof(T));
  if (streamData == 0x0)
  {
    cerr << "(Load_binArray) Error allocating memory to store " << numLines
        << " lines of " << sizeof(T) << "-byte elements!\n";
    return -1;
  }
  n = fread(streamData, sizeof(T), numLines * numCols, fp);
  fclose(fp);

  unsigned int streamIdx = 0;
  for (unsigned int i = 0; i < numLines; i++)
  {
    data[i][0] = *((T*)&streamData[streamIdx]);
    streamIdx += sizeof(T);
    data[i][1] = *((T*)&streamData[streamIdx]);
    streamIdx += sizeof(T);
  }

  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param numLinesCols
 *@return
 */
inline int Load_binArraySize(const char *filename,
                             TwoValues<unsigned int> &numLinesCols)
{
  FILE *fp;
  numLinesCols.set(0, 0);
  fp = fopen(filename, "r");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }
  TwoValues<int> sz(0, 0);

  unsigned int n;
  n = fread(&sz[0], sizeof(int), 1, fp);
  if ((sz[0] < 0) | (n != 1))
  {
    cerr << "Invalid number of lines (" << sz[0] << ") in " << filename << "\n";
    return -1;
  }
  n = fread(&sz[1], sizeof(int), 1, fp);
  if ((sz[1] < 0) | (n != 1))
  {
    cerr << "Invalid number of columns (" << sz[1] << ") in " << filename
        << "\n";
    return -1;
  }
  numLinesCols.set((unsigned int)sz[0], (unsigned int)sz[1]);

  fclose(fp);
  return numLinesCols[0];
}
//template<class T1,class T2,class T3> int Save_binListArray(char *filename, vector<T2> &data);

// T1 = element type, T2 = vector<T1>/list<T1>, T3 = T2::iterator

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T2, class T3> int Load_binListArray(const char *filename,
                                                             vector<T2> &data)
{
  FILE *fp;
  unsigned int numLists = 0, numElems = 0;
  unsigned int *index;
  unsigned int i, streamIdx;
  int n;

  fp = fopen(filename, "r");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  // Read number of lists
  n = fread(&numLists, sizeof(unsigned int), 1, fp);
  if (numLists <= 0 | n != 1)
  {
    cerr << "Invalid number of lists (" << numLists << ")\n";
    return -1;
  }

  // Read index
  index = (unsigned int*)malloc(sizeof(unsigned int) * numLists);
  n = fread(index, sizeof(unsigned int), numLists, fp);
  if (n != numLists)
  {
    cerr << "Invalid index!\n";
    free(index);
    return -1;
  }
  numElems = 0;
  for (i = 0; i < numLists; i++)
    numElems += index[i];

  // Read data
  T1 *streamData = (T1 *)malloc(sizeof(T1) * numElems);
  n = fread(streamData, sizeof(T1), numElems, fp);
  fclose(fp);
  if (n != (int)numElems)
  {
    cerr
        << "Error reading binary list array: not enough elements in the file.\n";
    free(streamData);
    free(index);
    return -1;
  }

  // Convert back to C++ structure
  data.resize(numLists);
  T3 iter;
  streamIdx = 0;
  for (i = 0; i < numLists; i++)
  {
    data[i].resize(index[i]);
    for (iter = data[i].begin(); iter != data[i].end(); iter++)
      (*iter) = streamData[streamIdx++];
  }

  free(index);
  free(streamData);
  return numLists;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T2, class T3> int Load_binListArray(const char *filename,
                                                             list<T2> &data)
{
  FILE *fp;
  unsigned int numLists = 0, numElems = 0;
  unsigned int *index;
  int n;

  fp = fopen(filename, "r");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  // Read number of lists
  n = fread(&numLists, sizeof(unsigned int), 1, fp);
  if (numLists <= 0 | n != 1)
  {
    cerr << "Invalid number of lists (" << numLists << ")\n";
    return -1;
  }

  // Read index
  index = (unsigned int*)malloc(sizeof(unsigned int) * numLists);
  n = fread(index, sizeof(unsigned int), numLists, fp);
  if (n != numLists)
  {
    cerr << "Invalid index!\n";
    free(index);
    return -1;
  }
  numElems = 0;
  for (unsigned int i = 0; i < numLists; i++)
    numElems += index[i];

  // Read data
  T1 *streamData = (T1 *)malloc(sizeof(T1) * numElems);
  n = fread(streamData, sizeof(T1), numElems, fp);
  fclose(fp);
  if (n != (int)numElems)
  {
    cerr
        << "Error reading binary list array: not enough elements in the file.\n";
    free(streamData);
    free(index);
    return -1;
  }

  // Convert back to C++ structure
  data.resize(numLists);
  T3 iter;
  unsigned int i = 0;
  for (iter = data.begin(); iter != data.end(); iter++)
  {
    iter->resize(index[i]);
    for (unsigned int j = 0; j < index[i]; j++)
    {
      (*iter)[j] = streamData[j];
    }
    i++;
  }

  free(index);
  free(streamData);
  return numLists;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T1b, class T2, class T3> int Load_binListArray(const char *filename,
                                                                        vector<
                                                                            T2> &data)
{
  FILE *fp;
  unsigned int numLists = 0, numElems = 0;
  unsigned int *index;
  unsigned int i, streamIdx;
  int n;

  fp = fopen(filename, "r");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  // Read number of lists
  n = fread(&numLists, sizeof(unsigned int), 1, fp);
  if (numLists <= 0 | n != 1)
  {
    cerr << "Invalid number of lists (" << numLists << ")\n";
    return -1;
  }

  // Read index
  index = (unsigned int*)malloc(sizeof(unsigned int) * numLists);
  n = fread(index, sizeof(unsigned int), numLists, fp);
  if (n != numLists)
  {
    cerr << "Invalid index!\n";
    free(index);
    return -1;
  }
  numElems = 0;
  for (i = 0; i < numLists; i++)
  {
    numElems += index[i];
  }

  // Read data
  char *streamData = (char *)malloc((sizeof(T1) + sizeof(T1b)) * numElems);
  n = fread(streamData, sizeof(T1) + sizeof(T1b), numElems, fp);
  fclose(fp);
  if (n != (int)numElems)
  {
    cerr
        << "Error reading binary list array: not enough elements in the file.\n";
    free(streamData);
    free(index);
    return -1;
  }

  // Convert back to C++ structure
  data.resize(numLists);
  T3 iter;
  streamIdx = 0;
  for (i = 0; i < numLists; i++)
  {
    data[i].resize(index[i]);
    for (iter = data[i].begin(); iter != data[i].end(); iter++)
    {
      iter->first = *((T1*)&streamData[streamIdx]);
      streamIdx += sizeof(T1);
      iter->second = *((T1b*)&streamData[streamIdx]);
      streamIdx += sizeof(T1b);
    }
  }

  free(index);
  free(streamData);
  return numLists;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T2, class T3> int Save_binListArray(const char *filename,
                                                             vector<T2> &data)
{
  unsigned int numEntries1 = data.size(), numEntries2 = 0;
  unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)
      * numEntries1);
  for (unsigned int i = 0; i < numEntries1; i++)
  {
    index[i] = data[i].size();
    numEntries2 += index[i];
  }
  T1 *streamData = (T1*)malloc(sizeof(T1) * numEntries2);

  int dataArrayIdx = 0;
  T3 iter;
  for (unsigned int i = 0; i < numEntries1; i++)
    for (iter = data[i].begin(); iter != data[i].end(); iter++)
      streamData[dataArrayIdx++] = *iter;

  FILE *fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries1, sizeof(int), 1, fp); // Number of lists in data
  fwrite(index, sizeof(int), numEntries1, fp); // List sizes
  fwrite(streamData, sizeof(T1), numEntries2, fp); // List elements

  fclose(fp);
  free(index);
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T2, class T3> int Save_binListArray(const char *filename,
                                                             list<T2> &data)
{
  unsigned int numEntries1 = data.size(), numEntries2 = 0;
  unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)
      * numEntries1);
  typename list<T2>::iterator dataIter;
  //	for(unsigned int i=0;i<numEntries1;i++) { index[i]=data[i].size(); numEntries2+=index[i]; }
  unsigned int indexIdx = 0;
  for (dataIter = data.begin(); dataIter != data.end(); dataIter++, indexIdx++)
  {
    index[indexIdx] = dataIter->size();
    numEntries2 += index[indexIdx];
  }
  T1 *streamData = (T1*)malloc(sizeof(T1) * numEntries2);

  int dataArrayIdx = 0;
  T3 iter;
  for (dataIter = data.begin(); dataIter != data.end(); dataIter++)
    for (iter = dataIter->begin(); iter != dataIter->end(); iter++)
      streamData[dataArrayIdx++] = *iter;

  FILE *fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries1, sizeof(int), 1, fp); // Number of lists in data
  fwrite(index, sizeof(int), numEntries1, fp); // List sizes
  fwrite(streamData, sizeof(T1), numEntries2, fp); // List elements

  fclose(fp);
  free(index);
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T1b, class T2, class T3> int Save_binListArray(const char *filename,
                                                                        vector<
                                                                            list<
                                                                                pair<
                                                                                    T1,
                                                                                    T1b> > > &data)
{
  unsigned int numEntries1 = data.size(), numEntries2 = 0;
  unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)
      * numEntries1);
  for (unsigned int i = 0; i < numEntries1; i++)
  {
    index[i] = data[i].size();
    numEntries2 += index[i];
  }
  char * streamData = (char*)malloc((sizeof(T1) + sizeof(T1b)) * numEntries2);

  int dataArrayIdx = 0;
  T3 iter;
  for (unsigned int i = 0; i < numEntries1; i++)
  {
    for (iter = data[i].begin(); iter != data[i].end(); iter++)
    {
      strncpy(&streamData[dataArrayIdx], (char*)&iter->first, 4);
      dataArrayIdx += sizeof(T1);
      strncpy(&streamData[dataArrayIdx], (char*)&iter->second, 4);
      dataArrayIdx += sizeof(T1b);
    }
  }

  FILE *fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries1, sizeof(int), 1, fp); // Number of lists in data
  fwrite(index, sizeof(int), numEntries1, fp); // List sizes
  fwrite(streamData, sizeof(T1) + sizeof(T1b), numEntries2, fp); // List elements

  fclose(fp);
  free(index);
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T1, class T1b, class T2, class T3> int Save_binListArray(const char *filename,
                                                                        list<
                                                                            vector<
                                                                                pair<
                                                                                    T1,
                                                                                    T1b> > > &data)
{
  unsigned int numEntries1 = data.size(), numEntries2 = 0;
  unsigned int *index = (unsigned int*)malloc(sizeof(unsigned int)
      * numEntries1);
  typename list<T2>::iterator dataIter;
  unsigned int indexIdx = 0;
  for (dataIter = data.begin(); dataIter != data.end(); dataIter++, indexIdx++)
  {
    index[indexIdx] = dataIter->size();
    numEntries2 += index[indexIdx];
  }

  char * streamData = (char*)malloc((sizeof(T1) + sizeof(T1b)) * numEntries2);

  int dataArrayIdx = 0;
  T3 iter;
  for (dataIter = data.begin(); dataIter != data.end(); dataIter++)
  {
    for (iter = dataIter->begin(); iter != dataIter->end(); iter++)
    {
      strncpy(&streamData[dataArrayIdx], (char*)&iter.first, 4);
      dataArrayIdx += sizeof(T1);
      strncpy(&streamData[dataArrayIdx], (char*)&iter.second, 4);
      dataArrayIdx += sizeof(T1b);
    }
  }

  FILE *fp = fopen(filename, "w");
  if (fp == 0)
  {
    cerr << "ERROR: cannot open " << filename << "\n";
    return -1;
  }

  fwrite(&numEntries1, sizeof(int), 1, fp); // Number of lists in data
  fwrite(index, sizeof(int), numEntries1, fp); // List sizes
  fwrite(streamData, sizeof(T1) + sizeof(T1b), numEntries2, fp); // List elements

  fclose(fp);
  free(index);
  free(streamData);
  return 1;
}

/**
 * TODO: add description
 */
class BufferedLineReader
{
  char *lines;
  vector<unsigned int> linesIdx;
public:
  BufferedLineReader()
  {
    lines = NULL;
  }
  ~BufferedLineReader()
  {
    reset();
  }

  /**
   * TODO: add description
   */
  void reset()
  {
    if (lines)
      free(lines);
    lines = NULL;
    linesIdx.resize(0);
  }

  /**
   * TODO: add description
   *
   *@return the number of lines currently in the reader.
   */
  unsigned int size()
  {
    return linesIdx.size();
  }

  /**
   * TODO: add description
   *
   *@param lineNum the line number to get.
   *@return
   */
  char *getline(unsigned int lineNum)
  {
    if (lineNum >= linesIdx.size())
      return NULL;
    else
      return &lines[linesIdx[lineNum]];
  }

  /**
   * TODO: add description
   *
   *@param filename
   *@return
   */
  short Load(const char *filename);
};

/**
 * TODO: add description
 *
 *@param filename
 *@param data
 *@return
 */
template<class T> bool Load_binArray_multiple(const char *filename, vector<
    vector<T> > &data)
{
  BufferedLineReader blrNames;
  unsigned int fIdx, numFiles;
  data.resize(0);
  if (blrNames.Load(filename) <= 0)
  {
    cerr << "ERROR reading " << filename << "!\n";
    return false;
  }
  numFiles = blrNames.size();

  // Count total number of entries
  TwoValues<unsigned int> totLinesCols(0, 0), curLinesCols(0, 0);
  for (fIdx = 0; fIdx < numFiles; fIdx++)
  {
    if (strlen(blrNames.getline(fIdx)) == 0)
      continue;
    const char* fName = blrNames.getline(fIdx);
    if (Load_binArraySize(fName, curLinesCols) < 0)
      return false;
    if (fIdx == 0)
      totLinesCols[1] = curLinesCols[1];
    else if (totLinesCols[1] != curLinesCols[1])
    {
      cerr << "ERROR: Inconsistent number of columns in "
          << blrNames.getline(fIdx) << " (Load_binArray_multiple)\n";
      return false;
    }
    totLinesCols[0] += curLinesCols[0];
  }

  // Load data
  FILE *fp;
  vector < vector<T> > curData;
  data.resize(totLinesCols[0]);
  unsigned int pivot, idxData = 0;
  for (fIdx = 0; fIdx < numFiles; fIdx++)
  {
    if (strlen(blrNames.getline(fIdx)) == 0)
      continue;
    const char* fName = blrNames.getline(fIdx);
    if (Load_binArray(fName, curData) < 0)
      return false;
    if (idxData + curData.size() > data.size())
    {
      cerr << "ERROR: Too many lines in " << blrNames.getline(fIdx)
          << " (Load_binArray_multiple)\n";
      return false;
    }
    for (pivot = 0; pivot < curData.size(); pivot++)
      data[idxData++] = curData[pivot];
  }

  return true;
}

std::string get_time_string();

#endif
