/**
 * Helper functions related to Inspect result parsing
 *
 *  Created on: Aug 18, 2010
 *      Author: jsnedecor
 */

#ifndef INSPECT_PARSE_H_
#define INSPECT_PARSE_H_

#include "utils.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <cstdio>
#include <cstring>

class InspectResultsLine
{
public:
  bool parseFromFields(const vector<string>& fields,
                       const vector<string>& field_names);

  string SpectrumFile;
  int scan;
  string Annotation;
  string Protein;
  int Charge;
  float MQScore;
  float Score;
  int Length;
  float TotalPRMScore;
  float MedianPRMScore;
  float FractionY;
  float FractionB;
  float Intensity;
  int NTT;
  float p_value;
  float F_Score;
  float DeltaScore;
  float DeltaScoreOther;
  int RecordNumber;
  int DBFilePos;
  int SpecFilePos;
};

/**Inspect parsing and translation between SpecNets and Inspect formats
 */
class InspectResultsSet
{
public:
  map<int, InspectResultsLine> results; //scan number, InspectResultsLine

  void inspectToSpecNets(string &inspect, string &specnets);

  bool loadInspectResultsFile(const char *inspect_results_file);

  InspectResultsLine* getScan(int scan);
};

#endif /* INSPECT_PARSE_H_ */
