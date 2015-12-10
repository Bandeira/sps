#include "Tokenizer.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////
namespace spsReports {
////////////////////////////////////////////////////////////////////////////////

void StringExplode(string str, vector<string> &results, bool removeDelimiters, string separators, string delimiters)
{
  size_t found;

  string all = separators + delimiters;

  // find first of something
  found = str.find_first_of(all);

  while(found != string::npos) {

    // check if a delimiter was found.
    if(delimiters.find_first_of(str[found]) != string::npos) {
      // If so, find matching delimiter
      size_t found2 = str.find_first_of(str[found], found+1);
      // find next, if not yet at the end of the string
      if(found2 == string::npos) {
        // if at the end of the string, remove delimiter and break cycle
        if(str.length() > found)
          str.erase(found,1);
        break;
      }

      // auxiliary used for delimiter removal
      size_t found3 = str.find_first_of(all, found2+1);

      if(removeDelimiters) {
        int dec = 0;
        if(str.length() > found2) {
          str.erase(found2,1);
          dec++;
        }

        if(str.length() > found) {
          str.erase(found,1);
          dec++;
        }
        found3 -= dec;
      }

      found = found3;

    } else {

      if(found > 0) {
        results.push_back(str.substr(0,found));
      }

      if(str.length() > found)
        str = str.substr(found+1);
      else
        str = "";

      found = str.find_first_of(all);
    }
  }
  if(str.length() > 0){
    results.push_back(str);
  }
}

///////////////////////////////////////////////////////////////////////////////
bool stringSortCompare(const string &s1, const string &s2)
{
  if(s1.length() < s2.length())
    return true;
  else if(s1.length() > s2.length())
    return false;
  else
    return (s1.compare(s2) < 0);
}
///////////////////////////////////////////////////////////////////////////////
bool stringUniqueCompare(const string &s1, const string &s2)
{
  return (s1.compare(s2) == 0);
}
///////////////////////////////////////////////////////////////////////////////

};
///////////////////////////////////////////////////////////////////////////////
