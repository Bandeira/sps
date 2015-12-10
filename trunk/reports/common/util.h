#ifndef __UTIL_H__
#define __UTIL_H__


#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>



class BadConversion : public std::runtime_error {
public:
  BadConversion(std::string const& s)
   : std::runtime_error(s)
   { }
};

inline double convertToDouble(std::string const& s)
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
   throw BadConversion("convertToDouble(\"" + s + "\")");
  return x;
}

#endif
