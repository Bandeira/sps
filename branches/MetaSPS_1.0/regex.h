#ifndef REGEX_H
#define REGEX_H

#include <regex.h>


namespace sps
{

/**
    C++ regular expression wrapper.

    @warning  Does not handle erroneous patterns.
*/

class regex
{
  regex_t cr;

public:
  regex(const char * spattern)
  {
    regcomp(& cr, spattern, REG_EXTENDED | REG_NOSUB);
  }

  ~regex()
  {
    regfree(& cr);
  }

  friend bool operator == (const regex & cr, const char * sstring)
  {
    return regexec(& cr.cr, sstring, (size_t) 0, 0, 0) == 0;
  }

  friend bool operator != (const regex & cr, const char * sstring)
  {
    return ! operator == (cr, sstring);
  }
};

} // namespace sps

#endif
