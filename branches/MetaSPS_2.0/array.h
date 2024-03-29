#ifndef ARRAY_H
#define ARRAY_H

#include <algorithm>


namespace sps
{

template <typename T, std::size_t S>
  struct array
  {
    typedef T                   value_type;
    typedef value_type &        reference;
    typedef const value_type &  const_reference;
    typedef std::size_t         size_type;
    typedef std::ptrdiff_t      difference_type;

    value_type data[S];

    array() : data() {}

    reference operator [] (size_type ns)
    {
      return data[ns];
    }

    const_reference operator [] (size_type ns) const
    {
      return data[ns];
    }

    size_type size() const
    {
      return S;
    }
  };

} // namespace sps


#endif
