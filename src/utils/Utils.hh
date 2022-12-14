/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _UTILS_HH_
#define _UTILS_HH_

#include <vector>
#include <ostream>

namespace Amanzi {
namespace Utils {

template <typename M>
std::ostream&
dump_map(const M& map, std::ostream& stream)
{
  for (typename M::const_iterator it = map.begin(); it != map.end(); ++it) {
    stream << it->first << ": " << it->second << std::endl;
  }

  return stream;
}

template <typename V>
std::ostream&
dump_vector_as_map(const V& data, std::ostream& stream)
{
  int local = 0;
  for (typename V::const_iterator it = data.begin(); it != data.end(); ++it, ++local) {
    stream << local << ": " << *it << std::endl;
  }
  return stream;
}

} // namespace Utils
} // namespace Amanzi
#endif
