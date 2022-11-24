/*
   Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
   Amanzi is released under the three-clause BSD License.
   The terms of use and "as is" disclaimer for this license are
   provided in the top-level COPYRIGHT file.

   Author: Konstantin Lipnikov

   Useful non-member functions.
*/

#include <iostream>
#include <string>
#include <vector>

namespace Amanzi {

// remove leading and trailing spaces.
inline void
trim(std::string& str)
{
  size_t first = str.find_first_not_of("\n ");
  if (first != std::string::npos) {
    size_t last = str.find_last_not_of("\n ");
    str = str.substr(first, (last - first + 1));
  }
}


// out
inline void
replace_all(std::string& str, const std::string& sub_old, const std::string& sub_new)
{
  size_t pos(0);
  int n0 = sub_old.size();
  int n1 = sub_new.size();

  while (true) {
    pos = str.find(sub_old, pos);
    if (pos == std::string::npos) break;

    if (n0 == n1) {
      str.replace(pos, n0, sub_new);
    } else {
      str.insert(pos, sub_new, 0, n1);
      pos += n1;
      str.erase(pos, n0);
    }

    pos += n0;
  }
}


// split string using multiple delimeters and ignoring empty strings
inline std::vector<std::string>
split(const std::string& str, const std::string& delimiters)
{
  std::vector<std::string> result;
  size_t current, next(-1);
  do {
    next = str.find_first_not_of(delimiters, next + 1);
    if (next == std::string::npos) break;

    next -= 1;
    current = next + 1;
    next = str.find_first_of(delimiters, current);
    result.push_back(str.substr(current, next - current));
  } while (next != std::string::npos);

  return result;
}

} // namespace Amanzi
