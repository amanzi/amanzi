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


namespace Amanzi {

// remove leading and trailing spaces.
inline 
void trim(std::string& str)
{
  size_t first = str.find_first_not_of("\n ");
  if (first != std::string::npos) {
    size_t last = str.find_last_not_of("\n ");
    str = str.substr(first, (last - first + 1));
  }
}


// out
inline 
void replace_all_copies(std::string& str,
                        const std::string& sub_old, const std::string& sub_new)
{
  size_t pos(0);
  int n = sub_old.size();

  while (true) {
    pos = str.find(sub_old, pos);
    if (pos == std::string::npos) break;

    str.replace(pos, n, sub_new);
    pos += n;
  }
}

}  // namespace Amanzi

