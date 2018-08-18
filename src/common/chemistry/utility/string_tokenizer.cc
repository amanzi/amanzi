/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  This class is a simple string tokenizer.
*/

#include "string_tokenizer.hh"

#include <sstream>
#include <string>
#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

void StringTokenizer::tokenize(const std::string& source,
                               const std::string& delimiters) {
  clear();
  std::string::size_type spos(source.find_first_not_of(delimiters, 0));
  std::string::size_type epos(source.find_first_of(delimiters, spos));

  while (std::string::npos != epos || std::string::npos != spos) {
    push_back(source.substr(spos, epos - spos));
    spos = source.find_first_not_of(delimiters, epos);
    epos = source.find_first_of(delimiters, spos);
  }
}


void StringTokenizer::tokenize_leave_delimiters(const std::string& source,
                                               const std::string& delimiters) {
  clear();
  std::string::size_type spos(source.find_first_not_of(delimiters, 0));
  std::string::size_type epos(source.find_first_of(delimiters, spos));

  // add delimeter
  if (spos > 0) push_back(source.substr(0, spos));
  while (std::string::npos != epos || std::string::npos != spos) {
    push_back(source.substr(spos, epos - spos));
    // find position of delimiter
    spos = epos;
    epos = source.find_first_not_of(delimiters, spos);
    if (std::string::npos != epos) {
      // add delimeter
      push_back(source.substr(spos, epos - spos));
      // find position of non-delimiter
      spos = source.find_first_not_of(delimiters, epos);
      epos = source.find_first_of(delimiters, spos);
    }
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
