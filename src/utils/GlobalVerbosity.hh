/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_GLOBAL_VERBOSITY_
#define AMANZI_GLOBAL_VERBOSITY_

#include "Teuchos_VerboseObject.hpp"

namespace Amanzi {

class VerbosityLevel {
 public:
  static Teuchos::EVerbosityLevel level_;

  static const unsigned int kVerbositySize = 16;
  static std::string verbosityHeader(const std::string str)
  {
    std::string header(str);
    if (header.size() > kVerbositySize) {
      header.erase(kVerbositySize);
    } else if (header.size() < kVerbositySize) {
      header.append(kVerbositySize - header.size(), ' ');
    }
    return header;
  };
};

} // namespace Amanzi

#endif
