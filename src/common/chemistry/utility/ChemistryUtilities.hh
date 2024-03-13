/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson <jnjohnson@lbl.gov>
*/

/*
  Chemistry

  Common stand alone utility functions
*/

#ifndef AMANZI_CHEMISTRY_UTILITIES_HH_
#define AMANZI_CHEMISTRY_UTILITIES_HH_

#include <vector>
#include <sstream>
#include <string>
#include <iomanip>

#include "VerboseObject.hh"

namespace Amanzi {
namespace AmanziChemistry {

namespace utilities {

inline double
log_to_ln(double d)
{
  return d * 2.30258509299405;
}
inline double
ln_to_log(double d)
{
  return d * 0.434294481903252;
}

template <typename T>
void
PrintVector(const std::string& name,
            const std::vector<T>& data,
            const Teuchos::RCP<VerboseObject>& vo,
            const int precision = -1,
            const bool comma_seperated = false)
{
  std::stringstream output;
  if (precision > 0) { output << std::setprecision(precision); }
  output << name << " : { ";
  for (auto i = data.begin(); i != data.end(); ++i) {
    output << *i;
    if (i != --data.end()) {
      if (comma_seperated) {
        output << ", ";
      } else {
        output << "  ";
      }
    }
  }
  output << " }\n";
  vo->Write(Teuchos::VERB_HIGH, output.str());
}

} // namespace utilities
} // namespace AmanziChemistry
} // namespace Amanzi
#endif
