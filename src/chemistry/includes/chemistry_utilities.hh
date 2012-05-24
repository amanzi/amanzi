/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_UTILITIES_HH_
#define AMANZI_CHEMISTRY_UTILITIES_HH_

#include <vector>
#include <sstream>
#include <string>

//
// Common stand alone utility functions
//

#include "chemistry_output.hh"

namespace amanzi {
namespace chemistry {

extern ChemistryOutput* chem_out;

namespace utilities {

/*******************************************************************************
 **
 **  Custom comparison operators
 **
 ******************************************************************************/
bool CaseInsensitiveStringCompare(const std::string& string1, 
                                  const std::string& string2);
bool CompareFabs(const double& a, const double& b);

/*******************************************************************************
 **
 **  String conversion utilities
 **
 ******************************************************************************/
void LowerCaseString(const std::string& in, std::string* out);
void RemoveLeadingAndTrailingWhitespace(std::string* line);

/*******************************************************************************
 **
 **  Math conversion utilities
 **
 ******************************************************************************/
double log_to_ln(double d);
double ln_to_log(double d);

/*******************************************************************************
 **
 **  Print Utilities
 **
 ******************************************************************************/
template <typename T>
void PrintVector(const std::string& name, 
                 const std::vector<T>& data) {
  std::stringstream output;
  output << name << " : { ";
  for (typename std::vector<T>::const_iterator i = data.begin();
       i != data.end(); ++i) {
    output << *i << "  ";
  }
  output << " }\n";
  chem_out->Write(kVerbose, output);
}  // end PrintVector


}  // namespace utilities
}  // namespace chemistry
}  // namespace amanzi
#endif  // AMANZI_CHEMISTRY_UTILITIES_HH_
