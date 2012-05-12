/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "chemistry_utilities.hh"

#include <cmath>

#include <string>
#include <sstream>
#include <vector>

#include "chemistry_output.hh"

namespace amanzi {
namespace chemistry {
namespace utilities {

/*******************************************************************************
 **
 **  Custom comparison operators
 **
 ******************************************************************************/
bool CaseInsensitiveStringCompare(const std::string& string1, 
                                  const std::string& string2) {
  // really CaseInsensitiveStringsEqual...

  // if not the same length, not the same. don't bother checking
  // individual characters
  if (string1.size() != string2.size()) {
    return false;
  }
  // loop through each string, check each character individually.
  std::string::const_iterator char1, char2;
  for (char1 = string1.begin(), char2 = string2.begin(); 
       char1 != string1.end(); ++char1, ++char2) {
    if (std::tolower(*char1) != std::tolower(*char2)) {
      return false;
    }
  }
  return true;
}  // end CaseInsensitiveStringCompare()

      // std::stringstream output;
      // output << "CICS : strings are not the same length\n"
      //        << "  string one : " << string1
      //        << "\n  string two : " << string2 << std::endl;
      // chem_out.Write(kVerbose, output.str());

      // std::stringstream output;
      // output << "CICS : strings differ at position "
      //        << std::distance(string1.begin(), char1)
      //        << "\n  string one value: " << *char1
      //        << "\n  string two value: " << *char2 << std::endl;
      // chem_out.Write(kVerbose, output.str());


bool CompareFabs(const double& a, const double& b) {
  /* for use with stl algorithms such as max element. Must "return
   * true if the first argument is to be considered less than the
   * second argument, and false otherwise" */
  return std::fabs(a) < std::fabs(b);
}  // end CompareFabs()


/*******************************************************************************
 **
 **  Conversion routines
 **
 ******************************************************************************/

/*
** convert string to lower case
*/
void LowerCaseString(const std::string& in, std::string* out) {
  *out = in;
  for (std::string::iterator c = out->begin(); c != out->end(); ++c) {
    *c = std::tolower(*c);
  }
}

/*
** ln(x) = log10(x) * ln(10)
**
** input: log10(x)
** output: ln(x)
*/ 
double log_to_ln(double log10_x) {
  return log10_x * 2.30258509299;
}

/*
** log10(x) = ln(x) * ln(10)
**
** input: ln(x)
** return: log10(x)
*/
double ln_to_log(double ln_x) {
  return ln_x * 0.434294481904;
}

}  // namespace utilities
}  // namespace chemistry
}  // namespace amanzi
