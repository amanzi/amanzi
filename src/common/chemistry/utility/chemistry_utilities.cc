/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#include "chemistry_utilities.hh"

#include <cmath>

#include <string>
#include <sstream>
#include <vector>

namespace Amanzi {
namespace AmanziChemistry {
namespace utilities {

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

void RemoveLeadingAndTrailingWhitespace(std::string* line) {
  std::string whitespace(" \t\f\v\n\r");
  size_t start = line->find_first_not_of(whitespace);
  if (start != std::string::npos) {
    line->erase(0, start);
  } else if (start == std::string::npos) {
    // entire line is blank
    line->erase(0);
  }
  size_t end = line->find_last_not_of(whitespace);
  if (end != std::string::npos) {
    ++end;  // find returned the last non-whitespace character....
    line->erase(end);
  }
}  // end RemoveLeadingAndTrailingWhitespace()


/*
** ln(x) = log10(x) * ln(10)
**
** input: log10(x)
** output: ln(x)
*/ 
double log_to_ln(double log10_x) {
  return log10_x * 2.30258509299405;
}

/*
** log10(x) = ln(x) * ln(10)
**
** input: ln(x)
** return: log10(x)
*/
double ln_to_log(double ln_x) {
  return ln_x * 0.434294481903252;
}

}  // namespace utilities
}  // namespace AmanziChemistry
}  // namespace Amanzi
