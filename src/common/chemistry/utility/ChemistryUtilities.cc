/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson <jnjohnson@lbl.gov>

  Common stand alone utility functions
*/

#include <cmath>
#include <string>
#include <vector>

#include "ChemistryUtilities.hh"

namespace Amanzi {
namespace AmanziChemistry {
namespace utilities {

/* ******************************************************************
* Remove spaces
****************************************************************** */
void RemoveLeadingAndTrailingWhitespace(std::string* line)
{
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
}

}  // namespace utilities
}  // namespace AmanziChemistry
}  // namespace Amanzi
