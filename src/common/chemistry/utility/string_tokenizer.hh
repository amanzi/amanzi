/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ben Andre

  This class is a simple string tokenizer.
   - Found July 2003, unencumbered, on the web at:
     http:// www.thecodezone.com/diary/archives/000057.html
 
   - no longer available from origional source, but archived at:
     http://web.archive.org/web/20030810163805/http://www.thecodezone.com/diary/archives/000057.html
*/

#ifndef AMANZI_CHEMISTRY_STRING_TOKENIZER_HH_
#define AMANZI_CHEMISTRY_STRING_TOKENIZER_HH_

#include <string>
#include <vector>

namespace Amanzi {
namespace AmanziChemistry {

class StringTokenizer : public std::vector<std::string> {
 public:
  StringTokenizer() {};
  StringTokenizer(const std::string& source,
                  const std::string& delimiters = " \t\n") { tokenize(source, delimiters); }

  void tokenize(const std::string& source,
                const std::string& delimiters = " \t\n");

  // the following tokenizes, but places delimiters in list too - geh
  void tokenize_leave_delimiters(const std::string& source,
                                 const std::string& delimiters = " \t\n");
};

}  // namespace AmanziChemistry
}  // namespace Amanzi
#endif
