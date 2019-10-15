/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "dbc.hh"

#include <sstream>
#include <cstring>

namespace DBC {

Assertion::Assertion(const char* assertion, const char* filename,
                     unsigned int line_number)
  : assertion_(assertion), filename_(filename), line_number_(line_number){};


const char*
Assertion::what() const noexcept
{
  std::ostringstream message;
  message << "Assertion: \"" << assertion_ << "\" failed in file: " << filename_
          << ", at line: " << line_number_ << std::endl;

  auto s = message.str();
  char* tmp = strcpy(new char[s.size() + 1], s.c_str());
  return tmp;
}


void
amanzi_assert(const char* cond, const char* file, unsigned int line)
{
  Exceptions::amanzi_throw(Assertion(cond, file, line));
}

} // namespace DBC
