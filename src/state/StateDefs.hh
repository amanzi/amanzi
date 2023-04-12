/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
   State

  Some basic typedefs for State and company.
*/

#ifndef AMANZI_STATE_DEFS_HH_
#define AMANZI_STATE_DEFS_HH_

#include <string>
#include <vector>
#include <set>

#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "Tag.hh"

namespace Amanzi {

typedef bool NullFactory; // placeholder object for no factory required

namespace Tags {
static const Tag DEFAULT("");
static const Tag CURRENT("current");
static const Tag INTER("inter");
static const Tag NEXT(""); // an alias used by ATS
static const Tag COPY("copy");
} // namespace Tags

enum class Evaluator_kind { PRIMARY, SECONDARY, INDEPENDENT, OTHER };

inline
std::string to_string(const Evaluator_kind kind) {
  switch(kind) {
    case(Evaluator_kind::PRIMARY): return "primary";
    case(Evaluator_kind::SECONDARY): return "secondary";
    case(Evaluator_kind::INDEPENDENT): return "independent";
    default: return "other";
  }
}


// Type enumerator for derivatives in models.
template <int>
struct Deriv {};

} // namespace Amanzi

#endif
