/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

This is a single, static list of abbreviations used to shorten stdout IO,
particularly for headers and variable names.  Feel free to add to this.

*/

#pragma once

#include <vector>

namespace Amanzi {
namespace Keys {

// note, this is a vector and not a map to keep the order as provided here
static std::vector<std::pair<std::string, std::string>> abbreviations = {
  { "surface_star", "ss" },
  { "surface_column", "surf_col" },
  { "column", "col" },
  { "temperature", "temp" },
  { "internal_energy", "ie" },
  { "molar", "mol"},
};

} // namespace Keys
} // namespace Amanzi
