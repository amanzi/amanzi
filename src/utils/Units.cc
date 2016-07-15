/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#include <string>

#include <boost/algorithm/string.hpp>

#include "Units.hh"

extern bool Amanzi::Utils::concentration_mol_liter = true;

namespace boost {
namespace units {

std::string symbol_string(const Amanzi::Utils::concentration_amanzi&) {
  return (Amanzi::Utils::concentration_mol_liter) ? "mol/L" : "mol/m^3";
}

std::string name_string(const Amanzi::Utils::concentration_amanzi&) {
  return "concentration";
}

}  // namespace units
}  // namespace boost


namespace Amanzi {
namespace Utils {

/* ******************************************************************
* Convertion routine: time. We assume that input is in seconds.
****************************************************************** */
double Units::convert_time(double t, const std::string& unit)
{
  std::string u(unit);
  boost::algorithm::to_lower(u);

  double val(t);
  if (u == "y" || u == "yr" || u == "year") {
    val /= (365.25 * 24.0 * 3600.0); 
  } else if (u == "d" || u == "day") {
    val /= (24.0 * 3600.0); 
  } else if (u == "m" || u == "month") {
    val /= (365.25 * 2.0 * 3600.0); 
  }

  return val;
}

}  // namespace Utils
}  // namespace Amanzi

