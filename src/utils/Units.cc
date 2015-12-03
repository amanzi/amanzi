/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#include <string>

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

