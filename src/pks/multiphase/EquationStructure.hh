/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Collection of names of evaluators f1_ik, f2_ik, g1_ik, g2_ik, 
  and h_k for each equation in a system of equations:

    d h_k / dt - Sum_i { s1_ik div [K f1_ik grad g1_ik] } 
               - Sum_i { s2_ik div [D f2_ik grad g2_ik] } = Q_k

  subject to the constraint min(F, G) = 0. Here index i corresponds to
  a phase and s1_ik, s2_ik are scalar factors. Temporarily, we separate
  advective fluxes (1-terms) from molecular diffusion fluxes (2-terms).
*/

#ifndef AMANZI_EQUATION_STRUCTURE_PK_HH_
#define AMANZI_EQUATION_STRUCTURE_PK_HH_

#include <string>
#include <vector>

// Amanzi
#include "Key.hh"

namespace Amanzi {
namespace Multiphase {

typedef std::pair<Key, Key> EquationTerm;

struct EquationStructure {
 public:
  std::vector<EquationTerm> advection, diffusion;
  std::vector<double> adv_factors, diff_factors;

  Key storage;
  EquationTerm constraint;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif

