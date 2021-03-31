/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include <cstdlib>
#include <string>
#include <vector>

#include "exceptions.hh"

#include "beaker.hh"

namespace Amanzi {
namespace AmanziChemistry {

int Beaker::EnforceConstraint(
    BeakerState* state, const BeakerParameters& parameters,
    const std::vector<std::string>& names,
    const std::vector<double>& values)
{
  ResetStatus();
  UpdateParameters(parameters, 1.0);

  // initial guess
  for (int i = 0; i < ncomp_; i++) {
    state->total.at(i) = 0.0;
    state->free_ion.at(i) = 1.0e-9;

    if (names[i] == "total concentration") {
      state->total.at(i) = values[i];
    } else if (names[i] == "charge balance") {
      state->free_ion.at(i) = values[i];
    } else if (names[i] == "pH") {
      state->free_ion.at(i) = std::pow(10.0, -values[i]);
    }
  }

  CopyStateToBeaker(*state);

  // initialize to a large number (not necessary, but safe)
  double max_residual, max_rel_change(1.0e20), tolerance(1.0e-12);
  int max_rel_index, num_iterations(0);

  do {
    UpdateActivityCoefficients();
    UpdateEquilibriumChemistry();
    CalculateDTotal();

    jacobian_.Zero();

    for (int i = 0; i < ncomp_; i++) {
      if (names[i] == "total concentration") {
        residual_[i] = total_.at(i) - state->total.at(i);

        for (int j = 0; j < ncomp_; ++j) {
          jacobian_(i, j) += dtotal_(i, j);
        }

      } else if (names[i] == "charge balance") {
        residual_[i] = 0.0;
        for (int j = 0; j < ncomp_; ++j) {
          residual_[i] += primary_species().at(j).charge() * total_[j];

          for (int k = 0; k < ncomp_; ++k) {
            jacobian_(i, j) += primary_species().at(k).charge() * dtotal_(k, j);
          }
        }

      } else if (names[i] == "pH") {
        residual_[i] = 0.0;
        jacobian_(i, i) = 1.0;

      } else {
        Exceptions::amanzi_throw(ChemistryInvalidInput("Unknown geochemical constraint"));
      }
    }

    // scale and solve
    rhs_ = residual_;

    for (int i = 0; i < ncomp_; i++) {
      jacobian_.ScaleColumn(i, primary_species().at(i).molality());
    }

    lu_solver_.Solve(&jacobian_, &rhs_);
    num_iterations++;

    // calculate update truncating at a maximum of 5 in log space
    UpdateMolalitiesWithTruncation(5.0);
    // calculate maximum relative change in concentration over all species
    CalculateMaxRelChangeInMolality(&max_rel_change, &max_rel_index);

    max_residual = 0.0;
    for (int i = 0; i < ncomp_; i++) {
      max_residual = std::max(max_residual, std::fabs(residual_.at(i)));
    }

  } while (max_rel_change > tolerance &&
           num_iterations < max_iterations_);

  // for now, initialize total sorbed concentrations based on the current free
  // ion concentrations
  UpdateEquilibriumChemistry();
  CopyBeakerToState(state);
  status_.num_newton_iterations = num_iterations;
  if (max_rel_change < tolerance_) {
    status_.converged = true;
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
