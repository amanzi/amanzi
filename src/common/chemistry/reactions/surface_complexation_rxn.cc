/*
  Chemistry 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Class for surface complexation reaction

  Notes:
  - Each instance of this class should contain a single unique
    surface site (e.g. >FeOH) and ALL surface complexes associated
    with that site!
*/

#include "surface_complexation_rxn.hh"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "matrix_block.hh"
#include "chemistry_exception.hh"

#include "VerboseObject.hh"
#include "exceptions.hh"

namespace Amanzi {
namespace AmanziChemistry {

SurfaceComplexationRxn::SurfaceComplexationRxn()
    : use_newton_solve_(false) {
  surface_site_.clear();
  surface_complexes_.clear();
  // dSx_dmi_.clear();
}

SurfaceComplexationRxn::SurfaceComplexationRxn(
    SurfaceSite* surface_sites,
    const std::vector<SurfaceComplex>& surface_complexes) {
  // surface site
  surface_site_.push_back(*surface_sites);

  // surface complexes
  for (auto it = surface_complexes.begin(); it != surface_complexes.end(); ++it) {
    surface_complexes_.push_back(*it);
  }
}


SurfaceComplexationRxn::SurfaceComplexationRxn(SurfaceSite surface_sites) {
  // surface site
  surface_site_.push_back(surface_sites);
  surface_complexes_.clear();
}


void SurfaceComplexationRxn::AddSurfaceComplex(SurfaceComplex surface_complex) {
  surface_complexes_.push_back(surface_complex);
}


void SurfaceComplexationRxn::SetNewtonSolveFlag() {
  std::vector<SurfaceComplex>::const_iterator srfcplx = surface_complexes_.begin();
  const double tolerance = 1.e-20;
  for (; srfcplx != surface_complexes_.end(); srfcplx++) {
    if (std::fabs(srfcplx->free_site_stoichiometry() - 1.) > tolerance) {
      use_newton_solve_ = true;
      break;
    }
  }
}


void SurfaceComplexationRxn::UpdateSiteDensity(double site_density) {
  surface_site_.at(0).UpdateSiteDensity(site_density);
}


void SurfaceComplexationRxn::Update(const std::vector<Species>& primarySpecies)
{
  // see pflotran source: surface_complexation.F90:694, subroutine RTotalSorbEqSurfCplx1
  const double site_density = (surface_site_[0]).SiteDensity();

  bool one_more = false;
  double damping_factor = 1.0;
  int iterations(0), max_iterations(100);

  while (iterations < max_iterations) {
    iterations++;

    // Initialize total to free site concentration
    double free_site_concentration = (surface_site_[0]).free_site_concentration();

    // Update surface complex concentrations; Add to total
    double total = free_site_concentration;
    for (auto it = surface_complexes_.begin(); it != surface_complexes_.end(); ++it) {
      it->Update(primarySpecies, (surface_site_[0]));
      total += it->free_site_stoichiometry() * it->surface_concentration();
    }

    if (one_more) break;

    if (/*use_newton_solve_*/true) {
      double residual = site_density - total;
      double dresidual_dfree_site_conc = 1.0;
      for (auto it = surface_complexes_.begin(); it != surface_complexes_.end(); ++it) {
        dresidual_dfree_site_conc += it->free_site_stoichiometry() *
            it->surface_concentration() / free_site_concentration;
      }

      double dfree_site_conc = residual / dresidual_dfree_site_conc;

      if (iterations > 100) {
        // excessive iterations, try damping the update
        damping_factor = 0.5;
      }
      free_site_concentration += damping_factor * dfree_site_conc;
      double tolerance = 1.e-12;
      double rel_change_in_free_site_conc = std::fabs(dfree_site_conc / free_site_concentration);
      if (rel_change_in_free_site_conc < tolerance) {
        one_more = true;
      }
    } else {
      total = total / free_site_concentration;
      free_site_concentration = site_density / total;
      one_more = true;
    }

    // update free site concentration
    (surface_site_[0]).set_free_site_concentration(free_site_concentration);
  }

  if (iterations == max_iterations) {
    std::ostringstream error_stream;
    error_stream << "SurfaceComplexationRxn::Update(): \n"
                 << "loop reached max_iterations: " << iterations << std::endl;
    Exceptions::amanzi_throw(ChemistryMaxIterationsReached(error_stream.str()));
  }
}


void SurfaceComplexationRxn::AddContributionToTotal(std::vector<double> *total) {
  // see pflotran source: surface_complexation.F90:825, subroutine RTotalSorbEqSurfCplx1
  for (auto it = surface_complexes_.begin(); it != surface_complexes_.end(); ++it) {
    it->AddContributionToTotal(total);
  }
}


void SurfaceComplexationRxn::AddContributionToDTotal(
    const std::vector<Species>& primarySpecies,
    MatrixBlock* dtotal) {
  // see pflotran source: surface_complexation.F90:773, subroutine RTotalSorbEqSurfCplx1

  // All referenced equations #s are from the pflotran chemistry implementation
  // document by Peter Lichtner

  // Eq. 2.3-47c
  unsigned int num_primary_species = primarySpecies.size();
  std::vector<double> nu_li_nu_i_Si(num_primary_species, 0.0);

  double sum_nu_i_sq_Si = 0.;
  for (auto it = surface_complexes_.begin(); it != surface_complexes_.end(); ++it) {
    double tempd = it->free_site_stoichiometry() * it->surface_concentration();
    for (int icomp = 0; icomp < it->ncomp(); icomp++) {
      // sum of nu_li * nu_i * S_i
      nu_li_nu_i_Si.at(it->species_id(icomp)) += it->stoichiometry(icomp) * tempd;
    }
    // sum of nu_i^2 * S_i
    sum_nu_i_sq_Si += it->free_site_stoichiometry() * tempd;  // (free_site_stoich*surf_conc)
  }

  // complete the denominator within the brackets
  sum_nu_i_sq_Si /= surface_site_.at(0).free_site_concentration();
  double Sx_plus_sum_nu_i_sq_Si = 1. + sum_nu_i_sq_Si;
  for (int i = 0; i < num_primary_species; ++i) {
    nu_li_nu_i_Si.at(i) *= -1.0;
    nu_li_nu_i_Si.at(i) /= Sx_plus_sum_nu_i_sq_Si;
    // convert from dlogm to dm
    nu_li_nu_i_Si.at(i) /= primarySpecies.at(i).molality();
  }

  for (auto srfcplx = surface_complexes_.begin();
       srfcplx != surface_complexes_.end(); ++srfcplx) {
    double surface_concentration = srfcplx->surface_concentration();
    double nui_Si_over_Sx =
        srfcplx->free_site_stoichiometry() * surface_concentration /
        surface_site_.at(0).free_site_concentration();
    for (int icomp = 0; icomp < srfcplx->ncomp(); icomp++) {
      int primary_species_id_i = srfcplx->species_id(icomp);
      // 2.3-47c converted to non-log form
      double dSi_mi = (srfcplx->stoichiometry(icomp) * surface_concentration /
                                primarySpecies.at(primary_species_id_i).molality()) +
          nu_li_nu_i_Si.at(icomp) * nui_Si_over_Sx;
      for (int jcomp = 0; jcomp < srfcplx->ncomp(); jcomp++) {
        // 2.3-48a converted to non-log form
        double dPsij_dmi = dSi_mi * srfcplx->stoichiometry(jcomp);
        // NOTE: is this the correct i,j indexing...?
        dtotal->AddValue(srfcplx->species_id(jcomp), primary_species_id_i, dPsij_dmi);
      }
    }
  }
}


void SurfaceComplexationRxn::DisplaySite(const Teuchos::Ptr<VerboseObject> vo) const {
  std::vector<SurfaceSite>::const_iterator site;
  for (site = surface_site_.begin(); site != surface_site_.end(); site++) {
    site->Display(vo);
  }
}


void SurfaceComplexationRxn::DisplayComplexes(const Teuchos::Ptr<VerboseObject> vo) const {
  for (auto it = surface_complexes_.begin(); it != surface_complexes_.end(); ++it) {
    it->Display(vo);
  }
}


void SurfaceComplexationRxn::Display(const Teuchos::Ptr<VerboseObject> vo) const {
  DisplaySite(vo);
  DisplayComplexes(vo);
}


void SurfaceComplexationRxn::display(const Teuchos::Ptr<VerboseObject> vo) const {
  DisplaySite(vo);
  DisplayComplexes(vo);
}


void SurfaceComplexationRxn::DisplayResultsHeader(const Teuchos::Ptr<VerboseObject> vo) const {
  std::stringstream message;
  message << std::setw(7) << "---" << std::endl;
  vo->Write(Teuchos::VERB_HIGH, message.str());
}


void SurfaceComplexationRxn::DisplayResults(const Teuchos::Ptr<VerboseObject> vo) const {
  surface_site_[0].DisplayResultsHeader(vo);
  std::vector<SurfaceSite>::const_iterator site;
  for (site = surface_site_.begin();
       site != surface_site_.end(); site++) {
    site->DisplayResults(vo);
  }
  vo->Write(Teuchos::VERB_HIGH, "\n");
  surface_complexes_[0].DisplayResultsHeader(vo);
  for (auto it = surface_complexes_.begin(); it != surface_complexes_.end(); ++it) {
    it->DisplayResults(vo);
  }
}

}  // namespace AmanziChemistry
}  // namespace Amanzi
