#include <cmath>

#include <iostream>
#include <iomanip>
#include <sstream>

#include "ChemistryException.hpp"
#include "SurfaceComplexationRxn.hpp"
#include "Block.hpp"

#include "exceptions.hh"

SurfaceComplexationRxn::SurfaceComplexationRxn()
    : use_newton_solve_(false)
                  
{
  surface_site_.clear();
  surface_complexes_.clear();
  dSx_dmi_.clear();
} 

SurfaceComplexationRxn::SurfaceComplexationRxn(
                            SurfaceSite *surface_sites,
                            std::vector<SurfaceComplex> surface_complexes)
{
  // surface site
  surface_site_.push_back(*surface_sites);

  // surface complexes
  for (std::vector<SurfaceComplex>::const_iterator i = surface_complexes.begin(); 
       i != surface_complexes.end(); i++) {
    surface_complexes_.push_back(*i);
  } 
}

SurfaceComplexationRxn::SurfaceComplexationRxn(SurfaceSite surface_sites)
{
  // surface site
  surface_site_.push_back(surface_sites);
  surface_complexes_.clear();
}

SurfaceComplexationRxn::~SurfaceComplexationRxn() 
{
}

void SurfaceComplexationRxn::AddSurfaceComplex(SurfaceComplex surface_complex)
{
  surface_complexes_.push_back(surface_complex);
}

void SurfaceComplexationRxn::SetNewtonSolveFlag(void) 
{
  std::vector<SurfaceComplex>::const_iterator srfcplx = 
    surface_complexes_.begin();
  const double tolerance = 1.e-20;
  for (; srfcplx != surface_complexes_.end(); srfcplx++) {
    if (fabs(srfcplx->free_site_stoichiometry() - 1.) > tolerance) {
      set_use_newton_solve(true);
      break;
    }
  }
} // end SetNewtonSolveFlag

void SurfaceComplexationRxn::Update(const std::vector<Species> primarySpecies) 
{
  const double site_density = (surface_site_[0]).SiteDensity();

  bool one_more = false;
  int max_iterations = 5000;
  int iterations = 0;
  while(iterations < max_iterations) {
    iterations++;
    // Initialize total to free site concentration
    double free_site_concentration = (surface_site_[0]).free_site_concentration();
    double total = free_site_concentration;
    // Update surface complex concentrations; Add to total

    for (std::vector<SurfaceComplex>::iterator srfcplx = 
      surface_complexes_.begin(); 
      srfcplx != surface_complexes_.end(); srfcplx++) {
      srfcplx->Update(primarySpecies,(surface_site_[0]));
      total += srfcplx->free_site_stoichiometry() *
               srfcplx->surface_concentration();
    }

    if (one_more) break;

    if (/*use_newton_solve()*/true) {
      double residual = site_density - total;
      double dresidual_dfree_site_conc = 1.;
      std::vector<SurfaceComplex>::iterator srfcplx = 
        surface_complexes_.begin();
      for (; srfcplx != surface_complexes_.end(); srfcplx++) {
        dresidual_dfree_site_conc += srfcplx->free_site_stoichiometry() *
                                srfcplx->surface_concentration() / 
                                free_site_concentration;
      }
      double dfree_site_conc = residual / dresidual_dfree_site_conc;
      free_site_concentration += dfree_site_conc;
      double tolerance = 1.e-12;
      if (fabs(dfree_site_conc / free_site_concentration) < tolerance) {
        one_more = true;
      }
    }
    else {
      total = total / free_site_concentration;
      free_site_concentration = site_density / total;
      one_more = true;
    }

    // update free site concentration
    (surface_site_[0]).set_free_site_concentration(free_site_concentration);

  }
  if (iterations == max_iterations) {
    std::ostringstream error_stream;
    error_stream << "ERROR: SurfaceComplexationRxn::Update(): \n";
    error_stream << "ERROR: loop reached max_iterations: " << iterations << std::endl;
    Exceptions::amanzi_throw(ChemistryException(error_stream.str()));
  }
  
} // end Update()

void SurfaceComplexationRxn::AddContributionToTotal(std::vector<double> *total) 
{
  
  for (std::vector<SurfaceComplex>::iterator srfcplx = 
       surface_complexes_.begin(); 
       srfcplx != surface_complexes_.end(); srfcplx++) {
    srfcplx->AddContributionToTotal(total);
  }
} // end AddContributionToTotal()

void SurfaceComplexationRxn::AddContributionToDTotal(
                                   const std::vector<Species> primarySpecies,
                                   Block *dtotal) 
{
  // All referenced equations #s are from the pflotran chemistry implementation
  // document by Peter Lichtner
  
  // Eq. 2.3-47c
  int num_primary_species = (int)primarySpecies.size();
  double *nu_li_nu_i_Si = new double[num_primary_species];
  for (int i = 0; i < num_primary_species; i++)
    nu_li_nu_i_Si[i] = 0.;
  double sum_nu_i_sq_Si = 0.;
  for (std::vector<SurfaceComplex>::iterator srfcplx = 
       surface_complexes_.begin(); 
       srfcplx != surface_complexes_.end(); srfcplx++) {
    double tempd = srfcplx->free_site_stoichiometry() *
                   srfcplx->surface_concentration();
    for (int icomp = 0; icomp < srfcplx->ncomp(); icomp++) {
      // sum of nu_li * nu_i * S_i
      nu_li_nu_i_Si[srfcplx->species_id(icomp)] += 
                                   srfcplx->stoichiometry(icomp) * tempd;
    }
    // sum of nu_i^2 * S_i
    sum_nu_i_sq_Si += srfcplx->free_site_stoichiometry() *
                      tempd; // (free_site_stoich*surf_conc)
  }
  // complete the denominator within the brackets
  double Sx_plus_sum_nu_i_sq_Si = 1. + sum_nu_i_sq_Si;

  for (std::vector<SurfaceComplex>::iterator srfcplx = 
       surface_complexes_.begin(); 
       srfcplx != surface_complexes_.end(); srfcplx++) {
    double surface_concentration = srfcplx->surface_concentration();
    for (int icomp = 0; icomp < srfcplx->ncomp(); icomp++) {
      int primary_species_id_i = srfcplx->species_id(icomp);
      // 2.3-47c converted to non-log form
      double dSi_mi = surface_concentration /
                      primarySpecies[primary_species_id_i].molality() *
                      (srfcplx->stoichiometry(icomp) -
                       (nu_li_nu_i_Si[icomp] / Sx_plus_sum_nu_i_sq_Si));
      for (int jcomp = 0; jcomp < srfcplx->ncomp(); jcomp++) {
        // 2.3-48a converted to non-log form
        double dPsij_dmi = dSi_mi * srfcplx->stoichiometry(jcomp);
        dtotal->addValue(srfcplx->species_id(jcomp), primary_species_id_i, dPsij_dmi);
      }
    }
  }
  delete [] nu_li_nu_i_Si;

} // end AddContributionToDTotal()

void SurfaceComplexationRxn::DisplaySite(void) const
{
  std::vector<SurfaceSite>::const_iterator site;
  for (site = surface_site_.begin(); site != surface_site_.end(); site++) {
    site->Display();
  }
} // end DisplaySite()

void SurfaceComplexationRxn::DisplayComplexes(void) const
{
  std::vector<SurfaceComplex>::const_iterator complex;
  for (complex = surface_complexes_.begin(); 
       complex != surface_complexes_.end(); complex++) {
    complex->Display();
  }
} // end DisplayComplexes()

void SurfaceComplexationRxn::Display(void) const
{
  DisplaySite();
  DisplayComplexes();
} // end Display()

void SurfaceComplexationRxn::display(void) const
{
  DisplaySite();
  DisplayComplexes();
} // end display()

void SurfaceComplexationRxn::DisplayResultsHeader(void) const
{
  std::cout << std::setw(15) << "---"
            << std::endl;
} // end DisplayResultsHeader()

void SurfaceComplexationRxn::DisplayResults(void) const
{

  surface_site_[0].DisplayResultsHeader();
  std::vector<SurfaceSite>::const_iterator site;
  for (site = surface_site_.begin(); 
       site != surface_site_.end(); site++) {
    site->DisplayResults();
  }

  surface_complexes_[0].DisplayResultsHeader();
  std::vector<SurfaceComplex>::const_iterator complex;
  for (complex = surface_complexes_.begin(); 
       complex != surface_complexes_.end(); complex++) {
    complex->DisplayResults();
  }

} // end DisplayResults()
