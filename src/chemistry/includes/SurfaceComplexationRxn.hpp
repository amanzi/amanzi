#ifndef __SurfaceComplexationRxn_hpp__
#define __SurfaceComplexationRxn_hpp__

#include <string>
#include <vector>
#include <cmath>

#include "Block.hpp"
#include "SurfaceComplex.hpp"
#include "SurfaceSite.hpp"

// Class for aqueous equilibrium complexation reaction

class SurfaceComplexationRxn {

 public:
  SurfaceComplexationRxn();
  SurfaceComplexationRxn(std::string s);
  SurfaceComplexationRxn(SurfaceSite *surface_sites,
                         std::vector<SurfaceComplex> surface_complexes);
  SurfaceComplexationRxn(SurfaceSite surface_sites);
  ~SurfaceComplexationRxn();

  // add complexes to the reaction
  void AddSurfaceComplex(SurfaceComplex surface_complex);

  // update sorbed concentrations
  void Update(const std::vector<Species>primarySpecies);
  // add stoichiometric contribution of complex to sorbed total
  void AddContributionToTotal(std::vector<double> *total);
  // add derivative of total with respect to free-ion to sorbed dtotal
  void AddContributionToDTotal(const std::vector<Species> primarySpecies,
                               Block *dtotal);
  // If the free site stoichiometry in any of the surface complexes
  // is not equal to 1., we must use Newton's method to solve for
  // the free site concentration.  This function determines if this
  // is the case.
  void SetNewtonSolveFlag(void);

  void display(void) const;
  void Display(void) const;
  void DisplaySite(void) const;
  void DisplayComplexes(void) const;
  void DisplayResultsHeader(void) const;
  void DisplayResults(void) const;

 protected:

   void set_use_newton_solve(const bool b) { this->use_newton_solve_ = b; };

   bool use_newton_solve(void) const { return this->use_newton_solve_; };

 private:
   std::vector<SurfaceComplex> surface_complexes_;
   std::vector<SurfaceSite> surface_site_;
   bool use_newton_solve_;

   std::vector<double> dSx_dmi_; // temporary storage for derivative calculations

};

#endif // __SurfaceComplexationRxn_hpp__
