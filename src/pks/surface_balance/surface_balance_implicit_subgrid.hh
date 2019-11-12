/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (coonet @ ornl.gov)
   
 ------------------------------------------------------------------------- */

//! ImplicitSubgrid: an implicit PK for surface balance with multiple subgrid patches.

/*!

Sets up a collection of patches, for portions of the column covered in snow,
ponded water, and vegetated/bare ground.  The surface energy balance on these
area weighted patches are individually calculated then averaged to form the
total quantities.  All down- and up-scaling of relevant quantities are done
through the area weighting, which is calculated by a minimum threshold in snow
and a depression depth/geometry-based approach for water.  All snow is assumed
to first cover water (likely ice), then cover land, as both water and snow
prefer low-lying depressions due to gravity- and wind-driven redistributions,
respectively.

*/


#ifndef PK_SURFACE_BALANCE_IMPLICIT_SUBGRID_HH_
#define PK_SURFACE_BALANCE_IMPLICIT_SUBGRID_HH_

#include "PK_Factory.hh"
#include "surface_balance_base.hh"


namespace Amanzi {
namespace SurfaceBalance {

class ImplicitSubgrid : public SurfaceBalanceBase {

public:

  ImplicitSubgrid(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S) override;

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u) override;


  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) override;
  
  // computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new,  const Teuchos::RCP<State>& S) override;

 protected:
  // multiple primary variables
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_snow_dens_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_snow_death_rate_;

  Key snow_dens_key_;
  Key snow_age_key_;
  Key new_snow_key_;
  Key snow_source_key_;
  Key snow_death_rate_key_;
  Key area_frac_key_;

 private:
  // factory registration
  static RegisteredPKFactory<ImplicitSubgrid> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
