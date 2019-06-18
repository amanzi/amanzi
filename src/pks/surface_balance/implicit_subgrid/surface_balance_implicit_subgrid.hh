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
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitState(double t_old, double t_new,  const Teuchos::RCP<State>& S) {}

  // -- Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}

  virtual void set_dt(double dt) {dt_ = dt;}

 protected:
  // multiple primary variables
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_snow_dens_;

  Key snow_dens_key_;
  Key snow_age_key_;
  Key new_snow_key_;
  Key snow_source_key_;

 private:
  // factory registration
  static RegisteredPKFactory<ImplicitSubgrid> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
