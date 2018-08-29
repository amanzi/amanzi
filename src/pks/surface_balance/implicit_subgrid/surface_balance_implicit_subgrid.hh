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
#include "pk_physical_bdf_default.hh"

namespace Amanzi {
namespace SurfaceBalance {

class ImplicitSubgrid : public PK_PhysicalBDF_Default {

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

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // error monitor
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitState(double t_old, double t_new,  const Teuchos::RCP<State>& S) {}

  // -- Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {}

  virtual void set_dt(double dt) {dt_ = dt;}

 protected:
  // multiple primary variables
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_esource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wsource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_w_sub_source_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_e_sub_source_;

  bool eval_derivatives_;
  bool longwave_input_;
  Key sw_incoming_key_;
  
  double min_wind_speed_;       // wind speed of 0, under this model, would have 0 latent or sensible heat?
  double wind_speed_ref_ht_;    // reference height of the met data
  double roughness_bare_ground_;
  double roughness_snow_covered_ground_; // fetch lengths? Or elevation differences?  Or some other smoothness measure? [m]

  double snow_ground_trans_;    // snow depth at which soil starts to appear
  double min_snow_trans_;       // snow depth at which snow reaches no area coverage

  double dessicated_zone_thickness_; // max thickness of the zone over which
                                     // evaporation dessicates the soil, and
                                     // therefore vapor diffusion must act to
                                     // bring evaporated water to the surface.
                                     // A limiter on evaporation as the water
                                     // table drops below the surface.

  Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh_;
  Key domain_ss_;

 private:
  // factory registration
  static RegisteredPKFactory<ImplicitSubgrid> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
