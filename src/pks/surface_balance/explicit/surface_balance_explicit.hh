/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Adam Atchley, Satish Karra

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */

#ifndef PK_SURFACE_BALANCE_EXPLICIT_HH_
#define PK_SURFACE_BALANCE_EXPLICIT_HH_

#include "PK_Factory.hh"
#include "pk_physical_default.hh"

namespace Amanzi {
namespace SurfaceBalance {

class SurfaceBalanceExplicit : public PK_Physical_Default {

public:
  SurfaceBalanceExplicit(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- provide a timestep size
  virtual double get_dt() {
    return dt_;
  }

  virtual void set_dt(double dt) {dt_ = dt;}

  virtual std::string name(){ return "SurfaceBalanceExplicit";}

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S){};

  // -- Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) {};

  // -- advance via one of a few methods
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

 protected:
  // multiple primary variables
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_esource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wsource_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_w_v_source_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pvfe_wtemp_;

  double dt_;
  double min_wind_speed_;
  double snow_ground_trans_;
  double min_snow_trans_;

  Teuchos::RCP<const AmanziMesh::Mesh> subsurf_mesh_;

 private:
  // factory registration
  static RegisteredPKFactory<SurfaceBalanceExplicit> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
