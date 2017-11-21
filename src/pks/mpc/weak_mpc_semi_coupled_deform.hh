/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  This is a semi weak mpc which couples subsurface columns (no lateral flow) to surface_star system.
  surface_star system distributes water and energy and updates columns
  
  Algorithm:
  1 - Advance surface_star system (Pres*, Temp*)
  2 - Update subsurface columns' pressures and temperatures (Pres, Temp)
  3 - Advance the subsurface columns. Just loop over all the columns.
  4 - Update the surface_star system
  5 - commit state if timestep was successfull
  6 - take next timestep and repeart (steps 1-5)
  
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef WEAK_MPC_SEMI_COUPLED_DEFORM_HH_
#define WEAK_MPC_SEMI_COUPLED_DEFORM_HH_

#include "pk_physical_bdf_default.hh"
#include "mpc.hh"
#include "PK.hh"

namespace Amanzi {
  
class WeakMPCSemiCoupledDeform : public MPC<PK> {
 public:
   
  WeakMPCSemiCoupledDeform(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& solution);

  virtual double get_dt();
  virtual void set_dt(double dt);

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit); //virtual bool advance (double dt);
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  bool CoupledSurfSubsurfColumns(double t_old, double t_new, bool reinit);
  
  double FindVolumetricHead(double d, double delta_max, double delta_ex);
  double VolumetricHead(double x, double a, double b, double d);

  
private :
  static RegisteredPKFactory<WeakMPCSemiCoupledDeform> reg_;
  unsigned numPKs_;
  static unsigned flag_star, flag_star_surf;
  Key coupling_key_ ;
  bool subcycle_key_ ;
  bool sg_model_, dynamic_sg_model_;
};

  
}



#endif
