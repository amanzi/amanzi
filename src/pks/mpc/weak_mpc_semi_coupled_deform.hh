#ifndef WEAK_MPC_SEMI_COUPLED_DEFORM_HH_
#define WEAK_MPC_SEMI_COUPLED_DEFORM_HH_

//#include "pk_physical_bdf_base.hh"
#include "pk_physical_bdf_default.hh"
//#include "weak_mpc.hh"
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

  //  void generalize_inputspec(const Teuchos::Ptr<State>& S);
  bool CoupledSurfSubsurf3D(double t_old, double t_new, bool reinit);

  bool CoupledSurfSubsurfColumns(double t_old, double t_new, bool reinit);
  
  double FindVolumetricHead(double d, double delta_max, double delta_ex);
  double VolumetricHead(double x, double a, double b, double d);

  
private :
  static RegisteredPKFactory<WeakMPCSemiCoupledDeform> reg_;
  unsigned numPKs_;
  static unsigned flag_star, flag_star_surf;
  Key coupling_key_ ;
  bool subcycle_key_ ;
  

  bool sg_model_;
};

  
}



#endif
