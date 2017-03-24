#ifndef WEAK_MPC_SEMI_COUPLED_HH_
#define WEAK_MPC_SEMI_COUPLED_HH_

//#include "pk_physical_bdf_base.hh"
#include "pk_physical_bdf_default.hh"
//#include "weak_mpc.hh"
#include "mpc.hh"
#include "PK.hh"

namespace Amanzi {
  
class WeakMPCSemiCoupled : public MPC<PK> {
public:
  
   WeakMPCSemiCoupled(Teuchos::ParameterList& FElist,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution) :
     PK(), MPC<PK>(),FElist_loc(FElist){
     
     plist_ = plist;
    
     generalize_inputspec(S.ptr());
     
     MPC<PK>::init_(FElist_loc,plist_,S, solution);
     
   };


  virtual double get_dt();
  virtual void set_dt(double dt);
  // virtual bool valid_step();
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit); //virtual bool advance (double dt);
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  void generalize_inputspec(const Teuchos::Ptr<State>& S);
  bool CoupledSurfSubsurf3D(double t_old, double t_new, bool reinit);

  bool CoupledSurfSubsurfColumns(double t_old, double t_new, bool reinit);
  
  double FindVolumetricHead(double d, double delta_max, double delta_ex);
  double VolumetricHead(double x, double a, double b, double d);

private :
  static RegisteredPKFactory<WeakMPCSemiCoupled> reg_;
  unsigned numPKs_;
  static unsigned flag_star, flag_star_surf;
  Key coupling_key_ ;
  bool subcycle_key_ ;
  Teuchos::ParameterList& FElist_loc;
  
  double min_dt_, surf_dt_, sync_time_;
  double delta_ex_, delta_max_; // subgrid model's parameters
  bool sg_model_;
};

  
}



#endif
