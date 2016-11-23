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
    PK(FElist, plist, S, solution),
    MPC<PK>(FElist, plist, S, solution),FElist_loc(FElist){
    
    plist_ = plist;
    S_loc = S.ptr();
    generalize_inputspec();
    //    MPC<PK>::init_(S,plist_,FElist_loc, soln);
    
  };


  virtual double get_dt();
  virtual void set_dt(double dt);
  // virtual bool valid_step();
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit); //virtual bool advance (double dt);
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  void generalize_inputspec();
  bool CoupledSurfSubsurf3D(double t_old, double t_new, bool reinit);

  bool CoupledSurfSubsurfColumns(double t_old, double t_new, bool reinit);
  

private :
  static RegisteredPKFactory<WeakMPCSemiCoupled> reg_;
  unsigned numPKs_;
  static unsigned flag_star, flag_star_surf;
  Key coupling_key_ ;
  bool subcycle_key_ ;
  Teuchos::ParameterList& FElist_loc;
  Teuchos::Ptr<State> S_loc;
  double min_dt_, surf_dt_, sync_time_;
};

  
}



#endif
