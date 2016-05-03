#ifndef WEAK_MPC_SEMI_COUPLED_HH_
#define WEAK_MPC_SEMI_COUPLED_HH_

#include "pk_physical_bdf_base.hh"
//#include "weak_mpc.hh"
#include "mpc.hh"
#include "PK.hh"

namespace Amanzi {
  
class WeakMPCSemiCoupled : public MPC<PK> {
public:
  WeakMPCSemiCoupled(Teuchos::Ptr<State> S,const Teuchos::RCP<Teuchos::ParameterList>& plist,
		     Teuchos::ParameterList& FElist,
		     const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(S, plist, FElist, soln),MPC<PK>(),FElist_loc(FElist){
    
    plist_ = plist;
    S_loc = S;
    generalize_inputspec();
    MPC<PK>::init_(S,plist_,FElist_loc, soln);
    
  };


  virtual double get_dt();
  virtual bool valid_step();
  virtual bool advance (double dt);
  virtual void setup(const Teuchos::Ptr<State>& S);
  virtual void initialize(const Teuchos::Ptr<State>& S);

  void generalize_inputspec();
  bool CoupledSurfSubsurf3D(double dt);

  bool CoupledSurfSubsurfColumns(double dt);
  

private :
  static RegisteredPKFactory<WeakMPCSemiCoupled> reg_;
  unsigned numPKs_;
  Key coupling_key_ ;
  Teuchos::ParameterList& FElist_loc;
  Teuchos::Ptr<State> S_loc;
};

  
}



#endif
