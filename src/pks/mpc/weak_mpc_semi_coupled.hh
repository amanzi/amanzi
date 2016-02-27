#ifndef WEAK_MPC_SEMI_COUPLED_HH_
#define WEAK_MPC_SEMI_COUPLED_HH_

#include "pk_physical_bdf_base.hh"
#include "weak_mpc.hh"
#include "mpc.hh"

namespace Amanzi {

class WeakMPCSemiCoupled : public WeakMPC {
 public:
  WeakMPCSemiCoupled(const Teuchos::RCP<Teuchos::ParameterList>& plist,
		     Teuchos::ParameterList& FElist,
		     const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),
    WeakMPC(plist, FElist, soln){}; 

  virtual bool advance (double dt);
  virtual void setup(const Teuchos::Ptr<State>& S);
  
 private :
  static RegisteredPKFactory<WeakMPCSemiCoupled> reg_;
  Key domain_ss, domain_sf;

};

  
}



#endif
