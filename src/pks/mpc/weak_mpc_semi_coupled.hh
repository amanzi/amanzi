

#ifndef WEAK_MPC_SEMI_COUPLED_HH_
#define WEAK_MPC_SEMI_COUPLED_HH_

#include "pk_physical_bdf_base.hh"
#include "weak_mpc.hh"
#include "mpc.hh"
namespace Amanzi {



namespace Flow {
  class PKPhysicalBDFBase;
}


  class WeakMPCSemiCoupled : public WeakMPC {
  public:
    WeakMPCSemiCoupled(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 Teuchos::ParameterList& FElist,
			  const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(plist, FElist, soln),
      WeakMPC(plist, FElist, soln)  {}; 
      //MPC<PK>(plist, FElist, soln) {};

    virtual ~WeakMPCSemiCoupled() {}
    
//    virtual double get_dt();
    virtual bool advance (double dt);
    
    virtual void setup(const Teuchos::Ptr<State>& S);
  //   virtual void initialize(const Teuchos::Ptr<State>& S);

 private :
    static RegisteredPKFactory<WeakMPCSemiCoupled> reg_;
  };

  
}



#endif
