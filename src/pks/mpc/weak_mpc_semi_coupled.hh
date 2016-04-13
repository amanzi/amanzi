#ifndef WEAK_MPC_SEMI_COUPLED_HH_
#define WEAK_MPC_SEMI_COUPLED_HH_

#include "pk_physical_bdf_base.hh"
//#include "weak_mpc.hh"
#include "mpc.hh"
#include "PK.hh"

namespace Amanzi {
  
class WeakMPCSemiCoupled : public MPC<PK> {
public:
  /*
    WeakMPCSemiCoupled(const Teuchos::RCP<Teuchos::ParameterList>& plist,
    Teuchos::ParameterList& FElist,
    const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),
    WeakMPC(plist, FElist, soln){}; 
  */
  WeakMPCSemiCoupled(const Teuchos::RCP<Teuchos::ParameterList>& plist,
		     Teuchos::ParameterList& FElist,
		     const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln),MPC<PK>(){
    plist_ = plist;
    Teuchos::RCP<Teuchos::ParameterList> pks_list;
    Teuchos::Array<std::string> pk_order = plist->get<Teuchos::Array<std::string> >("PKs order");

    Teuchos::Array<int> pk_proc = plist->get<Teuchos::Array<int> >("PKs process index");
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    Teuchos::Array<std::string> loc_list;
    int pk_start, pk_end;
   
   
    pk_start = pk_proc[rank];
    pk_end = pk_proc[rank+1];
   
    loc_list.push_back("PKG1");
    for(int i=pk_start; i<pk_end; i++){
      std::stringstream name;
      name << "PK" << i;
      loc_list.push_back(name.str());
    }
    /*    if(rank ==0)
      for(int i=0; i<11; i++){
	std::stringstream name;
	name << "PK" << i;
	loc_list.push_back(name.str());
      }
    else if(rank ==1)
      for(int i=11; i<21; i++){
        std::stringstream name;
        name << "PK" << i;
        loc_list.push_back(name.str());
      }
    */
    plist_->set("PKs order", loc_list);
    
    MPC<PK>::init_(plist_,FElist, soln);
    
    
  };
  virtual double get_dt();
  virtual bool valid_step();
  virtual bool advance (double dt);
  virtual void setup(const Teuchos::Ptr<State>& S);
  bool  CoupledSurfSubsurf3D(double dt);

  bool CoupledSurfSubsurfColumns(double dt);
private :
  static RegisteredPKFactory<WeakMPCSemiCoupled> reg_;
  unsigned numPKs_;
  Key coupling_key_ ;
  
};

  
}



#endif
