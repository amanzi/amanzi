/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  PK for coupling of Transport_PK and Chemestry_PK on multiple meshes(domains)
*/

#ifndef AMANZI_COUPLEDREACTIVETRANSPORT_PK_HH_
#define AMANZI_COUPLEDREACTIVETRANSPORT_PK_HH_

#include "Teuchos_RCP.hpp"

#include "PK.hh"
#include "PK_Factory.hh"
#include "PK_MPC.hh"
#include "PK_MPCAdditive.hh"
#include "ReactiveTransport_PK.hh"
#include "Transport_PK.hh"
#include "Chemistry_PK.hh"


namespace Amanzi {

class CoupledReactiveTransport_PK : public ReactiveTransport_PK {
 public:
  CoupledReactiveTransport_PK(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln);

  ~CoupledReactiveTransport_PK() {
  }

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt);

  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);

  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  std::string name() { return "coupled reactive transport";}

 protected:  

  bool AdvanceChemistry_(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                         const Epetra_MultiVector& mol_dens,
                         Teuchos::RCP<Epetra_MultiVector> tcc_copy,
                         double t_old, double t_new, bool reinit);
  
  void cast_sub_pks_();
  int num_sub_domains_;
  Teuchos::RCP<PK_MPC> tranport_pk_;
  Teuchos::RCP<PK_MPC> chemistry_pk_;
  std::vector<Teuchos::RCP<Transport::Transport_PK> > tranport_pk_sub_;
  std::vector<Teuchos::RCP<AmanziChemistry::Chemistry_PK> > chemistry_pk_sub_;
  
  // factory registration
  static RegisteredPKFactory<CoupledReactiveTransport_PK> reg_;
};

}  // namespace Amanzi
#endif
