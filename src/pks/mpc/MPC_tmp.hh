/* -------------------------------------------------------------------------
Amanzi

Licenses: see $ATS_DIR/COPYRIGHT, $ASCEM_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the Base MPC class.  A multi process coordinator (MPC)
is a PK which coordinates other PKs.  Each of these coordinated PKs
may be MPCs themselves, or physical PKs.  Note this does NOT provide a
full implementation of PK -- it does not supply the Advance() method.
Therefore this class cannot be instantiated, but must be inherited by
derived classes which finish supplying the functionality.  Instead,
this provides the data structures and methods (which may be overridden
by derived classes) for managing multiple PKs.

Most of these methods simply loop through the coordinated PKs, calling their
respective methods.
------------------------------------------------------------------------- */

#ifndef AMANZI_MPC_HH_
#define AMANZI_MPC_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {

template <class PK_t>
class MPCTmp : public PK {

public:
  MPCTmp(Teuchos::ParameterList& pk_tree,
         const Teuchos::RCP<Teuchos::ParameterList>& global_list,
         const Teuchos::RCP<State>& S,
         const Teuchos::RCP<TreeVector>& soln);

  ~MPCTmp(){ std::cout<<"Destructor "<<name()<<"\n";}

  // PK methods
  // -- sets up sub-PKs
  virtual void Setup();

  // -- calls all sub-PK initialize() methods
  virtual void Initialize();

  // -- loops over sub-PKs
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics();

  // -- identifier accessor
  std::string name() const { return name_; }

 protected:
  // identifier
  std::string name_;

  // list of the PKs coupled by this MPC
  typedef std::vector<Teuchos::RCP<PK_t> > SubPKList;
  SubPKList sub_pks_;

  // single solution vector for the coupled problem
  Teuchos::RCP<TreeVector> solution_;

  // plists
  Teuchos::RCP<Teuchos::ParameterList> global_list_;
  Teuchos::RCP<Teuchos::ParameterList> my_list_;
  Teuchos::ParameterList pk_tree_;

  // states
  Teuchos::RCP<State> S_;
  
};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
MPCTmp<PK_t>::MPCTmp(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& soln) :
    pk_tree_(pk_tree),
    global_list_(global_list),
    S_(S),
    solution_(soln)
{

  // name the PK
  name_ = pk_tree.name();

  const char* result = name_.data();
  while ((result = std::strstr(result, "->")) != NULL) {
    result += 2;
    name_ = result;   
  }

  //std::cout<<pk_tree<<"\n";
  //std::cout<<*global_list<<"\n";

  // get my ParameterList
  my_list_ = Teuchos::sublist(Teuchos::sublist(global_list_, "PKs"), name_);

  Teuchos::RCP<Teuchos::ParameterList> plist;
  //     Teuchos::sublist(Teuchos::sublist(global_list, "PKs"), name_);

  if (global_list_->isSublist("PKs")){
    plist =  Teuchos::sublist(global_list, "PKs");
  }

  std::vector<std::string> pk_name =  my_list_->get<Teuchos::Array<std::string> >("PKs order").toVector();

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  PKFactory pk_factory;
  // for (Teuchos::ParameterList::ConstIterator sub=plist->begin();
  //      sub!=plist->end(); ++sub) {
  for (int i=0; i<pk_name.size(); i++){
    //const std::string& sub_name = sub->first;
    const std::string& sub_name = pk_name[i];
    std::cout<<"sub_name "<<sub_name<<"\n";
    if (!plist->isSublist(sub_name)) {
      Errors::Message message("PK Tree: All entries in the PK tree must be ParameterLists");
      Exceptions::amanzi_throw(message);
    }
  }

  for (int i=0; i<pk_name.size(); i++){
    // Collect arguments to the constructor
    Teuchos::ParameterList& pk_sub_tree = pk_tree.sublist(pk_name[i]);
    Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
    solution_->PushBack(pk_soln);

    // create the PK
    Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(pk_sub_tree, global_list,
            S, pk_soln);
    Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype);
    sub_pks_.push_back(pk);
  }


}


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
void MPCTmp<PK_t>::Setup() {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->Setup();
  }
}


// -----------------------------------------------------------------------------
// Loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPCTmp<PK_t>::Initialize() {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->Initialize();
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit state method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPCTmp<PK_t>::CommitStep(double t_old, double t_new) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->CommitStep(t_old, t_new);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their CalculateDiagnostics method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPCTmp<PK_t>::CalculateDiagnostics() {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->CalculateDiagnostics();
  }
};


} // close namespace Amanzi

#endif
