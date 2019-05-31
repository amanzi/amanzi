/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a PK integrated using
Explicit.
------------------------------------------------------------------------- */

#ifndef ATS_PK_EXPLICIT_DEFAULT_HH_
#define ATS_PK_EXPLICIT_DEFAULT_HH_

#include "Teuchos_TimeMonitor.hpp"
#include "Explicit_TI_RK.hh"
#include "PK_Explicit.hh"
#include "TreeVector.hh"
#include "Epetra_Vector.h"

namespace Amanzi {

class TreeVector;
class PK_Explicit_Default: public PK_Explicit<TreeVector>{

public:
PK_Explicit_Default(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& glist,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& solution):
  PK(pk_tree, glist, S, solution)
  {
      // name the PK
  name_ = pk_tree.name();
  boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(name_,"->");
  if (res.end() - name_.end() != 0) boost::algorithm::erase_head(name_, res.end() - name_.begin());


  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(glist, "PKs");

  if (pks_list->isSublist(name_)) {
    plist_ = Teuchos::sublist(pks_list, name_); 
  }else{
    std::stringstream messagestream;
    messagestream << "There is no sublist for PK "<<name_<<"in PKs list\n";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // THIS MAY BE CALLED MORE THAN ONCE!
  //name_ = plist_->get<std::string>("PK name");

  // set up the VerboseObject
  vo_ = Teuchos::rcp(new VerboseObject(name_, *plist_));

  }

// Virtual destructor
virtual ~PK_Explicit_Default(){};

    // Default implementations of PK methods.
  // -- setup
    virtual void Setup(const Teuchos::Ptr<State>& S);

    // -- initialize
    virtual void Initialize(const Teuchos::Ptr<State>& S);

    // -- Choose a time step compatible with physics.
    virtual double get_dt();

    virtual void set_dt(double dt);

    // -- Advance from state S0 to state S1 at time S0.time + dt.
    virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  protected: //data  timestep control
    double dt_;
    Teuchos::RCP<Explicit_TI::RK<TreeVector> > time_stepper_;

    // timing
    Teuchos::RCP<Teuchos::Time> step_walltime_;

    // solution at the old timestep
    Teuchos::RCP<TreeVector> solution_old_;

};



  
} // namespace

#endif
