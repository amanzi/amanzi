/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_BDF_BASE_HH_
#define AMANZI_PK_PHYSICAL_BDF_BASE_HH_


#include "PK_default_base.hh"
#include "PK_physical_base.hh"
#include "FnTimeIntegratorPK.hh"

namespace Amanzi {

class PKPhysicalBDFBase : public FnTimeIntegratorPK, public PKPhysicalBase {

public:
  PKPhysicalBDFBase(){};

  PKPhysicalBDFBase(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& glist,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(pk_tree, glist, S, soln),
    PKPhysicalBase(pk_tree, glist, S, soln),
    FnTimeIntegratorPK(pk_tree, glist, S, soln){};
    
// Virtual destructor
  virtual ~PKPhysicalBDFBase(){};


};

} // namespace

#endif

