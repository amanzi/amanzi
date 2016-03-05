/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_BDF_BASE_HH_
#define AMANZI_PK_PHYSICAL_BDF_BASE_HH_


#include "PK_Default.hh"
#include "PK_Physical.hh"
#include "PK_BDF.hh"

namespace Amanzi {

class PK_PhysicalBDF : public PK_BDF, public PK_Physical {

public:
  PK_PhysicalBDF(){};

  PK_PhysicalBDF(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln) :
    PK_Default(pk_tree, glist, S, soln),
    PK_Physical(pk_tree, glist, S, soln),
    PK_BDF(pk_tree, glist, S, soln){};

  PK_PhysicalBDF(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<TreeVector>& solution):
    PK_Default(plist, FElist, solution),
    PK_Physical(plist, FElist, solution),
    PK_BDF(plist, FElist, solution){};
    
// Virtual destructor
  virtual ~PK_PhysicalBDF(){};


};

} // namespace

#endif

