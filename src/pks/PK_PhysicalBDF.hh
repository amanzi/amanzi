/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_BDF_HH_
#define AMANZI_PK_PHYSICAL_BDF_HH_

#include "PK_BDF.hh"
#include "PK_Physical.hh"
#include "primary_variable_field_evaluator.hh"
#include "Debugger.hh"

namespace Amanzi {

class PK_PhysicalBDF : virtual public PK_Physical, public PK_BDF {

public:
  PK_PhysicalBDF() {};

  PK_PhysicalBDF(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln):
      PK_Physical(pk_tree, glist, S, soln) {};
   
  // Virtual destructor
  virtual ~PK_PhysicalBDF() {};
};

} // namespace Amanzi

#endif

