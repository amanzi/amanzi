/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Daniil Svyatsky, Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */
#ifndef AMANZI_PK_PHYSICAL_EXPLICIT_HH_
#define AMANZI_PK_PHYSICAL_EXPLICIT_HH_

#include "PK_Explicit.hh"
#include "PK_Physical.hh"

namespace Amanzi {

template <class Vector>
class PK_PhysicalExplicit : virtual public PK_Physical, public PK_Explicit<Vector> {

public:
  PK_PhysicalExplicit() {};

  PK_PhysicalExplicit(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln):
      PK_Physical(pk_tree, glist, S, soln) {};
   
  // Virtual destructor
  virtual ~PK_PhysicalExplicit() {};
};

} // namespace Amanzi

#endif
