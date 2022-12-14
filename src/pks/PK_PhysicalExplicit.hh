/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky
      Ethan Coon
*/

/*
  Process Kernels

  Default base with a few methods implemented in standard ways.
*/

#ifndef AMANZI_PK_PHYSICAL_EXPLICIT_HH_
#define AMANZI_PK_PHYSICAL_EXPLICIT_HH_

#include "PK_Explicit.hh"
#include "PK_Physical.hh"
#include "TreeVector.hh"

namespace Amanzi {

template <class Vector>
class PK_PhysicalExplicit : public PK_Physical, public PK_Explicit<Vector> {
 public:
  PK_PhysicalExplicit() : PK(), PK_Physical(), PK_Explicit<Vector>(){};

  PK_PhysicalExplicit(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& glist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln)
    : PK(pk_tree, glist, S, soln),
      PK_Physical(pk_tree, glist, S, soln),
      PK_Explicit<Vector>(pk_tree, glist, S, soln){};

  // Virtual destructor
  virtual ~PK_PhysicalExplicit(){};
};

} // namespace Amanzi

#endif
