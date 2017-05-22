/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   MatrixMFD factory for self-registering MatrixMFDs.

   See a more thorough factory discussion in $ATS_DIR/src/factory/factory.hh.

   Simplest usage:

   // pk_implementation.hh
   #include "pk.hh"
   class DerivedMatrixMFD : public MatrixMFD {
     DerivedMatrixMFD(Teuchos::ParameterList& plist, const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution);
     ...
   private:
     static RegisteredMatrixMFDFactory<MatrixMFD,DerivedMatrixMFD> factory_; // my factory entry
     ...
   };

   ------------------------------------------------------------------------- */

#include "MatrixMFD_Factory.hh"
#include "MatrixMFD_Surf_ScaledConstraint.hh"
#include "MatrixMFD_Surf.hh"
#include "MatrixMFD_TPFA_ScaledConstraint.hh"
#include "MatrixMFD_TPFA.hh"
#include "Matrix_TPFA.hh"
#include "Matrix_TPFA_Surf.hh"
#include "MatrixMFD_ScaledConstraint.hh"
#include "MatrixMFD.hh"
#include "MatrixMFD_Coupled_Surf.hh"
#include "MatrixMFD_Coupled_TPFA.hh"
#include "MatrixMFD_Coupled.hh"

namespace Amanzi {
namespace Operators {

Teuchos::RCP<MatrixMFD>
CreateMatrixMFD(Teuchos::ParameterList& plist,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {

  bool surf = plist.get<bool>("coupled to surface", false);
  bool tpfa = plist.get<bool>("TPFA", false);
  bool fv = plist.get<bool>("FV", false);
  bool scaled_constraint = plist.get<bool>("scaled constraint equation", false);

  if (surf && tpfa) {
    Errors::Message msg("MatrixMFD Factory: both Surf and TPFA cannot be specified.");
    Exceptions::amanzi_throw(msg);
  }

  if (surf) {
    if (scaled_constraint) {
      return Teuchos::rcp(new MatrixMFD_Surf_ScaledConstraint(plist, mesh));
    } else {
      return Teuchos::rcp(new MatrixMFD_Surf(plist, mesh));
    }
  } else if (tpfa) {
    if (scaled_constraint) {
      return Teuchos::rcp(new MatrixMFD_TPFA_ScaledConstraint(plist, mesh));
    } else {
      return Teuchos::rcp(new MatrixMFD_TPFA(plist, mesh));     
    }
  } else if (fv){
    if (surf){
      return Teuchos::rcp(new Matrix_TPFA_Surf(plist, mesh));  
    }
    else {
      return Teuchos::rcp(new Matrix_TPFA(plist, mesh));  
    }
  } else {
    if (scaled_constraint) {
      return Teuchos::rcp(new MatrixMFD_ScaledConstraint(plist, mesh));
    } else {
      return Teuchos::rcp(new MatrixMFD(plist, mesh));
    }
  }

  ASSERT(0);
  return Teuchos::null;
}


Teuchos::RCP<MatrixMFD_Coupled>
CreateMatrixMFD_Coupled(Teuchos::ParameterList& plist,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
  bool surf = plist.get<bool>("coupled to surface", false);
  bool tpfa = plist.get<bool>("TPFA", false);

  if (surf && tpfa) {
    Errors::Message msg("MatrixMFD_Coupled Factory: both Surf and TPFA cannot be specified.");
    Exceptions::amanzi_throw(msg);
  }

  if (surf) {
    return Teuchos::rcp(new MatrixMFD_Coupled_Surf(plist, mesh));
  } else if (tpfa) {
    return Teuchos::rcp(new MatrixMFD_Coupled_TPFA(plist, mesh));
  } else {
    return Teuchos::rcp(new MatrixMFD_Coupled(plist, mesh));
  }

  ASSERT(0);
  return Teuchos::null;
}

} // namespace
} // namespace
