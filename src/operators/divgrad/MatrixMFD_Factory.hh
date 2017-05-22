/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   MatrixMFD factory.

   See a more thorough factory discussion in $ATS_DIR/src/factory/factory.hh.

   Simplest usage:

   // pk_implementation.hh
   #include "pk.hh"
   class DerivedMatrixMFD : public MatrixMFD {
     DerivedMatrixMFD(Teuchos::ParameterList& plist,
               const Teuchos::RCP<TreeVector>& solution);
     ...
   private:
     static RegisteredMatrixMFDFactory<MatrixMFD,DerivedMatrixMFD> factory_; // my factory entry
     ...
   };

   ------------------------------------------------------------------------- */

#ifndef ATS_MatrixMFD_FACTORY_HH_
#define ATS_MatrixMFD_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Mesh.hh"

namespace Amanzi {
namespace Operators {

class MatrixMFD;
class MatrixMFD_Coupled;

Teuchos::RCP<MatrixMFD>
CreateMatrixMFD(Teuchos::ParameterList& plist,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

Teuchos::RCP<MatrixMFD_Coupled>
CreateMatrixMFD_Coupled(Teuchos::ParameterList& plist,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);


} // namespace
} // namespace

#endif
