/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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

namespace Amanzi {
namespace Operators {

MatrixMFD_Factory::map_type* MatrixMFD_Factory::map_;

} // namespace
} // namespace
