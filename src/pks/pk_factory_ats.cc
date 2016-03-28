/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   PK factory for self-registering PKs.

   See a more thorough factory discussion in $ATS_DIR/src/factory/factory.hh.

   Simplest usage:

   // pk_implementation.hh
   #include "pk.hh"
   class DerivedPK : public PK {
     DerivedPK(Teuchos::ParameterList& plist, const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution);
     ...
   private:
     static RegisteredPKFactory<PK,DerivedPK> factory_; // my factory entry
     ...
   };

   ------------------------------------------------------------------------- */

#include "pk_factory_ats.hh"
#include "PK_Factory.hh"

namespace Amanzi {

PKFactory_ATS::map_type* PKFactory_ATS::map_ats_;
PKFactory::map_type* PKFactory::map_;

} // namespace
