/*
  Amanzi

  License: see $AMANZI_DIR/COPYRIGHT
  Author: Ethan Coon

  PK factory for self-registering PKs.

  See a more thorough factory discussion in src/utils/Factory.hh.

  Usage:

  Add a private, static member of type RegisteredPKFactory to the class
  declaration, and a special _reg.hh file that instantiates the static
  registry.

  Example:

  // pk_implementation.hh
  #include "PK.hh"
  #include "PK_Factory.hh"
  class DerivedPK : public Amanzi::PK {
    ...
   private:
    static Amanzi::RegisteredPKFactory<DerivedPK> factory_;
    ...
  };

  // pk_implementation_reg.hh
  #include "pk_implementation.hh"
  template<>
  Amanzi::RegisteredPKFactory<DerivedPK> DerivedPK::factory_("pk unique id");
  
*/

#include "PK_Factory.hh"

namespace Amanzi {

PKFactory::map_type* PKFactory::map_;

} // namespace
