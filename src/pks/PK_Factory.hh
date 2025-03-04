/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Process kernel factory for creating PKs from the global input spec.
/*
  Developer notes:

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
    static Amanzi::RegisteredPKFactory<DerivedPK> reg_;
    ...
  };

  // pk_implementation_reg.hh
  #include "pk_implementation.hh"
  template<>
  Amanzi::RegisteredPKFactory<DerivedPK> DerivedPK::reg_("pk unique id");


  ETC: this was somewhat a poor interface decision.  The pk_tree list should
  likely be passed by shared pointer instead of by reference, as that is
  already the internal storage within a PList.
*/

#ifndef AMANZI_PK_FACTORY_HH_
#define AMANZI_PK_FACTORY_HH_

#include <iostream>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "StateDefs.hh"
#include "PK.hh"


namespace Amanzi {

class TreeVector;
class State;

class PKFactory {
 public:
  // Method to create a PK
  Teuchos::RCP<PK> CreatePK(std::string pk_name,
                            Teuchos::ParameterList& pk_tree,
                            const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                            const Teuchos::RCP<State>& state,
                            const Teuchos::RCP<TreeVector>& soln);

  // typedef describing a map with string keys and PK constructor values, this
  // stores constructors for all known PK classes
  typedef std::map<std::string,
                   PK* (*)(Teuchos::ParameterList&,
                           const Teuchos::RCP<Teuchos::ParameterList>&,
                           const Teuchos::RCP<State>&,
                           const Teuchos::RCP<TreeVector>&)>
    map_type;

 protected:
  static map_type* GetMap()
  {
    if (!map_) map_ = new map_type;
    return map_;
  }

 public:
  static int num_pks;
  static std::string list_pks;

 private:
  static map_type* map_;
};


template <typename T>
PK*
CreateT(Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
        const Teuchos::RCP<State>& state,
        const Teuchos::RCP<TreeVector>& soln)
{
  return new T(pk_tree, global_list, state, soln);
}


template <typename T>
class RegisteredPKFactory : public PKFactory {
 public:
  // Constructor for the registered factory.  Needs some error checking in
  // case a name s is already in the map? (i.e. two implementations trying to
  // call themselves the same thing) --etc
  RegisteredPKFactory(const std::string& s)
  {
    GetMap()->insert(std::pair<std::string,
                               PK* (*)(Teuchos::ParameterList&,
                                       const Teuchos::RCP<Teuchos::ParameterList>&,
                                       const Teuchos::RCP<State>&,
                                       const Teuchos::RCP<TreeVector>&)>(s, &CreateT<T>));
  }
};

} // namespace Amanzi

#endif
