/*
  This is the PKs component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

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

#ifndef AMANZI_PK_FACTORY_HH_
#define AMANZI_PK_FACTORY_HH_

#include <iostream>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "PK.hh"

namespace Amanzi {

class TreeVector;
class State;

class PKFactory {
 public:
  typedef std::map<std::string,
                   PK* (*)(Teuchos::ParameterList&,
                           const Teuchos::RCP<Teuchos::ParameterList>&,
                           const Teuchos::RCP<State>&,
                           const Teuchos::RCP<TreeVector>&)> map_type;

  static Teuchos::RCP<PK> CreatePK(Teuchos::ParameterList& pk_tree,
          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
          const Teuchos::RCP<State>& state,
          const Teuchos::RCP<TreeVector>& soln) {

    if (!global_list->isSublist("PKs")) {
      Errors::Message message("PK_Factory: Missing sublist \"PKs\" in global list.");
      Exceptions::amanzi_throw(message);
    }
   
    std::string pk_type;
    if (pk_tree.isParameter("PK type")){
      pk_type = pk_tree.get<std::string>("PK type");
      // if (!global_list->sublist("PKs").isSublist(pk_tree.name(pk_ite)) {
    }
    else{
      std::stringstream errmsg;
      errmsg << "PK_Factory: Missing PK type item \n";
      Errors::Message message(errmsg.str());
      Exceptions::amanzi_throw(message);
    }

    //std::string s = global_list->sublist("PKs").sublist(pk_name).get<std::string>("PK type");

    map_type::iterator iter = GetMap()->find(pk_type);
    if (iter == GetMap()->end()) {
      std::stringstream errmsg;
      errmsg << "PK Factory: cannot find PK type \"" << pk_type << "\"";
      for (map_type::iterator iter=GetMap()->begin();
           iter!=GetMap()->end(); ++iter) {
        errmsg << std::endl << "  option: " << iter->first;
      }
      Errors::Message message(errmsg.str());
      Exceptions::amanzi_throw(message);
    }
    return Teuchos::rcp(iter->second(pk_tree, global_list, state, soln));
  }

  // static Teuchos::RCP<PK> CreatePK(const Teuchos::RCP<Teuchos::ParameterList>& plist,
  //         Teuchos::ParameterList& FElist,
  //         const Teuchos::RCP<TreeVector>& soln) {

  //   std::string s = plist->get<std::string>("PK type");
  //   map_type::iterator iter = GetMap()->find(s);
  //   if (iter == GetMap()->end()) {
  //     std::stringstream errmsg;
  //     errmsg << "PK Factory: cannot find PK type \"" << s << "\"";
  //     for (map_type::iterator iter=GetMap()->begin();
  //          iter!=GetMap()->end(); ++iter) {
  //       errmsg << std::endl << "  option: " << iter->first;
  //     }
  //     Errors::Message message(errmsg.str());
  //     Exceptions::amanzi_throw(message);
  //   }
  //   return Teuchos::rcp(iter->second(plist, FElist, soln));
  // }


 private:
  static map_type* map_;

 protected:
  static map_type* GetMap() {
    if (!map_) map_ = new map_type;
    return map_;
  }
};


template<typename T> PK* CreateT(Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
        const Teuchos::RCP<State>& state,
        const Teuchos::RCP<TreeVector>& soln) {
  return new T(pk_tree, global_list, state, soln);
}



template<typename T>
class RegisteredPKFactory : public PKFactory {
public:
  // Constructor for the registered factory.  Needs some error checking in
  // case a name s is already in the map? (i.e. two implementations trying to
  // call themselves the same thing) --etc
  RegisteredPKFactory(const std::string& s) {
    GetMap()->insert(std::pair<std::string, PK* (*)(Teuchos::ParameterList&,
                     const Teuchos::RCP<Teuchos::ParameterList>&,
                     const Teuchos::RCP<State>&,
                     const Teuchos::RCP<TreeVector>&)>(s, &CreateT<T>));
  }
};

}  // namespace Amanzi

#endif
