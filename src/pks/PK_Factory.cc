/*
  Process Kernels

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

#include "PK_Factory.hh"
#include "State.hh"

namespace Amanzi {

Teuchos::RCP<PK>
PKFactory::CreatePK(std::string pk_name,
                     Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                     const Teuchos::RCP<State>& state,
                     const Teuchos::RCP<TreeVector>& soln)
{
  // make sure we can find PKs
  if (!global_list->isSublist("PKs")) {
    Errors::Message message("PK_Factory: Missing sublist \"PKs\" in global list.");
    Exceptions::amanzi_throw(message);
  }

  // This is the constructed PK's subtree of the full PK tree
  Teuchos::ParameterList pk_subtree;
  bool pk_subtree_found = false;
    
  if (pk_tree.isSublist(pk_name)) {
    pk_subtree = pk_tree.sublist(pk_name);
    pk_subtree_found = true;

  } else {
    // check for "flyweight" PKs, or PKs which share a common input spec with
    // other PKs, differing only in domain (and the domains of their children
    // in the case of MPCs)
    KeyTriple pk_triple;
    bool is_ds = Keys::splitDomainSet(pk_name, pk_triple);
    if (is_ds) {
      // flyweight PKs are defined across a domain_set
      Teuchos::Array<std::string> domain_sets;    
      domain_sets = global_list->sublist("state")
                    .get<Teuchos::Array<std::string> >("domain sets", domain_sets);
      for (auto ds : domain_sets) {
        if (ds == std::get<0>(pk_triple)) {
          // flyweight PK, alter the sublist and construct
          // -- get the domain name and base varname
          Key pk_flyweight = Keys::getKey(ds+"_*", std::get<2>(pk_triple));
          Teuchos::ParameterList pk_list_new = global_list->sublist("PKs").sublist(pk_flyweight);

          // -- overwrite the domain name
          Key new_domain = ds+"_"+std::get<1>(pk_triple);
          if (pk_list_new.isParameter("domain name"))
            pk_list_new.set("domain name", new_domain);

          // -- overwrite sub pks names with prepended domain
          if (pk_list_new.isParameter("PKs order")) {
            auto subpks_names = pk_list_new.get<Teuchos::Array<std::string> >("PKs order");
            for (auto& subpk_name : subpks_names) {
              KeyTriple subpk_triple;
              bool subpk_is_ds = Keys::splitDomainSet(subpk_name, subpk_triple);

              if (subpk_is_ds) {
                subpk_name = Keys::getKey(std::get<0>(subpk_triple)+"_"+std::get<1>(pk_triple),
                        std::get<2>(subpk_triple));
              }
            }
            pk_list_new.set("PKs order", subpks_names);
          }
            
          // push into the PKs list
          global_list->sublist("PKs").set(pk_name, pk_list_new);

          // get the flyweight subtree
          if (!pk_tree.isSublist(pk_flyweight)) {
            std::stringstream msg;
            msg << "PK_Factory: PK \"" << pk_name << "\" is a flyweight, but missing flyweight PK spec \""
                << pk_flyweight << "\"\n";
            Errors::Message message(msg.str());
            Exceptions::amanzi_throw(message);
          }            
          pk_subtree = pk_tree.sublist(pk_flyweight);
          pk_subtree.setName(pk_name);

          pk_subtree_found = true;
        }
      }
    }
  }

  // ensure we found the pk subtree
  if (!pk_subtree_found) {
    std::stringstream msg;
    msg << "PK_Factory: PK \"" << pk_name << "\" not a sublist of the provided PK tree:\n" << pk_tree;
    Errors::Message message(msg.str());
    Exceptions::amanzi_throw(message);
  }
    
  // get the PK type
  std::string pk_type;
  if (pk_subtree.isParameter("PK type")) {
    pk_type = pk_subtree.get<std::string>("PK type");
  } else {
    std::stringstream msg;
    msg << "PK_Factory: PK \"" << pk_name << "\" is missing a PK type:\n" << pk_subtree;
    Errors::Message message(msg.str());
    Exceptions::amanzi_throw(message);
  }

  // find the pk type
  map_type::iterator iter = GetMap()->find(pk_type);
  if (iter == GetMap()->end()) {
    std::stringstream message;
    message << "PK Factory: PK \"" << pk_name << "\" requested type \""
            << pk_type << "\" which is not a registered PK type.\n";
    
    for (map_type::iterator it = GetMap()->begin(); it != GetMap()->end(); ++it) {
      message  << std::endl << "  option: " << it->first;
    }
    Errors::Message msg(message.str());
    Exceptions::amanzi_throw(msg);
  }

  // construct the PK
  num_pks++;
  if (list_pks.size() < 1024) list_pks += "|" + pk_name;
  return Teuchos::rcp(iter->second(pk_subtree, global_list, state, soln));
}


std::string PKFactory::list_pks;
int PKFactory::num_pks = 0;

PKFactory::map_type* PKFactory::map_;

}  // namespace Amanzi
