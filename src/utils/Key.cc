/* 
   Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
   Amanzi is released under the three-clause BSD License. 
   The terms of use and "as is" disclaimer for this license are 
   provided in the top-level COPYRIGHT file.
   See $ATS_DIR/COPYRIGHT
   
   Author: Ethan Coon
*/


//! Keys are just strings.

/*
  Here we provide a bunch of useful string manipulation stuff that is
  specific to Amanzi style keys, domains, etc.
*/

#include "dbc.hh"
#include "Key.hh"

namespace Amanzi {
namespace Keys {

std::string
readKey(Teuchos::ParameterList& list, const Key& domain, const Key& basename, const Key& default_name) {
  std::string basename_key_arg = basename+" key";
  std::string basename_key_suffix_arg = basename+" key suffix";

  Key default_key;
  if (list.isParameter(basename_key_suffix_arg)) {
    default_key = getKey(domain, list.get<std::string>(basename_key_suffix_arg));
    return list.get<std::string>(basename_key_arg, default_key);
  } else if (!default_name.empty()) {
    default_key = getKey(domain, default_name);
    return list.get<std::string>(basename_key_arg, default_key);
  } else {
    return list.get<std::string>(basename_key_arg);
  }
}

bool
splitDomainSet(const Key& name, KeyTriple& result) {
  if (!in(name, dset_delimiter)) return false;
  Key domain;
  if (in(name, name_delimiter)) {
    auto domain_var = splitKey(name, name_delimiter);
    std::get<2>(result) = domain_var.second;
    domain = domain_var.first;
  } else {
    std::get<2>(result) = "";
    domain = name;
  }
  auto name_id = splitKey(domain, dset_delimiter);

  // if this assertion fails, the dset delimiter was in the varname!
  AMANZI_ASSERT(in(domain, dset_delimiter));
  std::get<0>(result) = name_id.first;
  std::get<1>(result) = name_id.second;
  return true;
}


} // namespace
} // namespace
