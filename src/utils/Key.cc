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

// Convenience function for requesting the name of a Key from an input spec.
//
// helper functions...
bool starts_with(const Key& key, const std::string& substr)
{
  return key.length() >= substr.length() && key.substr(0,substr.length()) == substr;
}

bool ends_with(const Key& key, const std::string& substr)
{
  return key.length() >= substr.length() &&
      key.substr(key.length()-substr.length(), key.length()) == substr;
}

bool in(const Key& key, const char& c)
{
  return key.find(c) != std::string::npos;
}

Key cleanPListName(std::string name)
{
  auto res = boost::algorithm::find_last(name,"->");
  if (res.end() - name.end() != 0) boost::algorithm::erase_head(name, res.end() - name.begin());
  return name;
}

// Keys are often a combination of a domain and a variable name.

// Generate a DOMAIN-VARNAME key.
Key getKey(const Key& domain, const Key& name, const char& delimiter)
{
  return (domain.empty() || domain == std::string("domain") ) ? name : domain+delimiter+name;
}

// Split a DOMAIN-VARNAME key.
KeyPair splitKey(const Key& name, const char& delimiter)
{
  std::size_t pos = name.find(delimiter);
  if (pos == std::string::npos)
    return std::make_pair(Key(""), name);
  else
    return std::make_pair(name.substr(0,pos), name.substr(pos+1,name.size()));
}

// Grab the domain prefix of a DOMAIN-VARNAME key.
Key getDomain(const Key& name)
{
  if (name.find(deriv_delimiter) == std::string::npos) {
    // not a derivative
    return splitKey(name).first;
  } else {
    auto split = splitKey(name);
    if (split.first.size() == 0) {
      return split.first;
    } else {
      // pop the initial d
      return splitKey(name).first.substr(1,std::string::npos);
    }
  }
}

// Grab the domain prefix of a DOMAIN-VARNAME Key, including the delimiter
Key getDomainPrefix(const Key& name)
{
  std::size_t pos = name.find(name_delimiter);
  return pos == std::string::npos ? Key("") : name.substr(0,pos+1);
}

// Grab the varname suffix of a DOMAIN-VARNAME Key
Key getVarName(const Key& name)
{
  return splitKey(name).second;
}

// Domain Sets are of the form NAME:ID, where ID is an integer or
// region string indexing the domain set.
Key getDomainInSet(const Key& ds_name, const Key& subdomain)
{
  return getKey(ds_name, subdomain, dset_delimiter);
}

Key getDomainInSet(const Key& ds_name, const int& subdomain)
{
  return getKey(ds_name, std::to_string(subdomain), dset_delimiter);
}

Key getDomainSetName(const Key& name_id)
{
  return splitKey(name_id, dset_delimiter).first;
}

bool isDomainSet(const Key& name)
{
  KeyTriple result;
  return splitDomainSet(name, result);
}

// reconstruct a key from components
Key getKey(const Key& ds_name, const Key& ds_id, const Key& varname)
{
  return getKey(getKey(ds_name, ds_id, dset_delimiter), varname);
}

// reconstruct a key from components
Key getKey(const Key& ds_name, const int& ds_id, const Key& varname)
{
  return getKey(ds_name, std::to_string(ds_id), varname);
}

// Check if a key, interpreted as a domain set, matches the domain-set name
bool matchesDomainSet(const Key& domain_set, const Key& name)
{
  KeyTriple result;
  return splitDomainSet(name, result) ? std::get<0>(result) == domain_set : false;
}

// Derivatives are of the form dKey|dKey.
Key getDerivKey(const Key& var, const Key& wrt)
{
  return std::string("d")+var+deriv_delimiter+"d"+wrt;
}


std::string
readKey(Teuchos::ParameterList& list, const Key& domain, const Key& basename, const Key& default_name)
{
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


bool splitDomainSet(const Key& name, KeyTriple& result)
{
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
