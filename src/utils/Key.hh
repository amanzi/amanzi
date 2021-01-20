/* -*-  mode: c++; indent-tabs-mode: nil -*- */
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

  Currently the following characters are reserved, and should not be used in
  names:

  * '-' is used between domain names and variable names, e.g. "surface-ponded_depth"
  * '|' is used in derivative names, e.g. "dwater_content|dpressure"
  * ':' is used in domain sets, between the set name and the ID, e.g. "column:0" or "subdomain:upstream"

  Note these are distinct and can be combined, e.g. the following are valid names:
  * 'column:0-pressure'
  * 'subdomain_surface:upstream_surface-ponded_depth'
  * 'dsubdomain_surface:upstream_surface-ponded_depth|dsubdomain_surface:upstream_surface-pressure'

  This can get quite ugly, but is needed for generic code.  Note that both ' '
  and '_' are valid characters in names, but users should prefer to use
  underscores as some tools can have difficulty with spaces.

*/


#ifndef UTILS_KEY_HH_
#define UTILS_KEY_HH_


#include <set>
#include "boost/algorithm/string.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

// Keys and containers
typedef std::string Key;
typedef std::set<Key> KeySet;
typedef std::vector<Key> KeyVector;

typedef std::pair<Key,Key> KeyPair;
typedef std::set<std::pair<Key, Key> > KeyPairSet;

typedef std::tuple<Key,Key,Key> KeyTriple;

namespace Keys {

static const char name_delimiter = '-';
static const char deriv_delimiter = '|';
static const char dset_delimiter = ':';

// Convenience function for requesting the name of a Key from an input spec.
//
// helper functions...
inline bool starts_with(const Key& key, const std::string& substr) {
  return key.length() >= substr.length() && key.substr(0,substr.length()) == substr;
}

inline bool ends_with(const Key& key, const std::string& substr) {
  return key.length() >= substr.length() &&
      key.substr(key.length()-substr.length(), key.length()) == substr;
}

inline bool in(const Key& key, const char& c) {
  return key.find(c) != std::string::npos;
}

// for parsing parameter lists
Key
readKey(Teuchos::ParameterList& list, const Key& domain, const Key& basename,
        const Key& default_name="");

inline Key
cleanPListName(std::string name)
{
  auto res = boost::algorithm::find_last(name,"->");
  if (res.end() - name.end() != 0) boost::algorithm::erase_head(name, res.end() - name.begin());
  return name;
}

// Keys are often a combination of a domain and a variable name.

// Generate a DOMAIN-VARNAME key.
inline Key
getKey(const Key& domain, const Key& name, const char& delimiter=name_delimiter)
{
  return (domain.empty() || domain == std::string("domain") ) ? name : domain+delimiter+name;
}

// Split a DOMAIN-VARNAME key.
inline KeyPair
splitKey(const Key& name, const char& delimiter=name_delimiter)
{
  std::size_t pos = name.find(delimiter);
  if (pos == std::string::npos)
    return std::make_pair(Key(""), name);
  else
    return std::make_pair(name.substr(0,pos), name.substr(pos+1,name.size()));
}

// Grab the domain prefix of a DOMAIN-VARNAME key.
inline Key
getDomain(const Key& name) {
  if (name.find(deriv_delimiter) == std::string::npos) {
    // not a derivative
    return splitKey(name).first;
  } else {
    // pop the initial d
    return splitKey(name).first.substr(1,std::string::npos);
  }
}

// Grab the domain prefix of a DOMAIN-VARNAME Key, including the delimiter
inline Key
getDomainPrefix(const Key& name) {
  std::size_t pos = name.find(name_delimiter);
  return pos == std::string::npos ? Key("") : name.substr(0,pos+1);
}

// Grab the varname suffix of a DOMAIN-VARNAME Key
inline Key
getVarName(const Key& name) {
  return splitKey(name).second;
}

// Domain Sets are of the form NAME:ID, where ID is an integer or
// region string indexing the domain set.
inline Key
getDomainInSet(const Key& ds_name, const Key& subdomain)
{
  return getKey(ds_name, subdomain, dset_delimiter);
}

inline Key
getDomainInSet(const Key& ds_name, const int& subdomain)
{
  return getKey(ds_name, std::to_string(subdomain), dset_delimiter);
}

inline Key
getDomainSetName(const Key& name_id)
{
  return splitKey(name_id, dset_delimiter).first;
}

template<typename Index=std::string>
Index getDomainSetIndex(const Key& name_id)
{
  return splitKey(name_id, dset_delimiter).second;
}

template<>
inline int getDomainSetIndex(const Key& name_id)
{
  return std::atoi(getDomainSetIndex<std::string>(name_id).c_str());
}

// Split a domain set into DOMAIN, *, VARNAME
bool splitDomainSet(const Key& name, KeyTriple& result);

inline bool
isDomainSet(const Key& name) {
  KeyTriple result;
  return splitDomainSet(name, result);
}

// reconstruct a key from components
inline Key
getKey(const Key& ds_name, const Key& ds_id, const Key& varname)
{
  return getKey(getKey(ds_name, ds_id, dset_delimiter), varname);
}

// reconstruct a key from components
inline Key
getKey(const Key& ds_name, const int& ds_id, const Key& varname)
{
  return getKey(ds_name, std::to_string(ds_id), varname);
}

// Check if a key, interpreted as a domain set, matches the domain-set name
inline bool
matchesDomainSet(const Key& domain_set, const Key& name) {
  KeyTriple result;
  return splitDomainSet(name, result) ? std::get<0>(result) == domain_set : false;
}

// Derivatives are of the form dKey|dKey.
inline Key
getDerivKey(const Key& var, const Key& wrt) {
  return std::string("d")+var+deriv_delimiter+"d"+wrt;
}

} // namespace Keys
} // namespace Amanzi

#endif




