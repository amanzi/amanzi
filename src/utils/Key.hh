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

// parameter list names take their parent history?
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
getKey(const Key& domain, const Key& name)
{
  return (domain.empty() || domain == std::string("domain") ) ? name : domain+"-"+name;
}

// Split a DOMAIN-VARNAME key.
inline KeyPair
splitKey(const Key& name, const char& delimiter='-')
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
  if (name.find('|') == std::string::npos) {
    // not a derivative
    return splitKey(name).first;
  } else {
    // pop the initial d
    return splitKey(name).first.substr(1,std::string::npos);
  }
}

// Grab the domain prefix of a DOMAIN-VARNAME Key, including the "-"
inline Key
getDomainPrefix(const Key& name) {
  std::size_t pos = name.find('-');
  return pos == std::string::npos ? Key("") : name.substr(0,pos+1);
}

// Grab the varname suffix of a DOMAIN-VARNAME Key
inline Key
getVarName(const Key& name) {
  return splitKey(name).second;
}

// Domain set keys are of the form DOMAIN_*-VARNAME, where * is an integer
// indexing the domain set.

// Split a domain set into DOMAIN, *, VARNAME
bool
splitDomainSet(const Key& name, KeyTriple& result);

// Check if a key, interpreted as a domain set, matches the domain-set name
inline bool
matchesDomainSet(const Key& domain_set, const Key& name) {
  KeyTriple result;
  return splitDomainSet(name, result) ? std::get<0>(result) == domain_set : false;
}

// Derivatives are of the form dKey|dKey.
inline Key
getDerivKey(Key var, Key wrt) {
  return std::string("d")+var+"|d"+wrt;
}

// Convenience function for requesting the name of a Key from an input spec.
Key
readKey(Teuchos::ParameterList& list, const Key& domain, const Key& basename,
        const Key& default_name="");

// helper functions...
inline bool starts_with(const std::string& key, const std::string& substr) {
  return key.length() >= substr.length() && key.substr(0,substr.length()) == substr;
}

inline bool ends_with(const std::string& key, const std::string& substr) {
  return key.length() >= substr.length() &&
      key.substr(key.length()-substr.length(), key.length()) == substr;
}

} // namespace Keys
} // namespace Amanzi

#endif




