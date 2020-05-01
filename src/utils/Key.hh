/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Keys are just strings.

/*
  Here we provide a bunch of useful string manipulation stuff that is
  specific to Amanzi style keys, domains, etc.
*/

#pragma once

#include <set>
#include "boost/algorithm/string.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

// Keys and containers
typedef std::string Key;
typedef std::set<Key> KeySet;
typedef std::vector<Key> KeyVector;

typedef std::pair<Key, Key> KeyPair;
typedef std::set<KeyPair> KeyPairSet;
typedef std::vector<KeyPair> KeyPairVector;

typedef std::tuple<Key, Key, Key> KeyTriple;
typedef std::set<KeyTriple> KeyTripleSet;

namespace Keys {

//
// A few helper string manipulation functions
//
inline bool
startsWith(const std::string& input, const std::string& start)
{
  return input.substr(0, start.size()) == start;
}

inline KeyPair
split(const std::string& input, const char& divider, const std::string& def, bool null_first)
{
  std::size_t pos = input.find(divider);
  if (pos == std::string::npos) {
    if (null_first) {
      return std::make_pair(def, input);
    } else {
      return std::make_pair(input, def);
    }
  } else {
    return std::make_pair(input.substr(0, pos), input.substr(pos+1, input.size()));
  }
}


//
// A fully resolved key is of the form  "DOMAIN-VARNAME:VARTAG"
//

// Generate a DOMAIN-VARNAME key.
inline Key
getKey(const Key& domain, const Key& name)
{
  return (domain.empty() || domain == std::string("domain")) ?
           name :
           domain + "-" + name;
}

// Split a DOMAIN-VARNAME key.
inline KeyPair
splitKey(const Key& name, const Key& def="domain")
{
  return split(name, '-', def, true);
}

// Grab the domain prefix of a DOMAIN-VARNAME key.
inline Key
getDomain(const Key& name)
{
  return splitKey(name).first;
}

// Grab the domain prefix of a DOMAIN-VARNAME Key, including the "-"
inline Key
getDomainPrefix(const Key& name)
{
  return getDomain(name)+"-";
}

// Grab the varname suffix of a DOMAIN-VARNAME Key
inline Key
getVarName(const Key& name)
{
  return splitKey(name).second;
}

// Domain set keys are of the form DOMAIN_*-VARNAME, where * is an integer
// indexing the domain set.

// Split a domain set into DOMAIN, *, VARNAME
bool
splitDomainSet(const Key& name, KeyTriple& result);

// Check if a key, interpreted as a domain set, matches the domain-set name
inline bool
matchesDomainSet(const Key& domain_set, const Key& name)
{
  KeyTriple result;
  return splitDomainSet(name, result) ? std::get<0>(result) == domain_set :
                                        false;
}


// Tag'd variables are of the form VARNAME:TAG
inline Key
getKeyTag(const Key& var, const Key& tag)
{
  return var + ":" + tag;
}

// Split a DOMAIN-VARNAME key.
inline KeyPair
splitKeyTag(const Key& name)
{
  return split(name, ':', "", false);
}


// Convenience function for requesting the name of a Key from an input spec.
Key
readKey(Teuchos::ParameterList& list, const Key& domain, const Key& basename,
        const Key& default_name="");

// Convenience function for requesting a list of names of Keys from an input
// spec.
Teuchos::Array<Key>
readKeys(Teuchos::ParameterList& list, const Key& domain, const Key& basename,
         Teuchos::Array<Key> const* const default_names=nullptr);


// Convenience function to see if a map (or map-like) object has a key.
template <typename T, typename K>
bool
hasKey(const T& container, const K& key)
{
  return container.count(key) > 0;
}

// parameter list names take their parent history?
inline Key
cleanPListName(std::string name)
{
  auto res = boost::algorithm::find_last(name, "->");
  if (res.end() - name.end() != 0)
    boost::algorithm::erase_head(name, res.end() - name.begin());
  return name;
}


} // namespace Keys
} // namespace Amanzi
