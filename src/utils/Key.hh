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

#pragma once

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
bool starts_with(const Key& key, const std::string& substr);
bool ends_with(const Key& key, const std::string& substr);
bool in(const Key& key, const char& c);

// for parsing parameter lists
Key readKey(Teuchos::ParameterList& list, const Key& domain, const Key& basename,
        const Key& default_name="");

Key cleanPListName(std::string name);

// Keys are often a combination of a domain and a variable name.

// Generate a DOMAIN-VARNAME key.
Key getKey(const Key& domain, const Key& name, const char& delimiter=name_delimiter);

// Split a DOMAIN-VARNAME key.
KeyPair splitKey(const Key& name, const char& delimiter=name_delimiter);

// Grab the domain prefix of a DOMAIN-VARNAME key.
Key getDomain(const Key& name);

// Grab the domain prefix of a DOMAIN-VARNAME Key, including the delimiter
Key getDomainPrefix(const Key& name);

// Grab the varname suffix of a DOMAIN-VARNAME Key
Key getVarName(const Key& name);

// Domain Sets are of the form NAME:ID, where ID is an integer or
// region string indexing the domain set.
Key getDomainInSet(const Key& ds_name, const Key& subdomain);

Key getDomainInSet(const Key& ds_name, const int& subdomain);

Key getDomainSetName(const Key& name_id);

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

bool isDomainSet(const Key& name);

// reconstruct a key from components
Key getKey(const Key& ds_name, const Key& ds_id, const Key& varname);

// reconstruct a key from components
Key getKey(const Key& ds_name, const int& ds_id, const Key& varname);

// Check if a key, interpreted as a domain set, matches the domain-set name
bool matchesDomainSet(const Key& domain_set, const Key& name);

// Derivatives are of the form dKey|dKey.
Key getDerivKey(const Key& var, const Key& wrt);


} // namespace Keys
} // namespace Amanzi
