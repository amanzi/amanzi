/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
   Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
   Amanzi is released under the three-clause BSD License.
   The terms of use and "as is" disclaimer for this license are
   provided in the top-level COPYRIGHT file.

   Author: Ethan Coon
*/
//! Keys are just strings.

/*

  Yet another string library...

  Here we provide a bunch of useful string manipulation stuff that is
  specific to Amanzi style keys, domains, etc.

  DOMAINs may be the default domain (empty string) which is synonymous for
  "domain"; a standard domain, e.g. "surface"; or a domain set, which is
  characterized by a domain set name and a domain index DSET:INDEX,
  e.g. "watershed:upstream" or "column:3".  Note that the subsurface domain is
  called "domain" and is typically considered the "primary" or "default"
  domain.  This domain is special in that most Keys on that domain omit the
  domain, e.g. "pressure" instead of "domain-pressure".

  VARIABLE names are often long and descriptive,
  e.g. "incoming_shortwave_radiation" and should avoid spaces.

  A KEY is a combination of a DOMAIN name and a VARIABLE name, DOMAIN-VARIABLE,
  e.g. "surface-temperature".

  A KEY may also be a derivative, which is given by the KEY being derived,
  followed by the KEY differentiated with respect to,
  e.g. "dsurface-energy|dsurface-temperature".  Note that the "d" prefix is for
  improved readibility and will be removed in all operations on derivatives.

  Currently the following characters are reserved, and should not be used in
  names:

  * '-' is used between domain names and variable names, e.g. "surface-ponded_depth"
  * ':' is used in domain sets, between the DSET name and the INDEX,
     e.g. "column:0" or "subdomain:upstream"
  * '|' is used in derivative names, e.g. "dwater_content|dpressure"
  * '.' is used in a fully qualified name, such as for components of a
     MultiVector or CompositeVector, e.g. saturation.cell.0 or
     free_ion_concentration.cell.NO3-
  * ' ' should be avoided as it causes confusion, but names with spaces are
     still valid and do not break the code (though they may break visualization
     and other scripts).  Prefer to use '_' instead, which may appear in any
     part of a name.

  All of the following are valid KEYs:

  * "column:0-pressure" : DOMAIN="column:0", VARIABLE="pressure"
  * "subdomain_surface:upstream_surface-ponded_depth" : DOMAIN="subdomain_surface:upstream_surface", VARIABLE="ponded_depth"
  * "dsubdomain_surface:upstream_surface-ponded_depth|dsubdomain_surface:upstream_surface-pressure"
     : the derivative of the above quantity with respect to pressure.

  This can get quite ugly, but is needed for generic code.

*/

#pragma once

#include <set>
#include "Teuchos_ParameterList.hpp"
#include "Tag.hh"
#include "FIFO_Set.hh"

namespace Amanzi {

// Keys and containers
typedef std::string Key;
typedef std::set<Key> KeySet;
typedef std::vector<Key> KeyVector;

typedef std::pair<Key,Key> KeyPair;
typedef std::set<std::pair<Key, Key> > KeyPairSet;
typedef std::vector<KeyPair> KeyPairVector;

typedef std::tuple<Key,Key,Key> KeyTriple;
typedef std::set<KeyTriple> KeyTripleSet;

typedef std::pair<Key,Tag> KeyTag;
typedef std::vector<KeyTag> KeyTagVector;
using KeyTagSet = Utils::FIFO_Set<KeyTag>;

typedef std::tuple<Key,Tag,Key> DerivativeTriple;
typedef std::set<DerivativeTriple> DerivativeTripleSet;


namespace Keys {

static const char name_delimiter = '-';
static const char deriv_delimiter = '|';
static const char dset_delimiter = ':';
static const char tag_delimiter = '@';
static const char phase_delimiter = '_';

//
// Utility functions
// -----------------------------------------------------------------------------
bool starts_with(const Key& key, const char& c);
bool starts_with(const Key& key, const std::string& substr);
bool ends_with(const Key& key, const char& c);
bool ends_with(const Key& key, const std::string& substr);
bool in(const Key& key, const char& c);
bool in(const Key& key, const std::string& c);

Key merge(const Key& domain, const Key& name, const char& delimiter);
KeyPair split(const Key& name, const char& delimiter);

template<class Container, class Key>
bool hasKey(const Container& c, const Key& key) {
  return (bool) c.count(key);
}

//
// Working with Keys
// -----------------------------------------------------------------------------
inline
Key standardize(const Key& other) {
  if (other.empty()) return "domain";
  if (other == "subsurface") return "domain";
  return other;
}

// is this valid?
bool validKey(const Key& key);

// Generate a DOMAIN-VARNAME key.
// Note, if DOMAIN == "domain" or "" this returns VARNAME
Key getKey(const Key& domain, const Key& name);

// Split a DOMAIN-VARNAME key.
//
// May return an empty domain string
KeyPair splitKey(const Key& name);

// Get the domain prefix of a DOMAIN-VARNAME key.  Does not return empty string.
Key getDomain(const Key& name);
Key getDomainPrefix(const Key& name);

// Grab the varname suffix of a DOMAIN-VARNAME Key
Key getVarName(const Key& name);

//
// Domain Sets
// -----------------------------------------------------------------------------
// Domain Sets are of the form NAME:ID, where ID is an integer or
// region string indexing the domain set.
Key getDomainInSet(const Key& ds_name, const Key& subdomain);

Key getDomainInSet(const Key& ds_name, const int& subdomain);

Key getDomainSetName(const Key& name_id);

template<typename Index=std::string>
Index getDomainSetIndex(const Key& name_id)
{
  return split(name_id, dset_delimiter).second;
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


//
// Working with tags
// -----------------------------------------------------------------------------
// Tag'd variables are of the form VARNAME:TAG
Key getKey(const Key& var, const Tag& tag);
Key getKey(const KeyTag& var_tag);

KeyTag splitKeyTag(const Key& name);

//
// Working with derivatives
// -----------------------------------------------------------------------------
// Derivatives are of the form dKey|dKey.
Key getDerivKey(const Key& var, const Key& wrt);
Key getDerivKey(const Key& var, const Tag& tag,
                const Key& wrt, const Tag& wrt_tag);
Key getDerivKey(const KeyTag& var, const KeyTag& wrt);

//
// Helper functions for reading keys and domains from parameter lists
// -----------------------------------------------------------------------------

// Trilinos ParameterList includes the full tree of the name
// (main->sublist->name), clean this and just return the final name.
Key cleanPListName(const std::string& name);

inline
Key cleanPListName(const Teuchos::ParameterList& plist) {
  return cleanPListName(plist.name());
}

// Read a domain name in a standard way, potentially with a dtype
Key readDomain(Teuchos::ParameterList& plist, const Key& dtype="");

// Read a domain name in a standard way, potentially with a dtype
Key readDomain(Teuchos::ParameterList& plist,
               const Key& dtype,
               const Key& default_name);

// Read a domain name, guessing default values given a related domain as a hint.
//
// For instance, given a hint of "surface_column:0", hint_dtype="surface" and
// dtype="snow", this will correctly guess "snow_column:0" as the default
// value.
Key readDomainHint(Teuchos::ParameterList& plist,
                   const Key& hint,
                   Key hint_dtype,
                   Key dtype);

// Most domains are of a common type, either "domain", "surface", "snow",
// "canopy", etc.  This tries to guess the type for use in above read calls.
Key guessDomainType(const Key& domain);


// Read a suffix to construct domain sets.
Key readSuffix(Teuchos::ParameterList& list,
               const Key& basename,
               const Key& default_name="");


// Read a Key from a parameter list, given some default information Here
// "domain" is the supposed domain of the Key, "basename" is the
// "user-friendly, readable" version of the name that will be searched for in
// the parameter list (e.g. "molar density", "canopy water content", etc, and
// "default_name" is the default value of the VARIABLE name (not the full Key!)
Key readKey(Teuchos::ParameterList& list,
            const Key& domain,
            const Key& basename,
            const Key& default_name="");

// Convenience function for requesting a list of names of Keys from an input
// spec.
Teuchos::Array<Key>
readKeys(Teuchos::ParameterList& list, const Key& domain, const Key& basename,
         Teuchos::Array<Key> const* const default_names=nullptr);


} // namespace Keys
} // namespace Amanzi
