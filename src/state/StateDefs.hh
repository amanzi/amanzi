/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Some basic typedefs for State and company.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_DEFS_HH_
#define AMANZI_STATE_DEFS_HH_

#include <string>
#include <vector>
#include <set>

namespace Amanzi {

// Keys and containers
typedef std::string Key;

inline Key
getKey(Key domain, Key name) { return (domain.empty() || domain == std::string("domain") ) ? name : domain+"-"+name; }

inline Key
getDomain(Key name) {
  std::size_t pos = name.find('-');
  return pos == std::string::npos ? Key("") : name.substr(0,pos);
}

inline Key
getDomainPrefix(Key name) {
  std::size_t pos = name.find('-');
  return pos == std::string::npos ? Key("") : name.substr(0,pos+1);
}

inline Key
getDerivKey(Key var, Key wrt) {
  return std::string("d")+var+"_d"+wrt;
}

typedef std::set<Key> KeySet;
typedef std::vector<Key> KeyVector;
typedef std::set<std::pair<Key, Key> > KeyPairSet;

// Fields
typedef enum { NULL_FIELD_TYPE, COMPOSITE_VECTOR_FIELD, CONSTANT_VECTOR, CONSTANT_SCALAR } FieldType;

} // namespace

#endif
