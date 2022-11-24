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
#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi {
namespace Keys {

//
// Utility functions
// -----------------------------------------------------------------------------
bool
starts_with(const Key& key, const char& c)
{
  return key.length() >= 1 && c == key.front();
}
bool
starts_with(const Key& key, const std::string& substr)
{
  return key.length() >= substr.length() && key.substr(0, substr.length()) == substr;
}

bool
ends_with(const Key& key, const char& c)
{
  return key.length() >= 1 && c == key.back();
}
bool
ends_with(const Key& key, const std::string& substr)
{
  return key.length() >= substr.length() &&
         key.substr(key.length() - substr.length(), key.length()) == substr;
}

bool
in(const Key& key, const char& c)
{
  return key.find(c) != std::string::npos;
}

bool
in(const Key& key, const std::string& substr)
{
  return key.find(substr) != std::string::npos;
}

Key
replace_all(Key key, const std::string& find_s, const std::string& replace_s)
{
  std::size_t index = 0;
  while (true) {
    index = key.find(find_s, index);
    if (index == std::string::npos) break;

    key.replace(index, find_s.size(), replace_s);

    // Advance index forward so the next iteration doesn't pick it up as well.
    // This avoids problems with e.g. abc --> abcabc
    index += replace_s.size();
  }
  return key;
}

Key
merge(const Key& domain, const Key& name, const char& delimiter)
{
  return domain + delimiter + name;
}

KeyPair
split(const Key& name, const char& delimiter)
{
  std::size_t pos = name.find(delimiter);
  if (pos == std::string::npos)
    return std::make_pair(Key(""), name);
  else
    return std::make_pair(name.substr(0, pos), name.substr(pos + 1, std::string::npos));
}

// creates a clean name to be used as a variable name, domain name, tag name,
// etc, that has no delimiters in it, no spaces (which make for uglier IO),
// etc.
Key
cleanName(const std::string& name, bool delimiters_ok)
{
  Key result(name);
  if (!delimiters_ok) {
    std::replace(result.begin(), result.end(), name_delimiter, '_');
    std::replace(result.begin(), result.end(), dset_delimiter, '_');
  }
  std::replace(result.begin(), result.end(), deriv_delimiter, '_');
  std::replace(result.begin(), result.end(), tag_delimiter, '_');
  std::replace(result.begin(), result.end(), ' ', '_');
  return result;
}


//
// Working with DOMAINs
// -----------------------------------------------------------------------------
// Generate a DOMAIN-VARNAME key.
// Note, if DOMAIN == "domain" or "" this returns VARNAME

bool
validKey(const Key& key)
{
  bool result = true;
  int deriv_del_count = std::count(key.begin(), key.end(), deriv_delimiter);
  if (deriv_del_count > 1) {
    result = false;
  } else if (deriv_del_count == 1) {
    auto split_deriv = split(key, deriv_delimiter);
    result = validKey(split_deriv.first) && validKey(split_deriv.second);
  } else {
    int name_count = std::count(key.begin(), key.end(), name_delimiter);
    if (name_count > 1) {
      result = false;
    } else {
      auto domain_var = split(key, name_delimiter);
      if (starts_with(domain_var.first, dset_delimiter) ||
          ends_with(domain_var.first, dset_delimiter))
        result = false;
      if (std::count(domain_var.first.begin(), domain_var.first.end(), dset_delimiter) > 1)
        result = false;
      if (std::count(domain_var.second.begin(), domain_var.second.end(), dset_delimiter))
        result = false;
    }
  }
  return result;
}

Key
getKey(const Key& domain, const Key& variable)
{
  if (variable.empty()) {
    Errors::Message msg("Keys::getKey() called with empty variable name.");
    Exceptions::amanzi_throw(msg);
  }
  if (in(domain, name_delimiter)) {
    Errors::Message msg("Keys::getKey() called with invalid domain name \"");
    msg << domain << "\"";
    Exceptions::amanzi_throw(msg);
  }
  // there are times this is valid, e.g if the merged thing is a derivative
  // that has the delimiter in the wrt variable.
  // if (in(variable, name_delimiter)) {
  //   Errors::Message msg("Keys::getKey() called with invalid variable \"");
  //   msg << variable << "\"";
  //   Exceptions::amanzi_throw(msg);
  // }
  if (standardize(domain) == "domain") return variable;
  return merge(domain, variable, name_delimiter);
}


KeyPair
splitKey(const Key& name)
{
  if (!validKey(name)) {
    Errors::Message msg("Keys::splitKey() called with invalid argument \"");
    msg << name << "\"";
    Exceptions::amanzi_throw(msg);
  }
  return split(name, name_delimiter);
}


// Split a DOMAIN-VARNAME key.
//
// If delimiter does not appear, the first string is empty.

// Grab the domain prefix of a DOMAIN-VARNAME key.
Key
getDomain(const Key& name)
{
  if (name.find(tag_delimiter) != std::string::npos) {
    // includes tag, split the tag
    KeyTag key_tag = splitKeyTag(name);
    Key domain_in_tag = getDomain(key_tag.second.get());
    if (domain_in_tag != "domain") return domain_in_tag;
    return getDomain(key_tag.first);
  }

  if (name.find(deriv_delimiter) == std::string::npos) {
    // not a derivative
    auto split = splitKey(name);
    if (split.first.empty())
      return "domain";
    else
      return split.first;

  } else {
    // is a derivative
    auto split = splitKey(name);
    if (split.first.size() == 0) {
      return "domain";
    } else {
      // pop the initial d
      return split.first.substr(1, std::string::npos);
    }
  }
}

// Gets a prefix that can simply be "+" with the varname.  This is deprecated.
Key
getDomainPrefix(const Key& name)
{
  auto domain = getDomain(name);
  if (domain == "domain")
    return "";
  else
    return merge(domain, "", name_delimiter);
}


// Grab the varname suffix of a DOMAIN-VARNAME Key
Key
getVarName(const Key& name)
{
  if (name.find(deriv_delimiter) == std::string::npos)
    return splitKey(name).second;
  else
    return splitKey(split(name, deriv_delimiter).first).second;
}

// abbreviate
Key
abbreviate(Key name, int max_len)
{
  for (const auto& abbvs : abbreviations) {
    name = replace_all(name, abbvs.first, abbvs.second);
    if (max_len > 0 && name.size() < max_len) break;
  }
  return name;
}


// Domain Sets are of the form NAME:ID, where ID is an integer or
// region string indexing the domain set.
Key
getDomainInSet(const Key& ds_name, const Key& subdomain)
{
  if (ds_name.empty()) {
    Errors::Message msg("Cannot Keys::getDomainInSet() with empty domain set name.");
    Exceptions::amanzi_throw(msg);
  }
  if (subdomain.empty()) {
    Errors::Message msg("Cannot Keys::getDomainInSet() with empty domain set index.");
    Exceptions::amanzi_throw(msg);
  }
  return merge(ds_name, subdomain, dset_delimiter);
}

Key
getDomainInSet(const Key& ds_name, const int& subdomain)
{
  return getDomainInSet(ds_name, std::to_string(subdomain));
}

Key
getDomainSetName(const Key& name_id)
{
  if (!in(name_id, dset_delimiter)) {
    Errors::Message msg("Keys::getDomainSetName() argument \"");
    msg << name_id << "\" is not from a domain set.";
    Exceptions::amanzi_throw(msg);
  }
  return split(name_id, dset_delimiter).first;
}

bool
isDomainSet(const Key& name)
{
  KeyTriple result;
  return splitDomainSet(name, result);
}

bool
splitDomainSet(const Key& name, KeyTriple& result)
{
  if (!in(name, dset_delimiter)) return false;
  Key domain;
  if (in(name, name_delimiter)) {
    auto domain_var = splitKey(name);
    std::get<2>(result) = domain_var.second;
    domain = domain_var.first;
  } else {
    std::get<2>(result) = "";
    domain = name;
  }
  auto name_id = split(domain, dset_delimiter);
  if (name_id.first.empty() || name_id.second.empty()) {
    Errors::Message msg("Keys::splitDomainSet() argument \"");
    msg << name << "\" is not a valid name.";
    Exceptions::amanzi_throw(msg);
  }
  std::get<0>(result) = name_id.first;
  std::get<1>(result) = name_id.second;
  return true;
}

bool
isDomainInSet(const Key& name)
{
  KeyTriple result;
  return splitDomainSet(name, result);
}

// reconstruct a key from components
Key
getKey(const Key& ds_name, const Key& ds_id, const Key& varname)
{
  return getKey(merge(ds_name, ds_id, dset_delimiter), varname);
}

// reconstruct a key from components
Key
getKey(const Key& ds_name, const int& ds_id, const Key& varname)
{
  return getKey(ds_name, std::to_string(ds_id), varname);
}

// Check if a key, interpreted as a domain set, matches the domain-set name
bool
matchesDomainSet(const Key& domain_set, const Key& name)
{
  KeyTriple result;
  return splitDomainSet(name, result) ? std::get<0>(result) == domain_set : false;
}

// tags
Key
getKey(const Key& var, const Tag& tag)
{
  if (tag.get() == "")
    return var;
  else
    return merge(var, tag.get(), tag_delimiter);
}

Key
getKey(const KeyTag& var_tag)
{
  return getKey(var_tag.first, var_tag.second);
}

KeyTag
splitKeyTag(const Key& name)
{
  std::size_t pos = name.find(tag_delimiter);
  if (pos == std::string::npos)
    return std::make_pair(name, Tag(""));
  else
    return std::make_pair(name.substr(0, pos), Tag(name.substr(pos + 1, name.size())));
}

// Derivatives are of the form dKey|dKey.
Key
getDerivKey(const Key& var, const Key& wrt)
{
  std::string d("d");
  return merge(d + var, d + wrt, deriv_delimiter);
}

Key
getDerivKey(const Key& var, const Tag& tag, const Key& wrt, const Tag& wrt_tag)
{
  return getDerivKey(Keys::getKey(var, tag), Keys::getKey(wrt, wrt_tag));
}

Key
getDerivKey(const KeyTag& var, const KeyTag& wrt)
{
  return getDerivKey(Keys::getKey(var), Keys::getKey(wrt));
}


//
// Helper functions for reading keys and domains from parameter lists
// -----------------------------------------------------------------------------
Key
cleanPListName(const std::string& name)
{
  auto pos = name.rfind("->");
  if (pos == name.size()) {
    return "";
  } else if (pos == std::string::npos) {
    return name;
  } else {
    return name.substr(pos + 2, std::string::npos);
  }
}


// Read a domain name in a standard way, potentially with a prefix
Key
readDomain(Teuchos::ParameterList& plist, const Key& prefix)
{
  if (prefix.empty() || prefix == "domain") {
    return plist.get<std::string>("domain name");
  } else {
    return plist.get<std::string>(prefix + " domain name");
  }
}

// Read a domain name in a standard way, potentially with a prefix
Key
readDomain(Teuchos::ParameterList& plist, const Key& prefix, const Key& default_name)
{
  if (prefix.empty() || prefix == "domain") {
    return plist.get<std::string>("domain name", default_name);
  } else {
    return plist.get<std::string>(prefix + " domain name", default_name);
  }
}

Key
readDomainHint(Teuchos::ParameterList& plist,
               const Key& hint,
               Key hint_domain_type,
               Key domain_type)
{
  std::string param;
  if (domain_type.empty())
    param = "domain name";
  else
    param = (domain_type + " domain name");

  domain_type = standardize(domain_type);
  hint_domain_type = standardize(hint_domain_type);

  if (standardize(hint) == hint_domain_type) {
    return standardize(plist.get<std::string>(param, domain_type));

  } else if (domain_type == hint_domain_type) {
    return standardize(plist.get<std::string>(param, hint));

  } else if (Keys::starts_with(hint, hint_domain_type)) {
    Key default_domain;
    if (domain_type == "domain") {
      default_domain = hint.substr(hint_domain_type.size(), std::string::npos);
      if (Keys::starts_with(default_domain, "_"))
        default_domain = default_domain.substr(1, std::string::npos);
      if (Keys::starts_with(default_domain, dset_delimiter))
        default_domain = domain_type + default_domain;
    } else {
      default_domain = domain_type + hint.substr(hint_domain_type.size(), std::string::npos);
    }
    return standardize(plist.get<std::string>(param, default_domain));
  } else if (hint_domain_type == "domain") {
    return standardize(plist.get<std::string>(param, domain_type + "_" + hint));
  }
  return standardize(plist.get<std::string>(param));
}


Key
guessDomainType(const Key& domain)
{
  for (const auto& guess : { "snow", "canopy", "surface" }) {
    if (in(domain, guess)) { return guess; }
  }
  return "domain";
}


Key
readSuffix(Teuchos::ParameterList& list, const Key& basename, const Key& default_name)
{
  std::string basename_key_suffix_arg = basename + " key suffix";

  Key default_key;
  if (default_name.empty()) {
    return list.get<std::string>(basename_key_suffix_arg);
  } else {
    return list.get<std::string>(basename_key_suffix_arg, default_name);
  }
}


Key
readKey(Teuchos::ParameterList& list,
        const Key& domain,
        const Key& basename,
        const Key& default_name)
{
  std::string basename_key_arg = basename + " key";
  std::string basename_key_suffix_arg = basename + " key suffix";

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

Teuchos::Array<Key>
readKeys(Teuchos::ParameterList& list,
         const Key& domain,
         const Key& basename,
         Teuchos::Array<Key> const* const default_names)
{
  std::string basename_key_arg = basename + " keys";
  std::string basename_key_suffix_arg = basename + " key suffixes";

  if (list.isParameter(basename_key_suffix_arg)) {
    Teuchos::Array<std::string> suffixes =
      list.get<Teuchos::Array<std::string>>(basename_key_suffix_arg);
    Teuchos::Array<std::string> defaults(suffixes.size());
    int i = 0;
    for (const auto& suffix : suffixes) {
      defaults[i] = getKey(domain, suffix);
      ++i;
    }
    return list.get<Teuchos::Array<std::string>>(basename_key_arg, defaults);
  } else if (default_names) {
    Teuchos::Array<Key> default_keys(default_names->size());
    int i = 0;
    for (const auto& def : *default_names) {
      default_keys[i] = getKey(domain, def);
      ++i;
    }
    return list.get<Teuchos::Array<std::string>>(basename_key_arg, default_keys);
  } else {
    return list.get<Teuchos::Array<std::string>>(basename_key_arg);
  }
}


Tag
readTag(Teuchos::ParameterList& list, const Tag& default_tag)
{
  return readTag(list, "dependency tag", default_tag);
}


Tag
readTag(Teuchos::ParameterList& list, const std::string& param, const Tag& default_tag)
{
  std::string tag_str = list.get<std::string>(param, default_tag.get());
  return Tag{ tag_str };
}

} // namespace Keys
} // namespace Amanzi
