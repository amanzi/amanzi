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

#include "Key.hh"

namespace Amanzi {
namespace Keys {

std::string
readKey(Teuchos::ParameterList& list, const Key& domain, const Key& basename,
        const Key& default_name)
{
  std::string basename_key_arg = basename + " key";
  std::string basename_key_suffix_arg = basename + " key suffix";

  Key default_key;
  if (list.isParameter(basename_key_suffix_arg)) {
    default_key =
      getKey(domain, list.get<std::string>(basename_key_suffix_arg));
    return list.get<std::string>(basename_key_arg, default_key);
  } else if (!default_name.empty()) {
    default_key = getKey(domain, default_name);
    return list.get<std::string>(basename_key_arg, default_key);
  } else {
    return list.get<std::string>(basename_key_arg);
  }
}

Teuchos::Array<Key>
readKeys(Teuchos::ParameterList& list, const Key& domain, const Key& basename,
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
    return list.get<Teuchos::Array<std::string>>(basename_key_arg,
                                                 default_keys);
  } else {
    return list.get<Teuchos::Array<std::string>>(basename_key_arg);
  }
}


bool
splitDomainSet(const Key& name,
               std::tuple<std::string, std::string, std::string>& result)
{
  auto domain_var = splitKey(name);
  std::get<2>(result) = domain_var.second;

  // make sure format is domain_set_name_#-var_name
  std::size_t pos = domain_var.first.find_last_of('_');
  if (pos == std::string::npos) { return false; }

  // check sure the # is a number
  std::string domain_int =
    domain_var.first.substr(pos + 1, domain_var.first.size());
  std::size_t not_translated;
  if (!(domain_int.size() == 1 && domain_int[0] == '*')) {
    if (!std::isdigit(domain_int[0])) return false;
    std::stoi(domain_int, &not_translated);
    if (not_translated != domain_int.size()) { return false; }
  }
  std::get<1>(result) = domain_int;
  std::get<0>(result) = domain_var.first.substr(0, pos);
  return true;
}


} // namespace Keys
} // namespace Amanzi
