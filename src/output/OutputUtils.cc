/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! OutputUtils: Utility functions for I/O.
#include "AmanziTypes.hh"
#include "AmanziVector.hh"

namespace Amanzi {
namespace OutputUtils {


Vector_type_<GO>
asVector(const Map_ptr_type& map)
{
  Vector_type_<GO> vec(map);
  {
    auto data = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
    for (int i = 0; i != data.size(); ++i) { data(i, 0) = map->getGlobalElement(i); }
  }
  return vec;
}


std::vector<std::string>
getNames(const Teuchos::ParameterList& attrs, std::size_t count)
{
  std::vector<std::string> names;
  bool always_write =
    attrs.isParameter("always write subfield dof") && attrs.get<bool>("always write subfield dof");
  if (attrs.isParameter("subfieldnames") &&
      attrs.get<Teuchos::Array<std::string>>("subfieldnames").size() == count) {
    auto subfield_names = attrs.get<Teuchos::Array<std::string>>("subfieldnames");
    for (int i = 0; i != count; ++i) { names.emplace_back(attrs.name() + "." + subfield_names[i]); }
  } else if (count > 1 || always_write) {
    for (int i = 0; i != count; ++i) { names.emplace_back(attrs.name() + "." + std::to_string(i)); }
  } else {
    names.emplace_back(attrs.name());
  }
  return names;
}


} // namespace OutputUtils
} // namespace Amanzi
