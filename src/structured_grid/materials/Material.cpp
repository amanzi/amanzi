/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <Material.H>

Material::Material(const std::string&    _name,
		   const Array<const Region*>& _regions,
		   const std::vector<Property*>& _properties)
  : name(_name)
{
  regions.resize(_regions.size());
  for (int i=0; i<regions.size(); ++i) {
    regions[i] = _regions[i];
  }
  ClearProperties();
  int nprop = _properties.size();
  properties.resize(nprop);
  for (int i=0; i<nprop; ++i) {
    properties[i] = _properties[i]->clone();
  }
}

const Property*
Material::Prop(const std::string& pname) const
{
  int ip=-1;
  for (int i=0; i<properties.size()&&ip<0; ++i) {
    if (properties[i]->Name() == pname) {
      return properties[i];
    }
  }
  std::string str = "No property registered with name: " + pname;
  BoxLib::Abort(str.c_str());
  return 0;
}

Array<std::string>
Material::PropertyNames() const
{
  int nprop = properties.size();
  Array<std::string> names(nprop);
  for (int i=0; i<nprop; ++i) {
    names[i] = properties[i]->Name();
  }
  return names;
}

void
Material::ClearProperties()
{
  for (int i=0, nprop = properties.size(); i<nprop; ++i) {
    delete properties[i];
  }
  properties.clear();
}

Material::Material(const Material& rhs)
{
  name = rhs.name;
  if (rhs.regions.size()>0) {
    regions.resize(rhs.regions.size());
    for (int i=0; i<regions.size(); ++i) {
      regions[i] = rhs.regions[i];
    }
  }
  ClearProperties();
  int nprop = rhs.properties.size();
  properties.resize(nprop);
  for (int i=0; i<nprop; ++i) {
    properties[i] = rhs.properties[i]->clone();
  }
}

Material::~Material()
{
  ClearProperties();
}
