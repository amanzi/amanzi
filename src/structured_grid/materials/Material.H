/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _Material_H_
#define _Material_H_

#include <PArray.H>
#include <Region.H>
#include <Property.H>

class Material
{
public:
  Material(const std::string&    name,
	   const Array<const Region*>& regions,
	   const std::vector<Property*>& properties = std::vector<Property*>());

  Material(const Material& rhs);
  ~Material();

  template <typename T>
  void setVal(BaseFab<T>& fab,T val,int comp, const Real* dx) const
    {
      int ng = 0;
      for (int k=0; k<regions.size(); ++k) {
        regions[k]->setVal(fab,val,comp,dx,ng);
      }
    }


  const Array<const Region*>& Regions() const {return regions;}
  const std::string& Name() const {return name;}
  const Property* Prop(const std::string& pname) const;
  Array<std::string> PropertyNames() const;

protected:
  void ClearProperties();

  std::string name;
  Array<const Region*> regions;
  std::vector<Property*> properties;
};
#endif
