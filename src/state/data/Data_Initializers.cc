//! Data_Initializers are functions that know how to initialize data.
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  Initialize from parameter list.
*/

#include "Data_Initializers.hh"

namespace Amanzi {
namespace Data_Initializers {

template <>
bool Initialize<double>(Teuchos::ParameterList &plist, const Teuchos::ParameterList& attrs, double &t)
{
  return InitializePrimitiveByValue(plist, attrs, t);
}
  
template <>
bool Initialize<int>(Teuchos::ParameterList &plist, const Teuchos::ParameterList& attrs, int &t)
{
  return InitializePrimitiveByValue(plist, attrs, t);
}

template <>
bool Initialize<std::string>(Teuchos::ParameterList &plist, const Teuchos::ParameterList& attrs, std::string &t)
{
  return InitializePrimitiveByValue(plist, attrs, t);
}



} // namespace
} // namespace
