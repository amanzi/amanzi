//! Interface for Input implementations.
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  Defines an interface for vis writing.
*/

#ifndef AMANZI_INPUT_HH_
#define AMANZI_INPUT_HH_

#include <string>
#include <vector>
#include "AmanziTypes.hh"
#include "MeshDefs.hh"

namespace Amanzi {

// Interface for generic visualization.
class Input {
 public:

  virtual ~Input() {}

  // read data to file
  virtual void ReadField(Vector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const = 0;
  virtual void ReadField(IntVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const = 0;
  virtual void ReadFields(MultiVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const = 0;
  virtual void ReadFields(MultiVector_type& vec, const std::string& name, const std::vector<std::string>& subfield_names, const AmanziMesh::Entity_kind& location) const = 0;

  // can we template this (not yet...)
  virtual void ReadAttribute(double& val, const std::string& name) const = 0;
  virtual void ReadAttribute(int& val, const std::string& name) const = 0;
  virtual void ReadAttribute(std::string& val, const std::string& name) const = 0;
};


} // namespace Amanzi

#endif
