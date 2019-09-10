//! Interface for Output implementations.
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

#ifndef AMANZI_OUTPUT_HH_
#define AMANZI_OUTPUT_HH_

#include <string>
#include <vector>
#include "AmanziTypes.hh"
#include "MeshDefs.hh"

namespace Amanzi {

// Interface for generic visualization.
class Output {
 public:

  virtual ~Output() {}

  // open and close files
  virtual void InitializeCycle(double time, int cycle) = 0;
  virtual void FinalizeCycle() = 0;

  // write data to file
  virtual void WriteField(const Vector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const = 0;
  virtual void WriteField(const IntVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const = 0;
  virtual void WriteFields(const MultiVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const = 0;
  virtual void WriteFields(const MultiVector_type& vec, const std::string& name, const std::vector<std::string>& subfield_names, const AmanziMesh::Entity_kind& location) const = 0;

  // can we template this (not yet...)
  virtual void WriteAttribute(const double& val, const std::string& name) const = 0;
  virtual void WriteAttribute(const int& val, const std::string& name) const = 0;
  virtual void WriteAttribute(const std::string& val, const std::string& name) const = 0;

  virtual std::string Filename() const = 0;
};


} // namespace Amanzi

#endif
