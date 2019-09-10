//! Input/Output via HDF5 file for checkpoint/restart.
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  Input/Output via HDF5 file for checkpoint/restart.
*/


#include "InputOutputHDF5.hh"

namespace Amanzi {

// open and close files
void InputOutputHDF5::InitializeCycle(double time, int cycle)
{
  file_ = Teuchos::rcp(new FileHDF5(comm_, Filename_(cycle), FILE_CREATE));
}

void InputOutputHDF5::FinalizeCycle()
{
  file_ = Teuchos::null;
}

// read data to file
void InputOutputHDF5::ReadField(Vector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const
{
  file_->ReadVector(vec, name);
}
      
void InputOutputHDF5::ReadField(IntVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const
{
  file_->ReadVector(vec, name);
}

void InputOutputHDF5::ReadFields(MultiVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const
{
  file_->ReadMultiVectorBlock(vec, name);
}

void InputOutputHDF5::ReadFields(MultiVector_type& vec, const std::string& name, const std::vector<std::string>& subfield_names, const AmanziMesh::Entity_kind& location) const
{
  Errors::Message message("InputOutputHDF5::Not implemented, use the Block version instead.");
  throw(message);
}

// can we template this (not yet...)
void InputOutputHDF5::ReadAttribute(double& val, const std::string& name) const
{
  val = file_->ReadAttribute<double>(name);
}
  void InputOutputHDF5::ReadAttribute(const int& val, const std::string& name) const
{
}
  void InputOutputHDF5::ReadAttribute(const std::string& val, const std::string& name) const
{
}

  // write data to file
  void InputOutputHDF5::WriteField(const Vector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const
{
}
  void InputOutputHDF5::WriteField(const IntVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const
{
}
  void InputOutputHDF5::WriteFields(const MultiVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const
{
}
  void InputOutputHDF5::WriteFields(const MultiVector_type& vec, const std::string& name, const std::vector<std::string>& subfield_names, const AmanziMesh::Entity_kind& location) const
{
}

  // can we template this (not yet...)
  void InputOutputHDF5::WriteAttribute(const double& val, const std::string& name) const
{
}
  void InputOutputHDF5::WriteAttribute(const int& val, const std::string& name) const
{
}
  void InputOutputHDF5::WriteAttribute(const std::string& val, const std::string& name) const
{
}

};


} // namespace Amanzi

