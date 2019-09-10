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

#ifndef AMANZI_INPUT_OUTPUT_HDF5_HH_
#define AMANZI_INPUT_OUTPUT_HDF5_HH_

#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "AmanziTypes.hh"
#include "Output.hh"
#include "Input.hh"
#include "FileHDF5.hh"


namespace Amanzi {

class InputOutputHDF5 : public Input, public Output {
 public:
  InputOutputHDF5(Teuchos::ParameterList& plist, const Comm_ptr_type& comm) :
      filename_base_(plist.get<std::string>("file name base","checkpoint")),
      filename_digits_(plist.get<int>("file name digits", 5)),
      comm_(comm),
      cycle_(-1)
  {}
  
  virtual ~InputOutputHDF5() = default;

  // 
  virtual void InitializeCycle(double time, int cycle) {
    cycle_ = cycle;
    file_ = Teuchos::rcp(new FileHDF5(comm_, Filename(), FILE_CREATE)); }    
  virtual void FinalizeCycle() {
    file_ = Teuchos::null;
  }


  // read data to file
  virtual void ReadField(Vector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const {
    file_->ReadVector(vec, name);
  }
  virtual void ReadField(IntVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const {
    file_->ReadVector(vec, name);
  }
  virtual void ReadFields(MultiVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const {
    file_->ReadMultiVectorBlock(vec, name);
  }
  virtual void ReadFields(MultiVector_type& vec, const std::string& name, const std::vector<std::string>& subfield_names, const AmanziMesh::Entity_kind& location) const {
    file_->ReadMultiVectorBlock(vec, name);
  }

  // can we template this (not yet...)
  virtual void ReadAttribute(double& val, const std::string& name) const {
    val = file_->ReadAttribute<double>(name);
  }
  virtual void ReadAttribute(int& val, const std::string& name) const {
    val = file_->ReadAttribute<int>(name);
  }    
  virtual void ReadAttribute(std::string& val, const std::string& name) const {
    val = file_->ReadAttribute<std::string>(name);
  }    

  // write data to file
  virtual void WriteField(const Vector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const {
    file_->WriteVector(vec, name);
  }
  virtual void WriteField(const IntVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const {
    file_->WriteVector(vec, name);
  }
  virtual void WriteFields(const MultiVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const {
    file_->WriteMultiVectorBlock(vec, name);
  }
  virtual void WriteFields(const MultiVector_type& vec, const std::string& name, const std::vector<std::string>& subfield_names, const AmanziMesh::Entity_kind& location) const {
    file_->WriteMultiVectorBlock(vec, name);
  }

  // can we template this (not yet...)
  virtual void WriteAttribute(const double& val, const std::string& name) const {
    file_->WriteAttribute(val, name);
  }
  virtual void WriteAttribute(const int& val, const std::string& name) const {
    file_->WriteAttribute(val, name);
  }
  virtual void WriteAttribute(const std::string& val, const std::string& name) const {
    file_->WriteAttribute(val, name);
  }

  virtual std::string Filename() const {
    // create the restart file
    std::stringstream oss;
    oss.flush();
    oss << filename_base_;
    oss.fill('0');
    oss.width(filename_digits_);
    oss << std::right << cycle_;
    return oss.str();
  }

 protected:
  std::string filename_base_;
  int filename_digits_;
  Comm_ptr_type comm_;
  int cycle_;
  
  Teuchos::RCP<FileHDF5> file_;
};


} // namespace Amanzi

#endif
