//! OutputXDMF: writes XDMF+H5 files for VisIt
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  XDMF implementation of an Output object.
*/

#ifndef AMANZI_OUTPUT_XDMF_HH_
#define AMANZI_OUTPUT_XDMF_HH_

#include <string>
#include <vector>

#include "Teuchos_ParameterList.hpp"

#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "errors.hh"
#include "dbc.hh"
#include "Mesh.hh"

#include "FileHDF5.hh"
#include "FileXDMF.hh"
#include "Output.hh"

namespace Amanzi {

inline int XDMFCellTypeID(AmanziMesh::Cell_type type)
{
  // cell type id's defined in Xdmf/include/XdmfTopology.h
  AMANZI_ASSERT (AmanziMesh::cell_valid_type(type));
  switch (type)
  {
    case AmanziMesh::POLYGON:
      return 3;
    case AmanziMesh::TRI:
      return 4;
    case AmanziMesh::QUAD:
      return 5;
    case AmanziMesh::TET:
      return 6;
    case AmanziMesh::PYRAMID:
      return 7;
    case AmanziMesh::PRISM:
      return 8; //wedge
    case AmanziMesh::HEX:
      return 9;
    case AmanziMesh::POLYHED:
      return 3; //for now same as polygon
    default:
      return 3; //unknown, for now same as polygon
  }
}


class OutputXDMF : public Output {

 public:
  OutputXDMF(Teuchos::ParameterList& plist,
	     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  virtual ~OutputXDMF() = default;
  
  // open and close files
  virtual void InitializeCycle(double time, int cycle) override;
  virtual void FinalizeCycle() override;
  virtual std::string Filename() const override;

  // write data to file
  virtual void WriteField(const Vector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const override {
    xdmf_->WriteField<Vector_type::scalar_type>(name, location);
    std::stringstream path;
    path << "/" << name << ".0/" << cycle_;
    h5_data_->WriteVector<Vector_type::scalar_type>(vec, path.str());
  }
  virtual void WriteField(const IntVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const override {
    xdmf_->WriteField<IntVector_type::scalar_type>(name, location);
    std::stringstream path;
    path << "/" << name << ".0/" << cycle_;
    h5_data_->WriteVector<IntVector_type::scalar_type>(vec, path.str());
  }
    
  virtual void WriteFields(const MultiVector_type& vec, const std::string& name, const AmanziMesh::Entity_kind& location) const override {
    xdmf_->WriteFields<MultiVector_type::scalar_type>(name, vec.getNumVectors(), location);
    std::vector<std::string> paths;
    for (int i=0; i!=vec.getNumVectors(); ++i) {
      std::stringstream path;
      path << "/" << name << "." << i << "/" << cycle_;
      paths.push_back(path.str());
    }    
    h5_data_->WriteMultiVector<MultiVector_type::scalar_type>(vec, paths);
  }
  virtual void WriteFields(const MultiVector_type& vec, const std::string& name, const std::vector<std::string>& subfield_names, const AmanziMesh::Entity_kind& location) const override {
    AMANZI_ASSERT(subfield_names.size() == vec.getNumVectors());
    xdmf_->WriteFields<MultiVector_type::scalar_type>(name, subfield_names, location);

    std::vector<std::string> paths;
    for (int i=0; i!=vec.getNumVectors(); ++i) {
      std::stringstream path;
      path << "/" << name << "." << subfield_names[i] << "/" << cycle_;
      paths.push_back(path.str());
    }    
    h5_data_->WriteMultiVector<MultiVector_type::scalar_type>(vec, paths);
  }

  // can we template this (not yet...)
  virtual void WriteAttribute(const double& val, const std::string& name) const override {
    h5_data_->WriteAttribute(val, name, "/"); }
  virtual void WriteAttribute(const int& val, const std::string& name) const override {
    h5_data_->WriteAttribute(val, name, "/"); }
  virtual void WriteAttribute(const std::string& val, const std::string& name) const override {
    h5_data_->WriteAttribute(val, name, "/"); }
    
 protected:
  std::tuple<int,int,int> WriteMesh_(int cycle);
  
 protected:
  bool is_dynamic_;
  bool init_;

  int cycle_;
  double time_;
  bool write_mesh_;

  std::string filenamebase_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  Teuchos::RCP<FileHDF5> h5_data_;
  Teuchos::RCP<FileHDF5> h5_mesh_;
  Teuchos::RCP<FileXDMF> xdmf_;
};
  
} // namespace Amanzi

#endif  // OUTPUT_XDMF_HH_
