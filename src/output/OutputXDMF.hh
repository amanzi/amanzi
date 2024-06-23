/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! OutputXDMF: writes XDMF+H5 files for visualization in VisIt
/*
  XDMF implementation of an Output object, can only work as a Vis object as it
  needs a mesh and cannot handle face DoFs.
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

inline int
XDMFCellTypeID(AmanziMesh::Cell_kind type)
{
  // cell type id's defined in Xdmf/include/XdmfTopology.h
  switch (type) {
  case AmanziMesh::Cell_kind::POLYGON:
    return 3;
  case AmanziMesh::Cell_kind::TRI:
    return 4;
  case AmanziMesh::Cell_kind::QUAD:
    return 5;
  case AmanziMesh::Cell_kind::TET:
    return 6;
  case AmanziMesh::Cell_kind::PYRAMID:
    return 7;
  case AmanziMesh::Cell_kind::PRISM:
    return 8; // wedge
  case AmanziMesh::Cell_kind::HEX:
    return 9;
  case AmanziMesh::Cell_kind::POLYHED:
    return 3; // for now same as polygon
  default:
    return 3; // unknown, for now same as polygon
  }
}


class OutputXDMF : public Output {
 public:
  OutputXDMF(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  virtual ~OutputXDMF() = default;

  // open and close files
  virtual void createTimestep(double time, int cycle) override;
  virtual void finalizeTimestep() override;
  virtual std::string getFilename(int cycle) const override { return h5_data_filename_; }

  // write data to file
  using Output::write;
  virtual void write(const Teuchos::ParameterList& attrs, const double& val) const override;
  virtual void write(const Teuchos::ParameterList& attrs, const int& val) const override;
  virtual void write(const Teuchos::ParameterList& attrs, const std::string& val) const override;
  virtual void write(const Teuchos::ParameterList& attrs, const Vector_type& vec) const override;
  virtual void write(const Teuchos::ParameterList& attrs, const IntVector_type& vec) const override;
  virtual void write(const Teuchos::ParameterList& attrs,
                     const CompositeVector_<double_type>& vec) const override;
  virtual void
  write(const Teuchos::ParameterList& attrs, const CompositeVector_<int_type>& vec) const override;

 protected:
  std::tuple<int, int, int> writeMesh_(int cycle);
  std::tuple<int, int, int> writeDualMesh_(int cycle);

 protected:
  bool is_dynamic_;
  bool init_;

  int cycle_;
  double time_;
  bool write_mesh_;

  std::string filenamebase_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  std::unique_ptr<FileHDF5> h5_data_;
  std::string h5_data_filename_;
  std::unique_ptr<FileHDF5> h5_mesh_;
  std::string h5_mesh_filename_;
  std::unique_ptr<FileXDMF> xdmf_;
  std::string xdmf_filename_;
};

} // namespace Amanzi


#endif // OUTPUT_XDMF_HH_
