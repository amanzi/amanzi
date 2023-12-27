/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! reads and writes XDMF files.
/*
  writes XDMF metadata for file i/o.
*/

#ifndef AMANZI_FILE_XDMF_HH_
#define AMANZI_FILE_XDMF_HH_


#include "AmanziTypes.hh"
#include "MeshDefs.hh"

namespace Amanzi {

class FileXDMF {
 public:
  FileXDMF(const std::string& filename, int n_nodes, int n_elems, int n_conn);

  void createTimestep(const double time, const int cycle, const bool write_mesh);
  void finalizeTimestep(const double time, const int cycle, const bool write_mesh);

  template <typename scalar_type>
  void writeField(const std::string& var_name, const AmanziMesh::Entity_kind& location);
  template <typename scalar_type>
  void writeFields(const std::string& var_name,
                   const std::vector<std::string>& subfield_names,
                   const AmanziMesh::Entity_kind& location);
  template <typename scalar_type>
  void
  writeFields(const std::string& var_name, int n_dofs, const AmanziMesh::Entity_kind& location);

 private:
  void writeFile_(const std::string& filename, const Teuchos::XMLObject& xml);
  Teuchos::XMLObject getMeshTimestep_(const double time, const int cycle);
  Teuchos::XMLObject getMeshVisit_();
  Teuchos::XMLObject getDataTimestep_(const double time, const int cycle, const bool write_mesh);
  Teuchos::XMLObject getDataVisit_();
  Teuchos::XMLObject getHeaderGlobal_();
  Teuchos::XMLObject
  getHeaderLocal_(const std::string name, const double time, const int cycle, bool write_mesh);
  Teuchos::XMLObject getTopo_(const int cycle, const bool write_mesh);
  Teuchos::XMLObject getGeo_(const int cycle, const bool write_mesh);
  Teuchos::XMLObject getGrid_(std::string grid_filename);
  Teuchos::XMLObject FindGridNode_(Teuchos::XMLObject xmlobject);
  Teuchos::XMLObject FindMeshNode_(Teuchos::XMLObject xmlobject);

  Teuchos::XMLObject getField_(const std::string& varname,
                               const std::string& location,
                               const std::string& data_type,
                               int length,
                               int cycle);

 private:
  std::string filename_prefix_;
  std::string basename_;
  int cycle_;
  int last_mesh_cycle_;

  // mesh objects
  int n_elems_;
  int n_nodes_;
  int n_conns_;

  // xmf objects
  Teuchos::XMLObject xmf_data_visit_;
  Teuchos::XMLObject xmf_data_timestep_;
  Teuchos::XMLObject xmf_mesh_visit_;
  Teuchos::XMLObject xmf_mesh_timestep_;

  const static std::string xdmfHeader_;

 private:
  // a compile-time map from type to ASCEMIO::datatype_t
  template <typename T>
  struct DatatypeMap;
};

template <>
struct FileXDMF::DatatypeMap<int> {
  static const std::string type;
};
template <>
struct FileXDMF::DatatypeMap<float> {
  static const std::string type;
};
template <>
struct FileXDMF::DatatypeMap<double> {
  static const std::string type;
};


template <typename scalar_type>
void
FileXDMF::writeField(const std::string& var_name, const AmanziMesh::Entity_kind& location)
{
  auto node = FindMeshNode_(xmf_data_timestep_);
  std::string data_type = FileXDMF::DatatypeMap<scalar_type>::type;

  if (location == AmanziMesh::Entity_kind::CELL) {
    node.addChild(getField_(var_name, "Cell", data_type, n_elems_, cycle_));
  } else if (location == AmanziMesh::Entity_kind::NODE) {
    node.addChild(getField_(var_name, "Node", data_type, n_nodes_, cycle_));
  } else {
    Errors::Message message;
    message << "FileXDMF: Cannot write location of type \"" << AmanziMesh::to_string(location)
            << "\"";
    Exceptions::amanzi_throw(message);
  }
}

template <typename scalar_type>
void
FileXDMF::writeFields(const std::string& varname,
                      const std::vector<std::string>& subfield_names,
                      const AmanziMesh::Entity_kind& location)
{
  for (auto subname : subfield_names) {
    writeField<scalar_type>(varname + "." + subname, location);
  }
}

template <typename scalar_type>
void
FileXDMF::writeFields(const std::string& varname,
                      int n_dofs,
                      const AmanziMesh::Entity_kind& location)
{
  for (int i = 0; i != n_dofs; ++i) {
    writeField<scalar_type>(varname + "." + std::to_string(i), location);
  }
}


} // namespace Amanzi


#endif
