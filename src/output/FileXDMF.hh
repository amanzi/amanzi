/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Reads and writes XDMF files.

/*
  Writes XDMF metadata for file i/o.
*/

#ifndef AMANZI_FILE_XDMF_HH_
#define AMANZI_FILE_XDMF_HH_


#include "AmanziTypes.hh"
#include "MeshDefs.hh"

namespace Amanzi {

class FileXDMF {
 public:
  FileXDMF(const std::string& filename, int n_nodes, int n_elems, int n_conn);

  void
  CreateTimestep(const double time, const int cycle, const bool write_mesh);
  void CloseTimestep(const double time, const int cycle, const bool write_mesh);

  template <typename scalar_type>
  void WriteField(const std::string& var_name,
                  const AmanziMesh::Entity_kind& location,
                  const std::string& dofnum = "0");
  template <typename scalar_type>
  void WriteFields(const std::string& var_name,
                   const std::vector<std::string>& subfield_names,
                   const AmanziMesh::Entity_kind& location);
  template <typename scalar_type>
  void WriteFields(const std::string& var_name, int n_dofs,
                   const AmanziMesh::Entity_kind& location);

 private:
  void WriteFile_(const std::string& filename, const Teuchos::XMLObject& xml);
  Teuchos::XMLObject GetMeshTimestep_(const double time, const int cycle);
  Teuchos::XMLObject GetMeshVisit_();
  Teuchos::XMLObject
  GetDataTimestep_(const double time, const int cycle, const bool write_mesh);
  Teuchos::XMLObject GetDataVisit_();
  Teuchos::XMLObject GetHeaderGlobal_();
  Teuchos::XMLObject GetHeaderLocal_(const std::string name, const double time,
                                     const int cycle, bool write_mesh);
  Teuchos::XMLObject GetTopo_(const int cycle, const bool write_mesh);
  Teuchos::XMLObject GetGeo_(const int cycle, const bool write_mesh);
  Teuchos::XMLObject GetGrid_(std::string grid_filename);
  Teuchos::XMLObject FindGridNode_(Teuchos::XMLObject xmlobject);
  Teuchos::XMLObject FindMeshNode_(Teuchos::XMLObject xmlobject);

  Teuchos::XMLObject
  GetField_(const std::string& varname, const std::string& location,
            const std::string& data_type, int length, int cycle);

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
FileXDMF::WriteField(const std::string& var_name,
                     const AmanziMesh::Entity_kind& location,
                     const std::string& dofnum)
{
  auto node = FindMeshNode_(xmf_data_timestep_);
  std::string data_type = FileXDMF::DatatypeMap<scalar_type>::type;

  if (location == AmanziMesh::Entity_kind::CELL) {
    node.addChild(
      GetField_(var_name + "." + dofnum, "Cell", data_type, n_elems_, cycle_));
  } else if (location == AmanziMesh::Entity_kind::NODE) {
    node.addChild(
      GetField_(var_name + "." + dofnum, "Node", data_type, n_nodes_, cycle_));
  } else {
    Errors::Message message;
    message << "FileXDMF: Cannot write location of type \""
            << AmanziMesh::entity_kind_string(location) << "\"";
    Exceptions::amanzi_throw(message);
  }
}

template <typename scalar_type>
void
FileXDMF::WriteFields(const std::string& varname,
                      const std::vector<std::string>& subfield_names,
                      const AmanziMesh::Entity_kind& location)
{
  for (auto subname : subfield_names) {
    WriteField<scalar_type>(varname, location, subname);
  }
}

template <typename scalar_type>
void
FileXDMF::WriteFields(const std::string& varname, int n_dofs,
                      const AmanziMesh::Entity_kind& location)
{
  for (int i = 0; i != n_dofs; ++i) {
    WriteField<scalar_type>(varname, location, std::to_string(i));
  }
}


} // namespace Amanzi


#endif
