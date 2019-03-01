/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins, others
*/

#ifndef AMANZI_MESH_FACTORY_HH_
#define AMANZI_MESH_FACTORY_HH_

#include <string>
#include <vector>
#include "AmanziComm.hh"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MeshException.hh"
#include "MeshFramework.hh"
#include "Mesh.hh"

#include "GeometricModel.hh"

namespace Amanzi {
namespace AmanziMesh {

// -------------------------------------------------------------
//  class MeshFactory
// -------------------------------------------------------------
class MeshFactory {
 protected:

  // The parallel environment
  Comm_ptr_type comm_;

  // A list of preferred mesh frameworks to consider
  FrameworkPreference my_preference;

 private:

  // private, undefined copy constructor to avoid unwanted copies
  MeshFactory(MeshFactory& old);

  // Object encoding the level of verbosity and output stream for
  // diagnostic messages

  Teuchos::RCP<const Teuchos::ParameterList> plist_;

 public:

  // Default constructor.
  explicit MeshFactory(const Comm_ptr_type& comm, 
                       const Teuchos::RCP<const Teuchos::ParameterList>& plist = Teuchos::null);

  // Destructor
  ~MeshFactory() = default;

  // Get the framework preference
  const FrameworkPreference& preference(void) const
  { return my_preference; }

  // Set the framework preference
  void set_preference(const FrameworkPreference& pref);

  // Create a mesh by reading the specified file (or set of files)
  Teuchos::RCP<Mesh> create(const std::string& filename, 
                            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null, 
                            const bool request_faces = true,
                            const bool request_edges = false);

  // Generate a mesh from parameters
  Teuchos::RCP<Mesh> create(Teuchos::ParameterList& gen_plist,
                            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null, 
                            const bool request_faces = true,
                            const bool request_edges = false);

  

  // Create a hexahedral mesh of the specified dimensions
  Teuchos::RCP<Mesh> create(double x0, double y0, double z0,
                            double x1, double y1, double z1,
                            int nx, int ny, int nz, 
                            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null, 
                            const bool request_faces = true,
                            const bool request_edges = false);

    
  // Create a quadrilateral mesh of the specified dimensions
  Teuchos::RCP<Mesh> create(double x0, double y0,
                            double x1, double y1,
                            int nx, int ny,
                            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null, 
                            const bool request_faces = true,
                            const bool request_edges = false);

  /*    
  // Create a quadrilateral/hexahedral mesh using the specified parameter list
  Teuchos::RCP<Mesh> create(Teuchos::ParameterList& parameter_list, 
                            const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm = Teuchos::null, 
                            const bool request_faces = true,
                            const bool request_edges = false);
  */

  // Create a mesh by extract subsets of entities from an existing mesh
  Teuchos::RCP<Mesh> create(const Teuchos::RCP<const Mesh>& parent_mesh,
                            const Entity_ID_List& setids,
                            const Entity_kind setkind,
                            const bool flatten = false,
                            const bool extrude = false,
                            const bool request_faces = true,
                            const bool request_edges = false);

  // Create a mesh by extract subsets of entities from an existing mesh
  Teuchos::RCP<Mesh> create(const Teuchos::RCP<const Mesh>& parent_mesh,
                            const std::vector<std::string>& setnames,
                            const Entity_kind setkind,
                            const bool flatten = false,
                            const bool extrude = false,
                            const bool request_faces = true,
                            const bool request_edges = false);
  
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif
