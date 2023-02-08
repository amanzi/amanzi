/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins, others
*/

#pragma once

#include <string>
#include <vector>

#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziTypes.hh"
#include "VerboseObject.hh"
#include "Mesh.hh"
#include "MeshDefs.hh"
#include "MeshFrameworkTraits.hh"

namespace Amanzi {

namespace AmanziGeometry {
class GeometricModel;
}

namespace AmanziMesh {
class MeshFramework;
class MeshColumn;

// -------------------------------------------------------------
// Factory for creating a MeshColumn object from a parent and a column ID
// -------------------------------------------------------------
Teuchos::RCP<Mesh>
createColumnMesh(const Teuchos::RCP<const Mesh>& parent_mesh,
                 int col_id,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);

// -------------------------------------------------------------
//  class MeshFrameworkFactory
// -------------------------------------------------------------
class MeshFrameworkFactory {
 public:
  explicit MeshFrameworkFactory(const Comm_ptr_type& comm,
                       const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm=Teuchos::null,
                       const Teuchos::RCP<Teuchos::ParameterList>& plist=Teuchos::null);

  // undefined copy constructor to avoid unwanted copies
  MeshFrameworkFactory(MeshFrameworkFactory& old) = delete;

  ~MeshFrameworkFactory() = default;

  // Get/set the framework preference, an ordered list of available frameworks.
  //
  // During construction, the first framework that provides the needed
  // capability is used (except in extraction -- see that documentation).
  void set_preference(const Preference& pref);
  const Preference& get_preference() const { return preference_; }

  // Get/set the geometric model
  Teuchos::RCP<const AmanziGeometry::GeometricModel> geometric_model() const { return gm_; }
  void set_geometric_model(const Teuchos::RCP<const AmanziGeometry::GeometricModel>& gm) { gm_ = gm; }

  // Get/set the parameter list
  Teuchos::RCP<const Teuchos::ParameterList> parameter_list() const { return plist_; }
  void set_parameter_list(const Teuchos::RCP<const Teuchos::ParameterList>& plist);

  // Creation methods
  // -- Create a mesh by reading the specified file (or set of files)
  Teuchos::RCP<MeshFramework> create(const std::string& filename);

  // -- Generate a hex mesh (3D)
  //    Generates a structured mesh covering [x0,x1] X [y0,y1] X [z0,z1] with (nx, ny, nz) cells.
  Teuchos::RCP<MeshFramework> create(const double x0, const double y0, const double z0,
          const double x1, const double y1, const double z1,
          const int nx, const int ny, const int nz);

  // -- Generate a quad mesh (2D)
  //    Generates a structured mesh covering [x0,x1] X [y0,y1] with (nx, ny) cells.
  Teuchos::RCP<MeshFramework> create(const double x0, const double y0,
          const double x1, const double y1,
          const int nx, const int ny);

  // -- Generate a mesh from parameters in a list
  Teuchos::RCP<MeshFramework> create(const Teuchos::ParameterList& gen_plist);

  // -- Extract a mesh from another mesh and a collection of entities
  //    Lifts setids of type setkind from the parent mesh and makes a new mesh
  //    out of those entities (and their closure).  If parent_mesh is 3D and
  //    setkind is FACE, the mesh may be flattened to form a 2D mesh.
  //
  //    Note that preference is ignored here in favor of the type of parent_mesh.
  //    Frameworks that extract often need to assume the parent_mesh is the same
  //    type as the constructed mesh in order to use internal information.  Note
  //    this could likely be relaxed, but is easiest to require for now.
  Teuchos::RCP<MeshFramework> create(const Teuchos::RCP<const Mesh>& parent_mesh,
          const Entity_ID_View& setids,
          const Entity_kind setkind,
          const bool flatten=false);

  // -- Extract a mesh from another mesh and a set in that mesh.
  //    If expert parameter "extraction method" is missing, then it lifts setids
  //    of type setkind from the parent mesh's sets named in setnames and
  //    make a new mesh out of those entities (and their closure).
  //    If parent_mesh is 3D and setkind is FACE, the mesh may be flattened to
  //    form a 2D mesh.
  //
  //    Note that preference is ignored here in favor of the type of parent_mesh.
  //    Frameworks that extract often need to assume the parent_mesh is the same
  //    type as the constructed mesh in order to use internal information.  Note
  //    this could likely be relaxed, but is easiest to require for now.
  //
  //    If expert paraeter "extraction method" is present, then the specified
  //    value is use to build a new mesh.
  Teuchos::RCP<MeshFramework> create(const Teuchos::RCP<const Mesh>& parent_mesh,
          const std::vector<std::string>& setnames,
          const Entity_kind setkind,
          const bool flatten=false);

 protected:
  // The parallel environment
  Comm_ptr_type comm_;

  // A list of preferred mesh frameworks to consider
  Preference preference_;

  // Object encoding the level of verbosity and output stream for diagnostic
  // messages, other control parameters that are NOT about construction.
  Teuchos::RCP<Teuchos::ParameterList> plist_;

  // The geometric model describing the space within which the mesh lives
  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_;

  // logging
  Teuchos::RCP<VerboseObject> vo_;

 private:
  std::string extraction_method_;
};

} // namespace AmanziMesh
} // namespace Amanzi


