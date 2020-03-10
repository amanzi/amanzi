/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

namespace Amanzi {
namespace Functions {

//
// Computes function on a patch.
//
//template<class Device>
inline void
computeMeshFunction(const MultiFunction& f, double time, Patch& p)
{
  const AmanziMesh::Mesh* mesh = p.space.mesh.get();
  Kokkos::View<double**> txyz("txyz", mesh->space_dimension()+1, p.size());

  AmanziMesh::Entity_ID_List ids_list;
  mesh->get_set_entities(p.space.region, p.space.entity_kind,
                         p.space.ghosted ? AmanziMesh::Parallel_type::ALL : AmanziMesh::Parallel_type::OWNED,
                         ids_list);

  { // context for views
    // note this is a workaround until we have a way of getting the view on device
    AMANZI_ASSERT(ids_list.size() > 0);
    Kokkos::View<int*> ids("ids", ids_list.size());
    {
      Kokkos::View<int*, AmanziDefaultHost, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ids_host(ids_list.data(), ids_list.size());
      Kokkos::deep_copy(ids, ids_host);
    }
  
    if (p.space.entity_kind == AmanziMesh::NODE) {
      Errors::Message msg("computeMeshFunction on NODE not yet implemented (Mesh::node_get_coordinates() is not accessible on device)");
      throw(msg);

      // Kokkos::parallel_for(
      //     "computeMeshFunction txyz init node",
      //     p.size(),
      //     KOKKOS_LAMBDA(const int& i) {
      //       txyz(0,i) = time;
      //       auto cc = mesh->node_coordinate(ids[i]);
      //       txyz(1,i) = cc[0];
      //       txyz(2,i) = cc[1];
      //       if (mesh->space_dimension() == 3)
      //         txyz(3,i) = cc[2];
      //     });

    } else if (p.space.entity_kind == AmanziMesh::CELL) {

      Kokkos::parallel_for(
          "computeMeshFunction txyz init cell",
          p.size(),
          KOKKOS_LAMBDA(const int& i) {
            txyz(0,i) = time;
            auto cc = mesh->cell_centroid(ids[i]);
            txyz(1,i) = cc[0];
            txyz(2,i) = cc[1];
            if (mesh->space_dimension() == 3)
              txyz(3,i) = cc[2];
          });

    } else if (p.space.entity_kind == AmanziMesh::FACE) {

      Kokkos::parallel_for(
          "computeMeshFunction txyz init face",
          p.size(),
          KOKKOS_LAMBDA(const int& i) {
            txyz(0,i) = time;
            auto cc = mesh->face_centroid(ids[i]);
            txyz(1,i) = cc[0];
            txyz(2,i) = cc[1];
            if (mesh->space_dimension() == 3)
              txyz(3,i) = cc[2];
          });
    }
  }
  f.apply(txyz, p.data);
}

//
// Compute functions on a multi-patch.
//
//template<class Device>
inline void
computeMeshFunction(const std::vector<Teuchos::RCP<const MultiFunction>>& f,
                    double time, MultiPatch& mp)
{
  AMANZI_ASSERT(f.size() == mp.size());
  for (size_t i=0; i!=f.size(); ++i) {
    computeMeshFunction(*f[i], time, mp[i]);
  }  
}

//
// Compute set of functions on CompositeVector
//
//template<class Device>
inline void
computeMeshFunction(const std::vector<Teuchos::RCP<const MultiFunction>>& f,
                    double time, const MultiPatchSpace& mps, CompositeVector& cv)
{
  AMANZI_ASSERT(f.size() == mps.size());
  for (int i=0; i!=mps.size(); ++i) {
    std::string compname = AmanziMesh::entity_kind_string(mps[i].entity_kind);
    if (cv.HasComponent(compname)) {
      Patch p(mps[i]);
      computeMeshFunction(*f[i], time, p);
      patchToCompositeVector(p, compname, cv);
    }
  }
}


} // namespace Functions
} // namespace Amanzi


