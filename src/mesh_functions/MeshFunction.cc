/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>


#include "MeshFunction.hh"

namespace Amanzi {
namespace Functions {

//
// Computes function on a patch.
//
//template<class Device>
void
computeMeshFunction(const MultiFunction& f, double time, Patch& p)
{
  const AmanziMesh::Mesh* mesh = p.space.mesh.get();
  Kokkos::View<double**> txyz("txyz", mesh->space_dimension()+1, p.size());

  AmanziMesh::Entity_ID_List ids_list;
  mesh->get_set_entities(p.space.region, p.space.entity_kind,
                         p.space.ghosted ? AmanziMesh::Parallel_type::ALL : AmanziMesh::Parallel_type::OWNED,
                         ids_list);

  { // context for views
    // if empty, nothing to do
    if (ids_list.size() == 0) return;

    // note this is a workaround until we have a way of getting the view on device through the mesh interface
    Kokkos::View<int*, DeviceOnlyMemorySpace> ids("ids", ids_list.size());
    {
      Kokkos::View<int*, DefaultHost, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ids_host(ids_list.data(), ids_list.size());
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
void
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
void
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




//
// process a list for regions, components, and functions
//
Spec
processSpecWithFunction(Teuchos::ParameterList& list,
                        std::string function_name)
{
  Teuchos::Array<std::string> regions;
  if (list.isParameter("regions")) {
    regions = list.get<Teuchos::Array<std::string> >("regions");
  } else {
    regions.push_back(list.get<std::string>("region", Keys::cleanPListName(list.name())));
  }

  Teuchos::Array<std::string> components;
  if (list.isParameter("component")) {
    components.push_back(list.get<std::string>("component"));
  } else {
    std::vector<std::string> defs{"cell", "face", "node"};
    Teuchos::Array<std::string> def(defs);
    components = list.get<Teuchos::Array<std::string> >("components", def);
  }

  Teuchos::ParameterList f_list = list.sublist(function_name, true);

  // Make the function.
  Teuchos::RCP<Function> f;
  FunctionFactory f_fact;
  try {
    f = Teuchos::rcp(f_fact.Create(f_list));
  } catch (std::exception& msg) {
    Errors::Message m;
    m << "error in sublist " << function_name << ": " << msg.what();
    throw(m);
  }

  // TODO: Currently this does not support multiple-DoF functions.   --etc
  Teuchos::RCP<MultiFunction> func = Teuchos::rcp(new MultiFunction(f));

  return std::make_tuple(std::move(regions), std::move(components), func);
}

//
// parse a list of region-based specs
std::pair<MultiPatchSpace,
          std::vector<Teuchos::RCP<const MultiFunction>>>
processListWithFunction(Teuchos::ParameterList& list,
                        std::string function_name)
{
  MultiPatchSpace space;
  std::vector<Teuchos::RCP<const MultiFunction>> functions;

  // All are expected to be sublists of identical structure.
  for (auto sublist : list) {
    std::string name = sublist.first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList spec_plist = list.sublist(name);


      Spec spec;
      try {
        spec = processSpecWithFunction(spec_plist, function_name);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist " << name << ": " << msg.what();
        throw(m);
      }

      for (const auto& region : std::get<0>(spec)) {
        for (const auto& comp : std::get<1>(spec)) {
          AmanziMesh::Entity_kind entity_kind;
          try {
            entity_kind = AmanziMesh::entity_kind(comp);
          } catch (Errors::Message& msg) {
            Errors::Message m;
            m << "component \"" << comp << "\" provided in sublist \""
              << function_name << "\" must be one of \"cell\", \"face\", or \"node\".";
            throw(m);
          }
          space.AddPatch(region, entity_kind, 1);
          functions.emplace_back(std::get<2>(spec));
        }
      }
    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter " << name << " is not a sublist";
      throw(m);
    }
  }
  return std::make_pair(std::move(space), std::move(functions));
}



// void ProcessListWithoutFunction_(const Teuchos::ParameterList&,
//         const Teuchos::RCP<Functions::BoundaryFunction>&) const;



} // namespace Functions
} // namesapce Amanzi
