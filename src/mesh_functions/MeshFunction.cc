/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Key.hh"
#include "MeshDefs.hh"
#include "DataStructuresHelpers.hh"
#include "MeshFunction.hh"

namespace Amanzi {
namespace Functions {

MeshFunction::MeshFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                           AmanziMesh::Entity_kind entity_kind)
  : mesh_(mesh), entity_kind_(entity_kind)
{}

MeshFunction::MeshFunction(Teuchos::ParameterList& list,
                           const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                           const std::string& function_name,
                           AmanziMesh::Entity_kind entity_kind,
                           int flag)
  : mesh_(mesh), entity_kind_(entity_kind), flag_(flag)
{
  // All are expected to be sublists of identical structure.
  for (const auto& sublist : list) {
    std::string name = sublist.first;
    if (list.isSublist(name)) {
      Teuchos::ParameterList& spec_plist = list.sublist(name);

      try {
        readSpec_(spec_plist, function_name);
      } catch (Errors::Message& msg) {
        Errors::Message m;
        m << "in sublist " << name << ": " << msg.what();
        throw(m);
      }

    } else { // ERROR -- parameter is not a sublist
      Errors::Message m;
      m << "parameter " << name << " is not a sublist";
      throw(m);
    }
  }
}


// add a spec -- others may inherit this and overload to do some checking?
void
MeshFunction::addSpec(const Spec& spec)
{
  if (mesh_ == Teuchos::null) setMesh(std::get<1>(spec)->mesh);
  AMANZI_ASSERT(std::get<1>(spec)->mesh == mesh_);
  spec_list_.push_back(spec);
}


// data creation
Teuchos::RCP<CompositeVectorSpace>
MeshFunction::createCVS(bool ghosted) const
{
  AMANZI_ASSERT(mesh_ != Teuchos::null);
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh_)->SetGhosted(ghosted);
  for (auto [compname, ps, functor] : *this) {
    cvs->AddComponent(compname, ps->entity_kind, ps->num_vectors);
  }
  return cvs;
};

Teuchos::RCP<MultiPatchSpace>
MeshFunction::createMPS(bool ghosted) const
{
  auto mps = Teuchos::rcp(new MultiPatchSpace(mesh_, ghosted));
  for (auto spec : *this) mps->addPatch(std::get<1>(spec));
  return mps;
}


void
MeshFunction::readSpec_(Teuchos::ParameterList& list, const std::string& function_name)
{
  Teuchos::Array<std::string> regions;
  if (list.isParameter("regions")) {
    regions = list.get<Teuchos::Array<std::string>>("regions");
  } else {
    regions.push_back(list.get<std::string>("region", Keys::cleanPListName(list.name())));
  }

  Teuchos::Array<std::string> components;
  if (list.isParameter("component")) {
    components.push_back(list.get<std::string>("component"));
  } else {
    Teuchos::Array<std::string> defs;
    if (entity_kind_ != AmanziMesh::Entity_kind::UNKNOWN) {
      defs = { to_string(entity_kind_) };
    } else {
      defs = { "cell", "face", "node" };
    }
    components = list.get<Teuchos::Array<std::string>>("components", defs);
  }

  Teuchos::Array<AmanziMesh::Entity_kind> entity_kinds;

  if (entity_kind_ == AmanziMesh::Entity_kind::UNKNOWN) {
    if (list.isParameter("entity kind")) {
      entity_kinds.push_back(AmanziMesh::createEntityKind(list.get<std::string>("entity kind")));
    } else if (list.isParameter("entity kinds")) {
      auto ekinds = list.get<Teuchos::Array<std::string>>("entity kinds");
      for (auto ekind : ekinds) entity_kinds.push_back(AmanziMesh::createEntityKind(ekind));
    } else {
      for (auto compname : components) {
        AmanziMesh::Entity_kind ekind;
        try {
          ekind = AmanziMesh::createEntityKind(compname);
        } catch (std::exception& msg) {
          Errors::Message m;
          m << "error in sublist " << function_name << ": " << msg.what();
          throw(m);
        }
        entity_kinds.push_back(ekind);
      }
    }

    if (entity_kinds.size() != components.size()) {
      Errors::Message m;
      m << "error in sublist " << function_name
        << ": \"components\" and \"entity kinds\" contain lists of differing lengths.";
      throw(m);
    }

  } else {
    entity_kinds.resize(components.size(), entity_kind_);
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

  // TODO: Currently this does not support multiple-DoF functions.   --ETC
  Teuchos::RCP<MultiFunction> func = Teuchos::rcp(new MultiFunction(f));

  for (int i = 0; i != components.size(); ++i) {
    auto comp = components[i];
    auto ekind = entity_kinds[i];
    for (auto region : regions) { addSpec(comp, ekind, 1, region, func); }
  }
}


void
MeshFunction::Compute(double time, MultiPatch<double>& mp)
{
  AMANZI_ASSERT(mp.size() == spec_list_.size());
  for (int i = 0; i != mp.size(); ++i) {
    auto [name, ps, multi_func_ptr] = spec_list_[i];
    Impl::computeFunction(*multi_func_ptr, time, mp[i]);
  }
}


namespace Impl {

//
// helper function to get coordinates, txyz
//
Kokkos::View<double**>
getMeshFunctionCoordinates(double time, const PatchSpace& ps)
{
  const AmanziMesh::Mesh* mesh = ps.mesh.get();
  int dim = mesh->getSpaceDimension();
  Kokkos::View<double**> txyz("txyz", dim + 1, ps.size());

  auto ids = ps.getIDs();

  // if empty, nothing to do
  if (ids.size() == 0) return Kokkos::View<double**>();

  if (ps.entity_kind == AmanziMesh::NODE) {
    Kokkos::parallel_for(
      "computeMeshFunction txyz init node", ps.size(), KOKKOS_LAMBDA(const int& i) {
        txyz(0, i) = time;
        auto cc = mesh->getNodeCoordinate(ids[i]);
        txyz(1, i) = cc[0];
        txyz(2, i) = cc[1];
        if (mesh->getSpaceDimension() == 3) txyz(3, i) = cc[2];
      });

  } else if (ps.entity_kind == AmanziMesh::CELL) {
    Kokkos::parallel_for(
      "computeMeshFunction txyz init cell", ps.size(), KOKKOS_LAMBDA(const int& i) {
        txyz(0, i) = time;
        auto cc = mesh->getCellCentroid(ids[i]);
        txyz(1, i) = cc[0];
        txyz(2, i) = cc[1];
        if (dim == 3) txyz(3, i) = cc[2];
      });

  } else if (ps.entity_kind == AmanziMesh::FACE) {
    Kokkos::parallel_for(
      "computeMeshFunction txyz init face", ps.size(), KOKKOS_LAMBDA(const int& i) {
        txyz(0, i) = time;
        auto cc = mesh->getFaceCentroid(ids[i]);
        txyz(1, i) = cc[0];
        txyz(2, i) = cc[1];
        if (dim == 3) txyz(3, i) = cc[2];
      });
  }
  return txyz;
}


//
// Computes function of t,x,y,{z} on a patch
//
// template<class Device>
void
computeFunction(const MultiFunction& f, double time, Patch<double>& p)
{
  auto txyz = getMeshFunctionCoordinates(time, *p.space);
  Kokkos::fence();
  f.apply(txyz, p.data);
}


//
// Computes function of t,x,y,{z} on patch, sticking the result into a vector
//
// template<class Device>
void
computeFunction(const MultiFunction& f, double time, const PatchSpace& ps, CompositeVector& cv)
{
  AMANZI_ASSERT(ps.mesh == cv.getMesh()); // precondition -- same mesh
  auto txyz = getMeshFunctionCoordinates(time, ps);

  auto ids = ps.getIDs();

  Kokkos::fence();
  Kokkos::View<double**, Kokkos::LayoutLeft> cv_v =
    cv.viewComponent(AmanziMesh::to_string(ps.entity_kind), ps.ghosted);
  f.apply(txyz, cv_v, &ids);
}


void
copyFlags(const PatchSpace& ps, CompositeVector_<int>& flag_vec)
{
  AMANZI_ASSERT(ps.mesh == flag_vec.getMesh()); // precondition -- same mesh
  auto ids = ps.getIDs();
  auto flag_v = flag_vec.viewComponent(AmanziMesh::to_string(ps.entity_kind), ps.ghosted);

  if (ps.flag_type == -1) {
    Kokkos::parallel_for(
      "MeshFunctions::copyFlags", ps.size(), KOKKOS_LAMBDA(const int& i) {
        flag_v(ids(i), 0) = ps.flags(i);
      });
  } else {
    Kokkos::parallel_for(
      "MeshFunctions::copyFlags", ps.size(), KOKKOS_LAMBDA(const int& i) {
        flag_v(ids(i), 0) = ps.flag_type;
      });
  }
}

void
copyFlags(const MultiPatchSpace& mps, CompositeVector_<int>& flag_vec)
{
  for (auto& ps : mps) { copyFlags(*ps, flag_vec); }
}

} // namespace Impl
} // namespace Functions
} // namespace Amanzi
