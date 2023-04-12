/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//!
#ifndef DATA_STRUCTURE_HELPERS_HH_
#define DATA_STRUCTURE_HELPERS_HH_

#include "AmanziTypes.hh"
#include "AmanziMap.hh"
#include "Patch.hh"
#include "CompositeVector.hh"

namespace Amanzi {

template<typename T>
void
copyPatchToCompositeVector(const Patch<T>& p, const std::string& component, CompositeVector_<T>& cv)
{
  auto cv_c = cv.viewComponent(component, p.space->ghosted);

  const auto& mesh = cv.getMesh();
  auto ids = mesh->getSetEntities(p.space->region,
          p.space->entity_kind,
          p.space->ghosted ? AmanziMesh::Parallel_kind::ALL : AmanziMesh::Parallel_kind::OWNED);

  if (component != "boundary_face") {
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> range({ 0, 0 }, { p.data.extent(0), p.data.extent(1) });
    Kokkos::parallel_for(
      "patchToCompositeVector", range, KOKKOS_LAMBDA(const int& i, const int& j) {
        cv_c(ids[i], j) = p.data(i, j);
      });
  } else {
    AMANZI_ASSERT(false && "Not yet implemented: patchToCompositeVector with boundary_face");
    // have to do some dancing here... this is not correct because p.data is
    // based on faces, but component is based on boundary faces.  Need to
    // either create temporary space, then import, or more likely, unpack the
    // mapping to make sure we only access cv_c on boundary faces.
    Kokkos::MDRangePolicy<Kokkos::Rank<2>> range({ 0, 0 }, { p.data.extent(0), p.data.extent(1) });
    Kokkos::parallel_for(
      "patchToCompositeVector boundary_face", range, KOKKOS_LAMBDA(const int& i, const int& j) {
        cv_c(ids[i], j) = p.data(i, j);
      });
  }
}

template<typename T>
void
copyPatchToCompositeVector(const Patch<T>& p,
                       const std::string& component,
                       CompositeVector_<T>& cv,
                       CompositeVector_<int>& flag_cv)
{
  auto ids = cv.getMesh()->getSetEntities(p.space->region,
          p.space->entity_kind,
          p.space->ghosted ? AmanziMesh::Parallel_kind::ALL :
          AmanziMesh::Parallel_kind::OWNED);

  if (component != "boundary_face") {
    // AMANZI_ASSERT(ids.extent(0) == p.data.extent(0));
    auto flag_type = p.space->flag_type;

    auto cv_c = cv.viewComponent(component, p.space->ghosted);
    auto flag_c = flag_cv.viewComponent(component, p.space->ghosted);

    Kokkos::parallel_for(
      "patchToCompositeVector", p.data.extent(0), KOKKOS_LAMBDA(const int& i) {
        cv_c(ids[i], 0) = p.data(i, 0);
        flag_c(ids[i], 0) = flag_type;
      });
  } else {
    AMANZI_ASSERT(false && "Not yet implemented: patchToCompositeVector with boundary_face");
    // have to do some dancing here... this is not correct because p.data is
    // based on faces, but component is based on boundary faces.  Need to
    // either create temporary space, then import, or more likely, unpack the
    // mapping to make sure we only access cv_c on boundary faces.
    // Kokkos::parallel_for(
    //     "patchToCompositeVector boundary_face",
    //     p.data.extent(0),
    //     KOKKOS_LAMBDA(const int& i) {
    //       cv_c(ids(i),0) = p.data(i,0);
    //       flag_c(ids(i), 0) = p.space->flag_type;
    //     });
  }
}


//
// Copies values from a set of patches into a vector.
//
template<typename T>
void
copyMultiPatchToCompositeVector(const MultiPatch<T>& mp, const std::string& component, CompositeVector_<T>& cv)
{
  for (const auto& p : mp) { copyPatchToCompositeVector<T>(p, component, cv); }
}

//
// Copies values and flag from a set of patches into a vector and a flag vector.
//
template<typename T>
void
copyMultiPatchToCompositeVector(const MultiPatch<T>& mp,
                                        const std::string& component,
                                        CompositeVector_<T>& cv,
                                        CompositeVector_<int>& flag)
{
  for (const auto& p : mp) { copyPatchToCompositeVector<T>(p, component, cv, flag); }
}


void
DeriveFaceValuesFromCellValues(CompositeVector&);


// Create a BFS-ordered list of TreeVector(Space) nodes.
template <class T>
void
recurseTreeVectorBFS(T& tv, std::vector<Teuchos::RCP<T>>& list)
{
  for (typename T::iterator it = tv.begin(); it != tv.end(); ++it) { list.push_back(*it); }

  for (typename T::iterator it = tv.begin(); it != tv.end(); ++it) {
    recurseTreeVectorBFS<T>(**it, list);
  }
}

template <class T>
void
recurseTreeVectorBFS_const(const T& tv, std::vector<Teuchos::RCP<const T>>& list)
{
  for (const auto& it : tv) { list.push_back(it); }

  for (const auto& it : tv) {
    recurseTreeVectorBFS_const<T>(*it, list);
  }
}

// Create a list of leaf nodes of the TreeVector(Space)
template <class T>
std::vector<Teuchos::RCP<T>>
collectTreeVectorLeaves(T& tv)
{
  std::vector<Teuchos::RCP<T>> list;
  list.push_back(Teuchos::rcpFromRef(tv));
  recurseTreeVectorBFS<T>(tv, list);

  std::vector<Teuchos::RCP<T>> leaves;
  for (const auto& it : tv) {
    if (it->getData() != Teuchos::null) { leaves.push_back(it); }
  }

  return leaves;
}

template <class T>
std::vector<Teuchos::RCP<const T>>
collectTreeVectorLeaves_const(const T& tv)
{
  std::vector<Teuchos::RCP<const T>> list;
  list.push_back(Teuchos::rcpFromRef(tv));
  recurseTreeVectorBFS_const<T>(tv, list);

  std::vector<Teuchos::RCP<const T>> leaves;
  for (const auto& it : list) {
    if (it->getData() != Teuchos::null) { leaves.push_back(it); }
  }

  return leaves;
}

template <class T>
int
getNumTreeVectorLeaves(const T& tv)
{
  return collectTreeVectorLeaves_const(tv).size();
}


void
copyMeshCoordinatesToVector(const AmanziMesh::Mesh& mesh,
                            AmanziMesh::Entity_kind kind,
                            CompositeVector& vec);
void
copyVectorToMeshCoordinates(const CompositeVector& vec,
                            AmanziMesh::Mesh& mesh);

} // namespace Amanzi

#endif
