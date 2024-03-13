/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Utility functions on a TreeVector
   ------------------------------------------------------------------------- */

#ifndef AMANZI_DATA_STRUCTURES_TREEVECTOR_UTILS_HH_
#define AMANZI_DATA_STRUCTURES_TREEVECTOR_UTILS_HH_

#include <algorithm>

namespace Amanzi {


// Create a BFS-ordered list of TreeVector(Space) nodes.
template <class T>
void
recurseTreeVectorBFS(T& tv, std::vector<Teuchos::RCP<T>>& list)
{
  for (auto it : tv) list.emplace_back(it);
  for (auto it : tv) recurseTreeVectorBFS<T>(*it, list);
}

template <class T>
void
recurseTreeVectorBFS_const(const T& tv, std::vector<Teuchos::RCP<const T>>& list)
{
  for (auto it : tv) list.push_back(it);
  for (auto it : tv) recurseTreeVectorBFS_const<T>(*it, list);
}

// Create a list of leaf nodes of the TreeVector(Space)
template <class T>
std::vector<Teuchos::RCP<T>>
collectTreeVectorLeaves(T& tv)
{
  std::vector<Teuchos::RCP<T>> list;
  list.emplace_back(Teuchos::rcpFromRef(tv));
  recurseTreeVectorBFS<T>(tv, list);

  std::vector<Teuchos::RCP<T>> leaves;
  for (auto it : list) {
    if (it->Data() != Teuchos::null) leaves.push_back(it);
  }
  return leaves;
}

template <class T>
std::vector<Teuchos::RCP<const T>>
collectTreeVectorLeaves_const(const T& tv)
{
  std::vector<Teuchos::RCP<const T>> list;
  list.emplace_back(Teuchos::rcpFromRef(tv));
  recurseTreeVectorBFS_const<T>(tv, list);

  std::vector<Teuchos::RCP<const T>> leaves;
  for (auto it : list) {
    if (it->Data() != Teuchos::null) { leaves.emplace_back(it); }
  }
  return leaves;
}

template <class T>
int
getNumTreeVectorLeaves(const T& tv)
{
  return collectTreeVectorLeaves_const(tv).size();
}


} // namespace Amanzi


#endif
