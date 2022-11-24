/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

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
  for (typename T::iterator it = tv.begin(); it != tv.end(); ++it) { list.push_back(*it); }

  for (typename T::iterator it = tv.begin(); it != tv.end(); ++it) {
    recurseTreeVectorBFS<T>(**it, list);
  }
}

template <class T>
void
recurseTreeVectorBFS_const(const T& tv, std::vector<Teuchos::RCP<const T>>& list)
{
  for (typename T::const_iterator it = tv.begin(); it != tv.end(); ++it) { list.push_back(*it); }

  for (typename T::const_iterator it = tv.begin(); it != tv.end(); ++it) {
    recurseTreeVectorBFS_const<T>(**it, list);
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
  for (typename std::vector<Teuchos::RCP<T>>::iterator it = list.begin(); it != list.end(); ++it) {
    if ((*it)->Data() != Teuchos::null) { leaves.push_back(*it); }
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
  for (typename std::vector<Teuchos::RCP<const T>>::const_iterator it = list.begin();
       it != list.end();
       ++it) {
    if ((*it)->Data() != Teuchos::null) { leaves.push_back(*it); }
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
