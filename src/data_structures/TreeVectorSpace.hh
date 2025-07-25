/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Space for a TreeVector on an Amanzi mesh.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_TREEVECTOR_SPACE_HH_
#define AMANZI_TREEVECTOR_SPACE_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "Mesh.hh"

namespace Amanzi {

class CompositeVectorSpace;

class TreeVectorSpace {
 public:
  // Constructor
  TreeVectorSpace()
    : comm_(Amanzi::getDefaultComm()) {};
  TreeVectorSpace(const Comm_ptr_type& comm)
    : comm_(comm) {};
  TreeVectorSpace(const Teuchos::RCP<const CompositeVectorSpace>& cvfac)
    : data_(cvfac), comm_(cvfac->Comm())
  {}
  TreeVectorSpace(const TreeVectorSpace& other);

  // Comparators
  bool SameAs(const TreeVectorSpace& other) const;
  bool SubsetOf(const TreeVectorSpace& other) const;

  // Set/get data space
  Teuchos::RCP<const CompositeVectorSpace> Data() const { return data_; }
  void SetData(const Teuchos::RCP<const CompositeVectorSpace>& data) { data_ = data; }

  // Access to ANY communicator (this may be ill-posed!)
  Comm_ptr_type Comm() const { return comm_; }

  // Access to SubVectors
  using iterator = std::vector<Teuchos::RCP<TreeVectorSpace>>::iterator;
  using const_iterator = Teuchos::RCP<const TreeVectorSpace> const* const;

  // this is a very poor-man's hacky replacement for an iterator_adaptor
  // which aims to make const_iterators iterate over pointers to const objects.
  const_iterator begin() const
  {
    return reinterpret_cast<const Teuchos::RCP<const TreeVectorSpace>*>(&*subvecs_.begin());
  }
  const_iterator end() const
  {
    return reinterpret_cast<const Teuchos::RCP<const TreeVectorSpace>*>(&*subvecs_.end());
  }

  iterator begin() { return iterator(subvecs_.begin()); }
  iterator end() { return iterator(subvecs_.end()); }
  size_t size() const { return subvecs_.size(); }

  // Get a pointer to the sub-vector by index
  Teuchos::RCP<const TreeVectorSpace> SubVector(int index) const;

  // Add a sub-vector as a child of this node.
  void PushBack(const Teuchos::RCP<TreeVectorSpace>& subvec);

  // I/O
  void Print(std::ostream& os) const;

 private:
  // private and unimplemented
  TreeVectorSpace& operator=(const TreeVectorSpace&);

 private:
  Teuchos::RCP<const CompositeVectorSpace> data_;
  std::vector<Teuchos::RCP<TreeVectorSpace>> subvecs_;
  Comm_ptr_type comm_;
};


// non-member functions
inline Teuchos::RCP<TreeVectorSpace>
CreateTVSwithOneLeaf(const CompositeVectorSpace& cvs)
{
  auto tvs = Teuchos::rcp(new TreeVectorSpace());
  tvs->SetData(Teuchos::rcpFromRef(cvs));
  return tvs;
}

} // namespace Amanzi

#endif
