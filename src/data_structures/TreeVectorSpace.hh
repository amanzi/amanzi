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

//!
#ifndef AMANZI_TREEVECTOR_SPACE_HH_
#define AMANZI_TREEVECTOR_SPACE_HH_

#include <vector>
#include "Teuchos_RCP.hpp"

#include "AmanziComm.hh"
#include "Iterators.hh"
#include "Mesh.hh"
#include "CompositeVectorSpace.hh"

namespace Amanzi {

class CompositeSpace;

class TreeVectorSpace {
 public:
  // Constructor
  TreeVectorSpace() : comm_(Amanzi::getDefaultComm()){};
  TreeVectorSpace(const Comm_ptr_type& comm) : comm_(comm){};
  explicit TreeVectorSpace(const Teuchos::RCP<const CompositeSpace>& cvfac);
  TreeVectorSpace(const TreeVectorSpace& other);

  // Comparators
  bool isSameAs(const TreeVectorSpace& other) const;
  bool locallySameAs(const TreeVectorSpace& other) const;
  bool isSubsetOf(const TreeVectorSpace& other) const;

  // Set/get data space
  Teuchos::RCP<const CompositeSpace> getData() const { return data_; }
  void setData(const Teuchos::RCP<const CompositeSpace>& data) { data_ = data; }

  // Access to ANY communicator (this may be ill-posed!)
  Comm_ptr_type getComm() const { return comm_; }

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
  Teuchos::RCP<const TreeVectorSpace> getSubVector(int index) const;

  // Add a sub-vector as a child of this node.
  void PushBack(const Teuchos::RCP<TreeVectorSpace>& subvec);

  // I/O
  void print(std::ostream& os) const;

 private:
  // private and unimplemented
  TreeVectorSpace& operator=(const TreeVectorSpace&);

 private:
  Teuchos::RCP<const CompositeSpace> data_;
  std::vector<Teuchos::RCP<TreeVectorSpace>> subvecs_;
  Comm_ptr_type comm_;
};


// non-member functions
Teuchos::RCP<TreeVectorSpace>
createTVSwithOneLeaf(const CompositeVectorSpace& cvs);

} // namespace Amanzi

#endif
