/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_TREEVECTOR_SPACE_HH_
#define AMANZI_TREEVECTOR_SPACE_HH_

#include <vector>
#include "Teuchos_RCP.hpp"

#include "AmanziComm.hh"
#include "Iterators.hh"
#include "Mesh.hh"

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
  bool SameAs(const TreeVectorSpace& other) const;
  bool LocallySameAs(const TreeVectorSpace& other) const;
  bool SubsetOf(const TreeVectorSpace& other) const;

  // Set/get data space
  Teuchos::RCP<const CompositeSpace> Data() const { return data_; }
  void SetData(const Teuchos::RCP<const CompositeSpace>& data) { data_ = data; }

  // Access to ANY communicator (this may be ill-posed!)
  Comm_ptr_type Comm() const { return comm_; }
  Comm_ptr_type getComm() const { return comm_; }

  // Access to SubVectors
  typedef std::vector<Teuchos::RCP<TreeVectorSpace>> SubVectorsContainer;
  typedef Utils::iterator<SubVectorsContainer, TreeVectorSpace> iterator;
  typedef Utils::const_iterator<SubVectorsContainer, const TreeVectorSpace>
    const_iterator;

  const_iterator begin() const { return const_iterator(subvecs_.begin()); }
  const_iterator end() const { return const_iterator(subvecs_.end()); }

  iterator begin() { return iterator(subvecs_.begin()); }
  iterator end() { return iterator(subvecs_.end()); }
  size_t size() const { return subvecs_.size(); }

  // Get a pointer to the sub-vector by index
  Teuchos::RCP<const TreeVectorSpace> SubVector(int index) const;

  // Add a sub-vector as a child of this node.
  void PushBack(const Teuchos::RCP<TreeVectorSpace>& subvec);

 private:
  // private and unimplemented
  TreeVectorSpace& operator=(const TreeVectorSpace&);

 private:
  Teuchos::RCP<const CompositeSpace> data_;
  std::vector<Teuchos::RCP<TreeVectorSpace>> subvecs_;
  Comm_ptr_type comm_;
};

} // namespace Amanzi

#endif
