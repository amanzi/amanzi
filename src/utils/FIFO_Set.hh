/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

A set that is first-in, first-out ordered.  This is reasonably efficient for
small sets, but is linear on insertion, containment checks, etc as data is
stored as a vector under the hood.

*/

#pragma once

namespace Amanzi {
namespace Utils {

template <typename T>
class FIFO_Set {
  using Container_type = std::vector<T>;

 public:
  FIFO_Set() {}

  void insert(const T& entry)
  {
    auto nentries = store_.size();
    for (int i = 0; i != nentries; ++i) {
      if (entry == store_[i]) return;
    }
    store_.emplace_back(entry);
  }

  void clear() { store_.clear(); }

  using const_iterator = typename Container_type::const_iterator;

  const_iterator begin() const { return store_.begin(); }

  const_iterator end() const { return store_.end(); }

  const T& front() const { return store_.front(); }

  const T& back() const { return store_.back(); }

  std::size_t size() const { return store_.size(); }

 protected:
  Container_type store_;
};


} // namespace Utils
} // namespace Amanzi
