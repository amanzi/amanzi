/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_CONTAINER_HH_
#define AMANZI_CONTAINER_HH_

#include "boost/iterator_adaptors.hpp"
#include "Teuchos_RCP.hpp"

namespace Amanzi {
namespace Utils {

template <class Container_t, class Contained_t>
struct iterator : boost::iterator_adaptor<iterator<Container_t, Contained_t>,
                                          typename Container_t::iterator> {
  iterator() {}
  iterator(typename Container_t::iterator i)
    : boost::iterator_adaptor<iterator<Container_t, Contained_t>,
                              typename Container_t::iterator>(i)
  {}
};

template <typename Container_t, typename Contained_t>
struct const_iterator
  : boost::iterator_adaptor<const_iterator<Container_t, Contained_t>,
                            typename Container_t::const_iterator,
                            Teuchos::RCP<const Contained_t>, boost::use_default,
                            Teuchos::RCP<const Contained_t>> {
  const_iterator() {}
  const_iterator(iterator<Container_t, Contained_t> i)
    : boost::iterator_adaptor<
        const_iterator<Container_t, Contained_t>,
        typename Container_t::const_iterator, Teuchos::RCP<const Contained_t>,
        boost::use_default, Teuchos::RCP<const Contained_t>>(i.base())
  {}
  const_iterator(typename Container_t::const_iterator i)
    : boost::iterator_adaptor<
        const_iterator<Container_t, Contained_t>,
        typename Container_t::const_iterator, Teuchos::RCP<const Contained_t>,
        boost::use_default, Teuchos::RCP<const Contained_t>>(i)
  {}
};

} // namespace Utils
} // namespace Amanzi

#endif
