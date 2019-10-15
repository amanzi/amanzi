/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>

/*

Implementation for CompositeVector, an implementation of a block Tpetra
vector which spans multiple simplices and knows how to communicate
itself.

CompositeVectors are a collection of vectors defined on a common mesh and
communicator.  Each vector, or component, has a name (used as a key), a mesh
Entity_kind (CELL, FACE, NODE, or BOUNDARY_FACE), and a number of degrees of
freedom (dofs).  This, along with the Map_type provided from the mesh on a
given Entity_kind, is enough to create a Vector

Note that construction of the CompositeVector does not allocate the
Tpetra_Vector.  CreateData() must be called before usage.

Access using operator() is slow, and should only be used for debugging.
Prefer to use the ViewComponent() accessors.

DOCUMENT VANDELAY HERE! FIX ME --etc

*/

#ifndef AMANZI_COMPOSITEVECTOR_IMPL_HH_
#define AMANZI_COMPOSITEVECTOR_IMPL_HH_

#include "Mesh.hh"
#include "CompositeSpace.hh"

namespace Amanzi {


// Constructor
template <typename Scalar>
CompositeVector_<Scalar>::CompositeVector_(
  const Teuchos::RCP<const CompositeSpace>& map, InitMode mode)
  : BlockVector<Scalar>(map, mode), cvs_(map)
{}


// template<typename Scalar>
// CompositeVector_<Scalar>::CompositeVector_(const Teuchos::RCP<const
// CompositeSpace>& map, bool ghosted)
//     : map_(map),
//       mastervec_(Teuchos::rcp(new CompVector<Scalar>(map->getMap(false)))),
//       ghostvec_(Teuchos::rcp(new CompVector<Scalar>(map->getMap(true)))),
//       ghosted_(ghosted)
// {
//   CreateData_();
// }


template <typename Scalar>
CompositeVector_<Scalar>::CompositeVector_(
  const CompositeVector_<Scalar>& other, Teuchos::DataAccess access,
  InitMode mode)
  : BlockVector<Scalar>(other, access, mode), cvs_(other.getMap())
{}

// template<typename Scalar>
// CompositeVector_<Scalar>::CompositeVector_(const CompositeVector_<Scalar>&
// other, bool ghosted, InitMode mode)
//     : map_(other->getMap()),
//       mastervec_(Teuchos::rcp(new
//       CompVector<Scalar>(other.getMap()->getMap(false)))),
//       ghostvec_(Teuchos::rcp(new
//       CompVector<Scalar>(other.getMap()->getMap(true)))), ghosted_(ghosted)
// {
//   CreateData_();
//   InitData_(other, mode);
// }

template <typename Scalar>
CompositeVector_<Scalar>&
CompositeVector_<Scalar>::operator=(const CompositeVector_<Scalar>& other)
{
  BlockVector<Scalar>::operator=(other);
  return *this;
}


// -- Access a view of a single component's data.
// Ghosted views are simply the vector itself, while non-ghosted views are
// lazily generated.
template <typename Scalar>
cMultiVector_ptr_type_<Scalar>
CompositeVector_<Scalar>::GetComponent_(const std::string& name,
                                        bool ghosted) const
{
  if (name == std::string("boundary_face") &&
      !this->HasComponent("boundary_face") && this->HasComponent("face")) {
    ApplyVandelay_();
    return vandelay_vector_;
  }
  return BlockVector<Scalar>::GetComponent_(name, ghosted);
};


template <typename Scalar>
MultiVector_ptr_type_<Scalar>
CompositeVector_<Scalar>::GetComponent_(const std::string& name, bool ghosted)
{
  if (name == std::string("boundary_face") &&
      !this->HasComponent("boundary_face") && this->HasComponent("face")) {
    ApplyVandelay_();
    return vandelay_vector_;
  }
  return BlockVector<Scalar>::GetComponent_(name, ghosted);
};


//
// Vandelay functionality
// -----------------------------------------------------------------------------
template <typename Scalar>
void
CompositeVector_<Scalar>::CreateVandelay_() const
{
  vandelay_vector_ = Teuchos::rcp(new MultiVector_type_<Scalar>(
    Mesh()->exterior_face_map(false), this->getNumVectors("face"), false));
}


template <typename Scalar>
void
CompositeVector_<Scalar>::ApplyVandelay_() const
{
  if (vandelay_vector_ == Teuchos::null) CreateVandelay_();
  vandelay_vector_->doImport(*GetComponent_("face", false),
                             *Mesh()->exterior_face_importer(),
                             Tpetra::INSERT);
}


} // namespace Amanzi

#endif
