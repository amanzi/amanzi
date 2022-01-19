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

Interface for CompositeVector, an implementation of a block Tpetra
vector on multiple maps of a single mesh and knows how to communicate
itself.

CompositeVectors are a collection of vectors defined on a common mesh and
communicator.  Each vector, or component, has a name (used as a key), a mesh
Entity_kind (CELL, FACE, NODE, or BOUNDARY_FACE), and a number of degrees of
freedom (dofs).  This, along with the Map_type provided from the mesh on a
given Entity_kind, is enough to create a Vector

Note that construction of the CompositeVector does not allocate the
Tpetra_Vector.  CreateData() must be called before usage.

This vector provides the duck-type interface Vec and may be used with time
integrators/nonlinear solvers.

DOCUMENT VANDELAY HERE! FIX ME --etc

*/

#ifndef AMANZI_COMPOSITEVECTOR_DECL_HH_
#define AMANZI_COMPOSITEVECTOR_DECL_HH_

#include <vector>
#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"
#include "AmanziVector.hh"

#include "dbc.hh"
#include "CompositeSpace.hh"
#include "BlockVector.hh"
#include "DataStructuresHelpers.hh"

namespace Amanzi {

//
// Class interface
//
template <typename Scalar>
class CompositeVector_ : public BlockVector<Scalar> {
 public:
  // -- Constructors --
  // Constructor from a CompositeSpace (which is like a Map).
  CompositeVector_(const Teuchos::RCP<const CompositeSpace>& space,
                   InitMode mode = InitMode::ZERO);
  // CompositeVector_(const Teuchos::RCP<const CompositeSpace>& space, bool
  // ghosted);

  // Copy constructor.
  CompositeVector_(const CompositeVector_<Scalar>& other,
                   Teuchos::DataAccess access = Teuchos::DataAccess::Copy,
                   InitMode mode = InitMode::COPY);
  // CompositeVector_(const CompositeVector_<Scalar>& other,
  //                  bool ghosted,
  //                  InitMode mode=INIT_MODE_COPY);

  // Assignment operator.
  CompositeVector_<Scalar>& operator=(const CompositeVector_<Scalar>& other);
  void assign(const CompositeVector_<Scalar>& other) { *this = other; }

  // -- Accessors to meta-data --
  // Space/VectorSpace/Map accessor.
  Teuchos::RCP<const CompositeSpace> getMap() const { return cvs_; }
  Teuchos::RCP<const AmanziMesh::Mesh> Mesh() const { return getMap()->Mesh(); }

 protected:
  virtual cMultiVector_ptr_type_<Scalar>
  GetComponent_(const std::string& name, bool ghosted = false) const override;
  virtual MultiVector_ptr_type_<Scalar>
  GetComponent_(const std::string& name, bool ghosted = false) override;



  // The Vandelay is an Importer/Exporter which allows face unknowns
  // to be spoofed as boundary face unknowns.
  void CreateVandelay_() const;
  void ApplyVandelay_() const;

 protected:
  Teuchos::RCP<const CompositeSpace> cvs_;

  // importer and vector for boundary data
  mutable MultiVector_ptr_type_<Scalar> vandelay_vector_;
};

} // namespace Amanzi


#endif
