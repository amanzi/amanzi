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

This vector provides the duck-type interface Vec and may be used with time
integrators/nonlinear solvers.

DOCUMENT VANDELAY HERE! FIX ME --etc

*/

#include "AmanziTypes.hh"
#include "CompositeVector_decl.hh"
#include "CompositeVector_impl.hh"

#ifndef AMANZI_COMPOSITEVECTOR_HH_
#  define AMANZI_COMPOSITEVECTOR_HH_
namespace Amanzi {
using CompositeVector = CompositeVector_<double_type>;

void
DeriveFaceValuesFromCellValues(CompositeVector&);

} // namespace Amanzi
#endif
