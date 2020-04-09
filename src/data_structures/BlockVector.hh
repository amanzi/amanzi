/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

/*
  Interface for BlockVector, an implementation of a slightly improved
  Epetra_MultiVector which spans multiple simplices and knows how to
  communicate itself.
*/

#include "BlockVector_decl.hh"
#include "BlockVector_impl.hh"
