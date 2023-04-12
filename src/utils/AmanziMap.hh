/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Includes and a few helper functions and full definitions of Maps.
/*!

  Include this instead of directly including Epetra_Map.h to ease transition
  between linear algebra packages.

*/
#pragma once

#include "Teuchos_RCP.hpp"
#include "AmanziTypes.hh"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"

