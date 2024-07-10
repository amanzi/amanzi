/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Factory for Output objects.
/*

*/

#ifndef AMANZI_OUTPUT_FACTORY_HH_
#define AMANZI_OUTPUT_FACTORY_HH_

#include <memory>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziTypes.hh"
#include "Mesh.hh"
#include "Output.hh"

namespace Amanzi {
namespace OutputFactory {

// creates a formatter for the directory name of a multi-file checkpoint/vis
Output::FilenameFormatter
createDirectoryFormatter(Teuchos::ParameterList& plist);

std::unique_ptr<Output>
createForVis(Teuchos::ParameterList& plist, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

std::unique_ptr<Output>
createForCheckpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm);

} // namespace OutputFactory
} // namespace Amanzi

#endif
