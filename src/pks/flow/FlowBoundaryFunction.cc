/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
           Ethan Coon (version 2)
*/

#include "FlowBoundaryFunction.hh"

namespace Amanzi {
namespace Flow {

/* ****************************************************************
* Populates internal array with function values shifted by a
* face-based vector shift.
**************************************************************** */
void FlowBoundaryFunction::ComputeShift(double time, double* shift)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  if (unique_specs_.size() == 0) return;

  // create the input tuple
  int dim = mesh_->space_dimension();
  std::vector<double> args(1+dim);
  args[0] = time;

  // Loop over all FACE specs and evaluate the function at all IDs in the
  // list.
  for (UniqueSpecList::const_iterator it = unique_specs_[AmanziMesh::FACE]->begin();
       it != unique_specs_[AmanziMesh::FACE]->end(); ++it) {
    // Here we could specialize on the argument signature of the function:
    // time-independent functions need only be evaluated at each face on the
    // first call; space-independent functions need only be evaluated once per
    // call and the value used for all faces; etc. Right now we just assume
    // the most general case.
    Teuchos::RCP<MeshIDs> ids = (*it)->second;
    for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
      const AmanziGeometry::Point& xc = mesh_->face_centroid(*id);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      // Careful tracing of the typedefs is required here: it->first
      //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*it)->first->second)(args)[0] + shift[*id];
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi
