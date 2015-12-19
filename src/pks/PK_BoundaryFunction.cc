/*
  This is the process kernel component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "PK_BoundaryFunction.hh"

namespace Amanzi {

/* ******************************************************************
* TBW.
****************************************************************** */
void PK_BoundaryFunction::Define(const std::vector<std::string>& regions,
                                 const Teuchos::RCP<const MultiFunction>& f, 
                                 int method) 
{
  // Create the domain
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::FACE));

  // add the spec
  AddSpec(Teuchos::rcp(new Spec(domain, f)));

  for (std::vector<std::string>::const_iterator r = regions.begin(); r != regions.end(); ++r) {
    if (method != CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE) {
      CommonDefs::Action action(*r, method);
      actions_.push_back(action); 
    }
  }
}


/* ******************************************************************
* TBW.
****************************************************************** */
void PK_BoundaryFunction::Define(std::string& region,
                                 const Teuchos::RCP<const MultiFunction>& f,
                                 int method) 
{
  RegionList regions(1, region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::FACE));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));

  if (method != CommonDefs::BOUNDARY_FUNCTION_ACTION_NONE) {
    CommonDefs::Action action(region, method);
    actions_.push_back(action);   
  }
}


/* ******************************************************************
* Evaluate values at the given time.
****************************************************************** */
void PK_BoundaryFunction::Compute(double time) {
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  if (unique_specs_.size() == 0) return;

  // create the input tuple
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);
  args[0] = time;

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
      // is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*it)->first->second)(args)[0];
    }
  }
}


/* ****************************************************************
* Allocates data for mesh function by touching them.
**************************************************************** */
void PK_BoundaryFunction::Finalize() {
  finalized_ = true;
  if (unique_specs_.size() == 0) return;

  // Create the map of values, for now just setting up memory.
  for (UniqueSpecList::const_iterator it = unique_specs_[AmanziMesh::FACE]->begin();
       it != unique_specs_[AmanziMesh::FACE]->end(); ++it) {
    Teuchos::RCP<MeshIDs> ids = (*it)->second;
    for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
      value_[*id];
    };
  }

  int my_size = size();
  mesh_->get_comm()->SumAll(&my_size, &global_size_, 1);
}

}  // namespace Amanzi
