#include <algorithm>
#include "errors.hh"

#include "Domain.hh"
#include "FlowDomainFunction.hh"

namespace Amanzi {
namespace Functions {

/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void FlowDomainFunction::Define(const std::vector<std::string>& regions,
				const Teuchos::RCP<const MultiFunction>& f,
				int action) 
{
  // Create the domain
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  
  // add the spec
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
}


/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void FlowDomainFunction::Define(std::string region,
				const Teuchos::RCP<const MultiFunction>& f,
				int action) 
{
  RegionList regions(1,region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
}


/* ******************************************************************
* Compute and normalize the result, so far by volume
****************************************************************** */
void FlowDomainFunction::Compute(double time)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }
  
  if (specs_and_ids_.size() == 0) return;

  // create the input tuple
  int dim = mesh_->space_dimension();
  std::vector<double> args(1+dim);
  args[0] = time;
 
  for (SpecAndIDsList::const_iterator
           spec_and_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_and_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {
    
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;
    for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
      for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
      // Careful tracing of the typedefs is required here: spec_and_ids->first
      //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*spec_and_ids)->first->second)(args)[0];
    }
  }

}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
void FlowDomainFunction::ComputeDistribute(double time)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  int dim = (*mesh_).space_dimension();
  std::vector<double> args(1+dim);
  args[0] = time;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (SpecAndIDsList::const_iterator
           spec_and_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
           spec_and_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;

    // calculate physical volume of region.
    for (SpecIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
      if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
    }
    double volume_tmp = domain_volume;
    mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

    for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
      for (int i=0; i!=dim; ++i) args[i+1] = xc[i];
      // Careful tracing of the typedefs is required here: spec_and_ids->first
      //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] / domain_volume;
    }

  }
}


/* ******************************************************************
* Compute and distribute the result by specified weight if an action
* is set on. Otherwise, weight could be a null pointer.
****************************************************************** */
void FlowDomainFunction::ComputeDistribute(double t, double* weight)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  int dim = (*mesh_).space_dimension();
  std::vector<double> args(1+dim);
  args[0] = t;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  int i_action(0);
  for (SpecAndIDsList::const_iterator
           spec_and_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_and_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {

    int action = actions_[i_action];
    ++i_action;

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;

    if (action == DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {

      // calculate physical volume of region.
      for (SpecIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
	if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
	// Careful tracing of the typedefs is required here: spec_and_ids->first
	// is a RCP<Spec>, and the Spec's second is an RCP to the function.
	value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] / domain_volume;
      }      
    } else if (action == DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {

      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
	// Careful tracing of the typedefs is required here: spec_and_ids->first
	// is a RCP<Spec>, and the Spec's second is an RCP to the function.
	value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] * weight[*id] / domain_volume;
      }      
    } else {
      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
	// Careful tracing of the typedefs is required here: spec_and_ids->first
	// is a RCP<Spec>, and the Spec's second is an RCP to the function.
	value_[*id] = (*(*spec_and_ids)->first->second)(args)[0];
      }      
    }
  }
}


/* ******************************************************************
* Return all specified actions. 
****************************************************************** */
int FlowDomainFunction::CollectActionsList()
{
  int list(0);
  int nspec = spec_list_.size();
  for (int s = 0; s < nspec; s++) list |= actions_[s];
  return list;
}

}  // namespace Functions
}  // namespace Amanzi

