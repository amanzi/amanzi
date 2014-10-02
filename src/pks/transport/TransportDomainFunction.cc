/*
This is the transport component of the Amanzi code. 

License: see $AMANZI_DIR/COPYRIGHT
Author (v1): Neil Carlson
       (v2): Ethan Coon, Markus Berndt, Konstantin Lipnikov
*/


#include <algorithm>
#include "errors.hh"

#include "Domain.hh"
#include "TransportDomainFunction.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void TransportDomainFunction::Define(const std::vector<std::string>& regions,
                                     const Teuchos::RCP<const MultiFunction>& f,
                                     int action, const std::string& name)
{
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
  tcc_name_ = name;
}


/* ******************************************************************
* TBW.
****************************************************************** */
void TransportDomainFunction::Define(const std::string region,
                                     const Teuchos::RCP<const MultiFunction>& f,
                                     int action, const std::string& name)
{
  RegionList regions(1,region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
  tcc_name_ = name;
}


/* ******************************************************************
* Compute and normalize the result, so far by volume
****************************************************************** */
void TransportDomainFunction::Compute(double time)
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
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {
    
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;
    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
      // Careful tracing of the typedefs is required here: spec_ids->first
      // is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*c] = (*(*spec_ids)->first->second)(args)[0];
    }
  }
}



/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
void TransportDomainFunction::ComputeDistribute(double time)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  int dim = mesh_->space_dimension();
  std::vector<double> args(1+dim);
  args[0] = time;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (SpecAndIDsList::const_iterator
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;

    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
    }
    double volume_tmp = domain_volume;
    mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
      // Careful tracing of the typedefs is required here: spec_ids->first
      // is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*c] = (*(*spec_ids)->first->second)(args)[0] / domain_volume;
    }
  }
}


/* ******************************************************************
* Compute and distribute the result by specified weight if an action
* is set on. Otherwise, weight could be a null pointer.
****************************************************************** */
void TransportDomainFunction::ComputeDistribute(double t, double* weight)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  int dim = mesh_->space_dimension();
  std::vector<double> args(1+dim);
  args[0] = t;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int i_action(0);

  for (SpecAndIDsList::const_iterator
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {

    int action = actions_[i_action];
    ++i_action;

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;

    if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {

      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
	if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
	for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
	// Careful tracing of the typedefs is required here: spec_ids->first
	// is a RCP<Spec>, and the Spec's second is an RCP to the function.
	value_[*c] = (*(*spec_ids)->first->second)(args)[0] / domain_volume;
      }      
    } else if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {

      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
	if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
	for (int i=0; i!=dim; ++i) args[i+1] = xc[i];
	value_[*c] = (*(*spec_ids)->first->second)(args)[0] * weight[*c] / domain_volume;
      }      
    } else {
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
	for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
	value_[*c] = (*(*spec_ids)->first->second)(args)[0];
      }      
    }
  }
}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
void TransportDomainFunction::ComputeDistributeMultiValue(double t)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  int dim = mesh_->space_dimension();
  std::vector<double> args(1+dim);
  args[0] = t;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  value_.clear();

  for (SpecAndIDsList::const_iterator
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {
    
    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;

    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
    }

    double volume_tmp = domain_volume;
    mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
      // Careful tracing of the typedefs is required here: spec_ids->first
      // is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*c] = (*(*spec_ids)->first->second)(args)[0] / domain_volume;
    }
  }
}


/* ******************************************************************
 * * Compute and distribute the result by volume.
 * ****************************************************************** */
void TransportDomainFunction::ComputeDistributeMultiValue(double t, double* weight)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  int dim = mesh_->space_dimension();
  std::vector<double> args(1+dim);
  args[0] = t;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int i_action(0);
  value_.clear();

  for (SpecAndIDsList::const_iterator
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {

    int action = actions_[i_action];
    ++i_action;

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;

    if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
      }

      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
	// Careful tracing of the typedefs is required here: spec_ids->first
	// is a RCP<Spec>, and the Spec's second is an RCP to the function.
	value_[*c] = (*(*spec_ids)->first->second)(args)[0] / domain_volume;
      }      
    } else if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
      }

      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
	  value_[*c] = (*(*spec_ids)->first->second)(args)[0] * weight[*c] / domain_volume;
        }      
    } else {
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        AmanziGeometry::Point xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i+1] = xc[i];
        value_[*c] = (*(*spec_ids)->first->second)(args)[0];
      }        
    }
  }
}


/* ******************************************************************
* Return all specified actions. 
****************************************************************** */
int TransportDomainFunction::CollectActionsList()
{
  int list(0);
  for (int i = 0; i < actions_.size(); i++) list |= actions_[i];
  return list;
}


}  // namespace Transport
}  // namespace Amanzi
