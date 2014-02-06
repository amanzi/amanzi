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
namespace AmanziTransport {

/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void TransportDomainFunction::Define(const std::vector<std::string>& regions,
				const Teuchos::RCP<const MultiFunction>& f,
				int action,
                                const std::string& name)
{
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
  names_.push_back(name);
}


/* ******************************************************************
* TBW.
****************************************************************** */
void TransportDomainFunction::Define(const std::string region,
                                     const Teuchos::RCP<const MultiFunction>& f,
                                     int action,
                                     const std::string& name)
{
  RegionList regions(1,region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
  names_.push_back(name);
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
  double *args = new double[1+dim];
  double *xargs = args+1;
  args[0] = time;
 
  for (SpecAndIDsList::const_iterator
           spec_and_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
           spec_and_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {
    
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;
    for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
      for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
      // Careful tracing of the typedefs is required here: spec_and_ids->first
      //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*spec_and_ids)->first->second)(args)[0];
    }
  }

  delete [] args;
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

  int dim = (*mesh_).space_dimension();
  double* args = new double[1+dim];
  double *xargs = args+1;
  args[0] = time;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (SpecAndIDsList::const_iterator
           spec_and_ids=specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_and_ids!=specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;

    for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
    }
    double volume_tmp = domain_volume;
    mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

    for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
      for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
      // Careful tracing of the typedefs is required here: spec_and_ids->first
      //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
      value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] / domain_volume;
    }
  }

  delete [] args;
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

  int dim = (*mesh_).space_dimension();
  double* args = new double[1+dim];
  double *xargs = args+1;
  args[0] = t;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int i_action(0);

  for (SpecAndIDsList::const_iterator
           spec_and_ids=specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_and_ids!=specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {

    int action = actions_[i_action];
    ++i_action;

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;

    if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {

      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
	// Careful tracing of the typedefs is required here: spec_and_ids->first
	//  is a RCP<Spec>, and the Spec's second is an RCP to the function.
	value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] / domain_volume;
      }      
    } else if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {

      for (SpecIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
	if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	for (int i=0; i!=dim; ++i) xargs[i] = xc[i];
	// Careful tracing of the typedefs is required here: spec_and_ids->first
	//  is a RCP<Spec>, and the Spec's second is an RCP to the function.
	value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] * weight[*id] / domain_volume;
      }      
    } else {
      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
	// Careful tracing of the typedefs is required here: spec_and_ids->first
	//  is a RCP<Spec>, and the Spec's second is an RCP to the function.
	value_[*id] = (*(*spec_and_ids)->first->second)(args)[0];
      }      
    }
  }

  delete [] args;
}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
void TransportDomainFunction::ComputeDistributeMultiValue(double t, const std::string& name)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  int dim = (*mesh_).space_dimension();
  double* args = new double[1+dim];
  double *xargs = args+1;
  args[0] = t;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int i_name(0);
  value_.clear();

  for (SpecAndIDsList::const_iterator
           spec_and_ids=specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_and_ids!=specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {
    
    if (names_[i_name] == name) {
      double domain_volume = 0.0;
      Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;

      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
        if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
      }

      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
        AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
        for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
        // Careful tracing of the typedefs is required here: spec_and_ids->first
        //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
        value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] / domain_volume;
      }
    }
    i_name++;
  }

  delete [] args;
}


/* ******************************************************************
 * * Compute and distribute the result by volume.
 * ****************************************************************** */
void TransportDomainFunction::ComputeDistributeMultiValue(
    double t, const std::string& name, double* weight)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  int dim = (*mesh_).space_dimension();
  double* args = new double[1+dim];
  double *xargs = args+1;
  args[0] = t;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int i_action(0), i_name(0);
  value_.clear();

  for (SpecAndIDsList::const_iterator
           spec_and_ids=specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_and_ids!=specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_and_ids) {

    int action = actions_[i_action];
    ++i_action;

    if (names_[i_name] == name) {
      double domain_volume = 0.0;
      Teuchos::RCP<SpecIDs> ids = (*spec_and_ids)->second;

      if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {
        for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
          if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
        }

        double volume_tmp = domain_volume;
        mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

        for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
            AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	  for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
	  // Careful tracing of the typedefs is required here: spec_and_ids->first
	  //  is a RCP<Spec>, and the Spec's second is an RCP to the function.
	  value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] / domain_volume;
        }      
      } else if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
        for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
	  if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
        }

        double volume_tmp = domain_volume;
        mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

        for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
            AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	  for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
	  value_[*id] = (*(*spec_and_ids)->first->second)(args)[0] * weight[*id] / domain_volume;
        }      
      } else {
        for (SpecIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
            AmanziGeometry::Point xc = mesh_->cell_centroid(*id);
	  for (int i = 0; i != dim; ++i) xargs[i] = xc[i];
	  value_[*id] = (*(*spec_and_ids)->first->second)(args)[0];
        }        
      }
      i_name++;
    }
  }

  delete [] args;
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


}  // namespace Functions
}  // namespace Amanzi
