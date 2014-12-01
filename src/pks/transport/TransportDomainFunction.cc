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
void TransportDomainFunction::Define(
    const std::vector<std::string>& regions,
    const Teuchos::RCP<const MultiFunction>& f,
    int action, int submodel, const std::string& name)
{
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
  submodel_.push_back(submodel);
  tcc_name_ = name;
}


/* ******************************************************************
* TBW.
****************************************************************** */
void TransportDomainFunction::Define(
    const std::string region,
    const Teuchos::RCP<const MultiFunction>& f,
    int action, int submodel, const std::string& name)
{
  RegionList regions(1,region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
  submodel_.push_back(submodel);
  tcc_name_ = name;
}


/* ******************************************************************
* Compute and normalize the result, so far by volume
****************************************************************** */
void TransportDomainFunction::Compute(double t0, double t1)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }
  
  if (specs_and_ids_.size() == 0) return;

  double dt = t1 - t0;
  if (dt > 0.0) dt = 1.0 / dt;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);
 
  int n(0);
  for (SpecAndIDsList::const_iterator
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {
    
    args[0] = t1;
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;
    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      // spec_ids->first is a RCP<Spec>, Spec's second is an RCP to the function.
      value_[*c] = (*(*spec_ids)->first->second)(args)[0];
    }

    if (submodel_[n] == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
      args[0] = t0;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*c] -= (*(*spec_ids)->first->second)(args)[0];
        value_[*c] *= dt;
      }
    }

    n++;
  }
}



/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
void TransportDomainFunction::ComputeDistribute(double t0, double t1)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  double dt = t1 - t0;
  if (dt > 0.0) dt = 1.0 / dt;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  int n(0);
  for (SpecAndIDsList::const_iterator
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;

    args[0] = t1;
    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
    }
    double volume_tmp = domain_volume;
    mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      // spec_ids->first is a RCP<Spec>, Spec's second is an RCP to the function.
      value_[*c] = (*(*spec_ids)->first->second)(args)[0] / domain_volume;
    }

    if (submodel_[n] == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
      args[0] = t0;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*c] -= (*(*spec_ids)->first->second)(args)[0] / domain_volume;
        value_[*c] *= dt;
      }
    }

    n++;
  }
}


/* ******************************************************************
* Compute and distribute the result by specified weight if an action
* is set on. Otherwise, weight could be a null pointer.
****************************************************************** */
void TransportDomainFunction::ComputeDistribute(double t0, double t1, double* weight)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  double dt = t1 - t0;
  if (dt > 0.0) dt = 1.0 / dt;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int n(0);

  for (SpecAndIDsList::const_iterator
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {

    int action = actions_[n];
    int submodel = submodel_[n];

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;

    if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      args[0] = t1;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        // spec_ids->first is a RCP<Spec>, Spec's second is an RCP to the function.
        value_[*c] = (*(*spec_ids)->first->second)(args)[0] / domain_volume;
      }      

      if (submodel == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
        args[0] = t0;
        for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*c] -= (*(*spec_ids)->first->second)(args)[0] / domain_volume;
          value_[*c] *= dt;
        }
      }
    }
    else if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      args[0] = t1;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*c] = (*(*spec_ids)->first->second)(args)[0] * weight[*c] / domain_volume;
      }      

      if (submodel == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
        args[0] = t0;
        for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*c] -= (*(*spec_ids)->first->second)(args)[0] * weight[*c] / domain_volume;
          value_[*c] *= dt;
        }
      }
    }
    else {
      args[0] = t1;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*c] = (*(*spec_ids)->first->second)(args)[0];
      }      

      if (submodel == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
        args[0] = t0;
        for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*c] -= (*(*spec_ids)->first->second)(args)[0];
          value_[*c] *= dt;
        }
      }
    }

    n++;
  }
}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
void TransportDomainFunction::ComputeDistributeMultiValue(double t0, double t1)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  double dt = t1 - t0;
  if (dt > 0.0) dt = 1.0 / dt;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  value_.clear();

  int n(0);
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

    args[0] = t1;
    for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      // spec_ids->first is a RCP<Spec>, Spec's second is an RCP to the function.
      value_[*c] = (*(*spec_ids)->first->second)(args)[0] / domain_volume;
    }

    if (submodel_[n] == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
      args[0] = t0;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*c] -= (*(*spec_ids)->first->second)(args)[0] / domain_volume;
        value_[*c] *= dt;
      }
    }

    n++;
  }
}


/* ******************************************************************
 * * Compute and distribute the result by volume.
 * ****************************************************************** */
void TransportDomainFunction::ComputeDistributeMultiValue(double t0, double t1, double* weight)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  double dt = t1 - t0;
  if (dt > 0.0) dt = 1.0 / dt;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  int n(0);
  value_.clear();

  for (SpecAndIDsList::const_iterator
       spec_ids = specs_and_ids_[AmanziMesh::CELL]->begin();
       spec_ids != specs_and_ids_[AmanziMesh::CELL]->end(); ++spec_ids) {

    int action = actions_[n];
    int submodel = submodel_[n];

    double domain_volume = 0.0;
    Teuchos::RCP<SpecIDs> ids = (*spec_ids)->second;

    if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
      }

      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      args[0] = t1;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        // spec_ids->first is a RCP<Spec>, Spec's second is an RCP to the function.
        value_[*c] = (*(*spec_ids)->first->second)(args)[0] / domain_volume;
      }      

      if (submodel == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
        args[0] = t0;
        for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*c] -= (*(*spec_ids)->first->second)(args)[0] / domain_volume;
          value_[*c] *= dt;
        }
      }
    } 
    else if (action == TransportActions::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
      }

      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      args[0] = t1;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*c] = (*(*spec_ids)->first->second)(args)[0] * weight[*c] / domain_volume;
      }      

      if (submodel == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
        args[0] = t0;
        for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*c] -= (*(*spec_ids)->first->second)(args)[0] * weight[*c] / domain_volume;
          value_[*c] *= dt;
        }
      }
    }
    else {
      args[0] = t1;
      for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*c] = (*(*spec_ids)->first->second)(args)[0];
      }        

      if (submodel == TransportActions::DOMAIN_FUNCTION_SUBMODEL_INTEGRAND) {
        args[0] = t0;
        for (SpecIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*c] -= (*(*spec_ids)->first->second)(args)[0];
          value_[*c] *= dt;
        }
      }
    }

    n++;
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
