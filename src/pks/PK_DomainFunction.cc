/*
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include "errors.hh"

#include "Domain.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {

/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void PK_DomainFunction::Define(const std::vector<std::string>& regions,
                               const Teuchos::RCP<const MultiFunction>& f,
                               int action, int submodel)
{
  // Create the domain
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  
  // add the spec
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
  submodel_.push_back(submodel);
}


/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void PK_DomainFunction::Define(const std::string& region,
                               const Teuchos::RCP<const MultiFunction>& f,
                               int action, int submodel) 
{
  RegionList regions(1, region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
  actions_.push_back(action);
  submodel_.push_back(submodel);
}


/* ******************************************************************
* Evaluate function
****************************************************************** */
void PK_DomainFunction::Compute(double t0, double t1,
                                Teuchos::RCP<const Epetra_Vector> weight)
{
  int type = CollectActionsList();

  if (type & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    ASSERT(weight != Teuchos::null);
    ComputeIntegral_(t0, t1, weight->Values()); 
  } else if (type & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {
    ComputeIntegral_(t0, t1);
  } else {
    ComputeDensity_(t0, t1);
  }
}


/* ******************************************************************
* Compute and normalize the result, so far by volume
****************************************************************** */
void PK_DomainFunction::ComputeDensity_(double t0, double t1)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }
  
  if (unique_specs_.size() == 0) return;

  double dt = t1 - t0;
  if (dt > 0.0) dt = 1.0 / dt;
 
  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  int n(0);
  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::CELL]->begin();
       uspec != unique_specs_[AmanziMesh::CELL]->end(); ++uspec) {

    args[0] = t1;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    for (MeshIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      value_[*id] = (*(*uspec)->first->second)(args)[0];
    }
   
    if (submodel_[n] == CommonDefs::DOMAIN_FUNCTION_SUBMODEL_INTEGRAL) {
      args[0] = t0;
      for (MeshIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*id] -= (*(*uspec)->first->second)(args)[0];
        value_[*id] *= dt;
       }
    }
    n++;
  }
}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
void PK_DomainFunction::ComputeIntegral_(double t0, double t1)
{
  // lazily generate space for the values
  if (!finalized_) {
    Finalize();
  }

  double dt = t1 - t0;
  if (dt > 0.0) dt = 1.0 / dt;

   // create the input tuple (time + space)
  int dim = (*mesh_).space_dimension();
  std::vector<double> args(1 + dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  int n(0);  
  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::CELL]->begin();
       uspec != unique_specs_[AmanziMesh::CELL]->end(); ++uspec) {

    double domain_volume = 0.0;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    // calculate physical volume of region.
    for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
      if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
    }
    double volume_tmp = domain_volume;
    mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

    args[0] = t1;
    for (MeshIDs::const_iterator id = ids->begin(); id!=ids->end(); ++id) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      value_[*id] = (*(*uspec)->first->second)(args)[0] / domain_volume;
    }

    if (submodel_[n] == CommonDefs::DOMAIN_FUNCTION_SUBMODEL_INTEGRAL) {
      args[0] = t0;
      for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*id] -= (*(*uspec)->first->second)(args)[0] / domain_volume;
        value_[*id] *= dt;
      }
    }

    n++;
  }
}


/* ******************************************************************
* Compute and distribute the result by specified weight if an action
* is set on. Otherwise, weight could be a null pointer.
****************************************************************** */
void PK_DomainFunction::ComputeIntegral_(double t0, double t1, double* weight)
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

  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::CELL]->begin();
       uspec != unique_specs_[AmanziMesh::CELL]->end(); ++uspec) {

    int action = actions_[n];
    int submodel = submodel_[n];

    double domain_volume = 0.0;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    if (action == CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {
      for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
        if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id);
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      args[0] = t1;
      for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
        value_[*id] = (*(*uspec)->first->second)(args)[0] / domain_volume;
      }      

      if (submodel == CommonDefs::DOMAIN_FUNCTION_SUBMODEL_INTEGRAL) {
        args[0] = t0;
        for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*id] -= (*(*uspec)->first->second)(args)[0] / domain_volume;
          value_[*id] *= dt;
        }
      }
    }
    else if (action == CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
        if (*id < ncells_owned) domain_volume += mesh_->cell_volume(*id) * weight[*id];
      }
      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      args[0] = t1;
      for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*id] = (*(*uspec)->first->second)(args)[0] * weight[*id] / domain_volume;
      }      

      if (submodel == CommonDefs::DOMAIN_FUNCTION_SUBMODEL_INTEGRAL) {
        args[0] = t0;
        for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*id] -= (*(*uspec)->first->second)(args)[0] * weight[*id] / domain_volume;
          value_[*id] *= dt;
        }
      }
    }
    else {
      args[0] = t1;
      for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        value_[*id] = (*(*uspec)->first->second)(args)[0];
      }      

      if (submodel == CommonDefs::DOMAIN_FUNCTION_SUBMODEL_INTEGRAL) {
        args[0] = t0;
        for (MeshIDs::const_iterator id = ids->begin(); id != ids->end(); ++id) {
          const AmanziGeometry::Point& xc = mesh_->cell_centroid(*id);
          for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
          value_[*id] -= (*(*uspec)->first->second)(args)[0];
          value_[*id] *= dt;
        }
      }
    }

    n++;
  }
}


/* ******************************************************************
* Return all specified actions. 
****************************************************************** */
int PK_DomainFunction::CollectActionsList()
{
  int list(0);
  int nspec = spec_list_.size();
  for (int s = 0; s < nspec; s++) list |= actions_[s];
  return list;
}

}  // namespace Amanzi

