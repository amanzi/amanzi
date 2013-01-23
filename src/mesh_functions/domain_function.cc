/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "domain_function.hh"

#include <algorithm>
#include "errors.hh"

namespace Amanzi {

/* ******************************************************************
* Calculate pairs <list of cells, function>
****************************************************************** */
void DomainFunction::Define(const std::vector<std::string>& regions,
                            const Teuchos::RCP<const Function>& f, int action)
{
  // Form the set of cell indices that belong to any of the given regions.
  // We use a std::set here to filter out any duplicate indices.
  std::set<AmanziMesh::Entity_ID> this_domain;
  for (std::vector<std::string>::const_iterator r = regions.begin(); r != regions.end(); ++r) {
    if ((*mesh_).valid_set_name(*r, AmanziMesh::CELL)) {
      AmanziMesh::Entity_ID_List cell_list;
      (*mesh_).get_set_entities(*r, AmanziMesh::CELL, AmanziMesh::OWNED, &cell_list);
      this_domain.insert(cell_list.begin(), cell_list.end());
    } else {
      Errors::Message m;
      m << "Domain Function: unknown region: \"" << r->c_str() << "\"\n";
      Exceptions::amanzi_throw(m);
    }
  }
  
  // Verify that the cells in this_domain are disjoint from the previous domains.
  for (SpecList::const_iterator s = spec_list_.begin(); s != spec_list_.end(); ++s) {
    std::set<AmanziMesh::Entity_ID> overlap;
    std::set<AmanziMesh::Entity_ID>::iterator overlap_end;
    const std::set<AmanziMesh::Entity_ID>& prev_domain = s->first;
    std::set_intersection(prev_domain.begin(), prev_domain.end(),
                          this_domain.begin(), this_domain.end(),
                          std::inserter(overlap, overlap.end()));
    if (overlap.size() != 0) {
      Errors::Message m;
      m << "domain Fnction: conflicting definition\n";
      Exceptions::amanzi_throw(m);
    }
  }
  
  // Add this specification to the list.
  spec_list_.push_back(std::make_pair(this_domain, f));
  actions_.push_back(action);
  
  // Go ahead and register the cells in the value_ map; this enables one to
  // get at the cells in the domain function without computing its value.
  for (Domain::const_iterator d = this_domain.begin(); d != this_domain.end(); ++d) value_[*d];
}


/* ******************************************************************
* Compute and normalize the result, so far by volume
****************************************************************** */
void DomainFunction::Compute(double t)
{
  int dim = (*mesh_).space_dimension();
  double* args = new double[1+dim];
  args[0] = t;

  for (SpecList::const_iterator s = spec_list_.begin(); s != spec_list_.end(); ++s) {
    const Domain& domain = s->first;

    for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*d);
      for (int i = 0; i < dim; ++i) args[i+1] = xc[i];
      value_[*d] = (*(s->second))(args);
    }
  }
  delete [] args;
}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
void DomainFunction::ComputeDistribute(double t)
{
  int dim = (*mesh_).space_dimension();
  double* args = new double[1+dim];
  args[0] = t;

  for (SpecList::const_iterator s = spec_list_.begin(); s != spec_list_.end(); ++s) {
    const Domain& domain = s->first;

    double domain_volume = 0.0;
    for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
      domain_volume += mesh_->cell_volume(*d);
    }

    double volume_tmp = domain_volume;
    mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

    for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*d);
      for (int i = 0; i < dim; ++i) args[i+1] = xc[i];
      value_[*d] = (*(s->second))(args) / domain_volume;
    }
  }
  delete [] args;
}


/* ******************************************************************
* Compute and distribute the result by specified weight if an action
* is set on. Otherwise, weight could be a null pointer.
****************************************************************** */
void DomainFunction::ComputeDistribute(double t, double* weight)
{
  int dim = (*mesh_).space_dimension();
  double* args = new double[1+dim];
  args[0] = t;

  int nspec = spec_list_.size();
  for (int i = 0; i < nspec; i++) {
    Spec& s = spec_list_[i];
    const Domain& domain = s.first;

    int action = actions_[i];
    double domain_volume = 0.0;
    if (action == Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_VOLUME) {
      for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
        domain_volume += mesh_->cell_volume(*d);
      }

      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*d);
        for (int i = 0; i < dim; ++i) args[i+1] = xc[i];
        value_[*d] = (*(s.second))(args) / domain_volume;
      }

    } else if (action == Amanzi::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
        domain_volume += weight[*d] * mesh_->cell_volume(*d);
      }

      double volume_tmp = domain_volume;
      mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

      for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*d);
        for (int i = 0; i < dim; ++i) args[i+1] = xc[i];
        value_[*d] = (*(s.second))(args) * weight[*d] / domain_volume;
      }

    } else {
      for (Domain::const_iterator d = domain.begin(); d != domain.end(); ++d) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*d);
        for (int i = 0; i < dim; ++i) args[i+1] = xc[i];
        value_[*d] = (*(s.second))(args);
      }
    }
  }
  delete [] args;
}


/* ******************************************************************
* Return all specified actions. 
****************************************************************** */
int DomainFunction::ActionsList()
{
  int list = 0;
  int nspec = spec_list_.size();
  for (int s = 0; s < nspec; s++) {
    int action = actions_[s];
    list |= action;
  }
  return list;
}

}  // namespace Amanzi
