/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#include "StateDefs.hh"
#include "TransportBoundaryFunction_Coupler.hh"


namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
  void TransportBoundaryFunction_Coupler::Define(const std::vector<std::string> &regions,
                                                 Teuchos::ParameterList& spec)
{
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::FACE));
  spec_list_ = spec;
  Finalize_();
}


/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
  void TransportBoundaryFunction_Coupler::Define( std::string region,
                                                  Teuchos::ParameterList& spec)
{
  RegionList regions(1,region);
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::FACE));
  spec_list_ = spec;

  Finalize_();
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void TransportBoundaryFunction_Coupler::Compute(double time) {

  // create the input tuple
  // int dim = mesh_->space_dimension();
  // std::vector<double> args(1+dim);
  // args[0] = time;

  Key surface_tcc_key;
  if (spec_list_.isParameter("coupling_key")){
    surface_tcc_key = spec_list_.get<std::string>("coupling_key", "surface-total_component_concentration");
  }


  if (!S_->HasField(surface_tcc_key)) return;

  // Loop over side set specs and evaluate the function at all faces 
  // in the side set list.
  Teuchos::RCP<const Field> fld = S_->GetField(surface_tcc_key);

  if (!finalize_) Finalize_();

  Teuchos::RCP<const Field> copy;// = fld->GetCopy("subcycling");
  Teuchos::RCP<const CompositeVector> surf = fld->GetFieldData();
  const Epetra_MultiVector& surf_c = *surf->ViewComponent("cell",false);


  int i=0;
  for (unsigned int sc=0; sc!=surf_c.MyLength(); ++sc) {
    values_[i][0] = surf_c[0][sc];            
    i++;
  }

  

}


/* ******************************************************************
* Generate space for data (face ids and values).
****************************************************************** */
void TransportBoundaryFunction_Coupler::Finalize_() 
{

  std::vector<double> v;
  v.push_back(0.0);

  Key surface_tcc_key;
  if (spec_list_.isParameter("coupling_key")){
    surface_tcc_key = spec_list_.get<std::string>("coupling_key", "surface-total_component_concentration");
  }

  if (!S_->HasField(surface_tcc_key)) return;

  Teuchos::RCP<const CompositeVector> surf = S_->GetFieldData(surface_tcc_key);
  const Epetra_MultiVector& surf_c = *surf->ViewComponent("cell",false);


  for (unsigned int sc=0; sc!=surf_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
      surf->Mesh()->entity_get_parent(AmanziMesh::CELL, sc);

    faces_.push_back(f);
    values_.push_back(v);
  }
  finalize_ = true;


}

}  // namespace Transport
}  // namespace Amanzi

