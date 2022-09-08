/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky (dasvyatv@lanl.gov)
*/

#include "TransportSourceFunction_Concentration.hh"

#ifdef ALQUIMIA_ENABLED

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* Constructor of BCs for Alquimia.
****************************************************************** */
TransportSourceFunction_Concentration::TransportSourceFunction_Concentration(
    const Teuchos::ParameterList& plist,
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : mesh_(mesh)
{

  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();
  components_ = plist.get<Teuchos::Array<std::string> >("regions").toVector();
 
  Init_(regions);
}


/* ******************************************************************
* Delegating destructor.
****************************************************************** */
TransportSourceFunction_Concentration::~TransportSourceFunction_Concentration()
{}

/* ******************************************************************
* Internal subroutine that defines a source function.
****************************************************************** */
void TransportSourceFunction_Concentration::Init_(const std::vector<std::string>& regions)
{
  for (int i = 0; i < regions.size(); ++i) {
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(regions[i], AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL, &block);
    int nblock = block.size();

    // Now get the cells that are attached to these faces.
    for (int n = 0; n < nblock; ++n) {
      int c = block[n];
      value_[c].resize(components_.size());
    }
  }
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void TransportSourceFunction_Concentration::Compute(double t_old, double t_new) 
{
  // Loop over sides and evaluate values.
  for (auto it = begin(); it != end(); ++it) {
    int cell = it->first; 


}

}  // namespace Transport
}  // namespace Amanzi

#endif
