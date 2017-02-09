/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_TRANSPORT_SOURCE_FUNCTION_ALQUIMIA_HH_
#define AMANZI_TRANSPORT_SOURCE_FUNCTION_ALQUIMIA_HH_

#include <vector>
#include <map>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "TabularStringFunction.hh"
#include "TransportDomainFunction.hh"

#ifdef ALQUIMIA_ENABLED
#include "Alquimia_PK.hh"
#include "ChemistryEngine.hh"

namespace Amanzi {
namespace Transport {


template <class Value_Type>
class TransportSourceFunction_Alquimia : public TransportDomainFunction<Value_Type> {
 public:
  TransportSourceFunction_Alquimia(const Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                   Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk,
                                   Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine)
    : mesh_(mesh),
      chem_pk_(chem_pk),
      chem_engine_(chem_engine)
{
   // Check arguments.
  if (chem_engine_ != Teuchos::null) {
    chem_engine_->InitState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
    chem_engine_->GetPrimarySpeciesNames(TransportDomainFunction<Value_Type>::tcc_names_);
  } else {
    Errors::Message msg;
    msg << "Geochemistry is off, but a geochemical condition was requested.";
    Exceptions::amanzi_throw(msg); 
  }

  // Get the regions assigned to this geochemical condition. We do not
  // check for region overlaps here, since a better way is to derive from 
  // the generic mesh function.
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();
  std::vector<double> times = plist.get<Teuchos::Array<double> >("times").toVector();
  std::vector<std::string> conditions = plist.get<Teuchos::Array<std::string> >("geochemical conditions").toVector();

  // Function of geochemical conditions and the associates regions.
  f_ = Teuchos::rcp(new TabularStringFunction(times, conditions));
  Init_(regions);

}

  ~TransportSourceFunction_Alquimia(){
    chem_engine_->FreeState(alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);
  }
  
  void Compute(double t_old, double t_new){
    std::string cond_name = (*f_)(t_new);

    // Loop over sides and evaluate values.
    for (auto it = begin(); it != end(); ++it) {
      int cell = it->first; 

      // Dump the contents of the chemistry state into our Alquimia containers.
      chem_pk_->CopyToAlquimia(cell, alq_mat_props_, alq_state_, alq_aux_data_);

      // Enforce the condition.
      chem_engine_->EnforceCondition(cond_name, t_new, alq_mat_props_, alq_state_, alq_aux_data_, alq_aux_output_);

      // Move the concentrations into place.
      WhetStone::DenseVector& values = it->second;
      for (int i = 0; i < values.NumRows(); i++) {
        values(i) = alq_state_.total_mobile.data[i] / TransportDomainFunction<Value_Type>::domain_volume_;
      }
    }
  }

  // require by the case class
  virtual std::string name() const { return "volume"; } 

protected:

  std::map<int, Value_Type> value_;

  // string function of geochemical conditions
  Teuchos::RCP<TabularStringFunction> f_;

  // Chemistry state and engine.
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> chem_pk_;
  Teuchos::RCP<AmanziChemistry::ChemistryEngine> chem_engine_;

  // Containers for interacting with the chemistry engine.
  AlquimiaState alq_state_;
  AlquimiaMaterialProperties alq_mat_props_;
  AlquimiaAuxiliaryData alq_aux_data_;
  AlquimiaAuxiliaryOutputData alq_aux_output_;


 private:

  void Init_(const std::vector<std::string> &regions){
    for (int i = 0; i < regions.size(); ++i) {
      AmanziMesh::Entity_ID_List block;
      mesh_->get_set_entities(regions[i], AmanziMesh::CELL, AmanziMesh::USED, &block);
      int nblock = block.size();

      // Now get the cells that are attached to these faces.
      for (int n = 0; n < nblock; ++n) {
        int c = block[n];
        value_[c] = WhetStone::DenseVector(chem_engine_->NumPrimarySpecies());
      }
    }
  };


 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // iterator methods
  typename std::map<int, Value_Type>::iterator begin() { return value_.begin(); }
  typename std::map<int, Value_Type>::iterator end() { return value_.end(); }
  typename std::map<int, Value_Type>::size_type size() { return value_.size(); }


};

}  // namespace Transport
}  // namespace Amanzi

#endif


#endif
