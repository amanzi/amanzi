/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Rel perm( pc ( sat ) ).
*/

#include "Mesh_Algorithms.hh"
#include "RelPermEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
                                   const Teuchos::Ptr<State>& S,
                                   const Teuchos::RCP<WRMPartition>& wrm) :
    SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm) {
  InitializeFromPlist_(S);
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrm_(other.wrm_),
    pressure_key_(other.pressure_key_),
    patm_(other.patm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void RelPermEvaluator::InitializeFromPlist_(const Teuchos::Ptr<State>& S)
{
  // my keys is for rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("relative permeability key");
  }

  // my dependency is pressure.
  std::string domain = Keys::getDomain(my_key_);
  pressure_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(pressure_key_);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  patm_ = *S->GetScalarData("atmospheric_pressure");
  const auto& pres = S->GetFieldData(pressure_key_); 

  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    const auto& pres_c = *pres->ViewComponent(*comp);
    auto& result_c = *result->ViewComponent(*comp);

    int nids = pres_c.MyLength();
    if (*comp == "cell") {
      for (int c = 0; c != nids; ++c) {
        result_c[0][c] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm_ - pres_c[0][c]);
      }
    } else if (*comp == "boundary_face") {
      for (int bf = 0; bf != nids; ++bf) {
        int c = AmanziMesh::getBoundaryFaceInternalCell(*pres->Mesh(), bf);
        result_c[0][bf] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm_ - pres_c[0][bf]);
      }
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  patm_ = *S->GetScalarData("atmospheric_pressure");
  const auto& pres = S->GetFieldData(pressure_key_); 

  if (wrt_key == pressure_key_) {
    for (auto comp = result->begin(); comp != result->end(); ++comp) {
      const auto& pres_c = *pres->ViewComponent(*comp);
      auto& result_c = *result->ViewComponent(*comp);

      int nids = pres_c.MyLength();
      if (*comp == "cell") {
        for (int c = 0; c != nids; ++c) {
          // Negative sign indicates that dKdP = -dKdPc.
          result_c[0][c] = -wrm_->second[(*wrm_->first)[c]]->dKdPc(patm_ - pres_c[0][c]);
        }
      } else if (*comp == "boundary_face") {
        for (int bf = 0; bf != nids; ++bf) {
          int c = AmanziMesh::getBoundaryFaceInternalCell(*pres->Mesh(), bf);
          result_c[0][bf] = -wrm_->second[(*wrm_->first)[c]]->dKdPc(patm_ - pres_c[0][bf]);
        }
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace Flow
}  // namespace Amanzi
