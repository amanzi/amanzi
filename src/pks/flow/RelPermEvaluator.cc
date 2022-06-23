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
                                   const Teuchos::RCP<WRMPartition>& wrm)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
      wrm_(wrm) {
  InitializeFromPlist_(S);
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
      wrm_(other.wrm_),
      pressure_key_(other.pressure_key_),
      patm_(other.patm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator> RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void RelPermEvaluator::InitializeFromPlist_(const Teuchos::Ptr<State>& S)
{
  // my keys is for rel perm.
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("relative permeability key"), Tags::DEFAULT));
  }

  // my dependency is pressure.
  std::string domain = Keys::getDomain(my_keys_[0].first);
  pressure_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(std::make_pair(pressure_key_, Tags::DEFAULT));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  patm_ = S.Get<double>("atmospheric_pressure");
  const auto& pres = S.Get<CompositeVector>(pressure_key_); 

  for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
    const auto& pres_c = *pres.ViewComponent(*comp);
    auto& result_c = *results[0]->ViewComponent(*comp);

    int nids = pres_c.MyLength();
    if (*comp == "cell") {
      for (int c = 0; c != nids; ++c) {
        result_c[0][c] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm_ - pres_c[0][c]);
      }
    } else if (*comp == "boundary_face") {
      for (int bf = 0; bf != nids; ++bf) {
        int c = AmanziMesh::getBoundaryFaceInternalCell(*pres.Mesh(), bf);
        result_c[0][bf] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm_ - pres_c[0][bf]);
      }
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  patm_ = S.Get<double>("atmospheric_pressure");
  const auto& pres = S.Get<CompositeVector>(pressure_key_); 

  if (wrt_key == pressure_key_) {
    for (auto comp = results[0]->begin(); comp != results[0]->end(); ++comp) {
      const auto& pres_c = *pres.ViewComponent(*comp);
      auto& result_c = *results[0]->ViewComponent(*comp);

      int nids = pres_c.MyLength();
      if (*comp == "cell") {
        for (int c = 0; c != nids; ++c) {
          // Negative sign indicates that dKdP = -dKdPc.
          result_c[0][c] = -wrm_->second[(*wrm_->first)[c]]->dKdPc(patm_ - pres_c[0][c]);
        }
      } else if (*comp == "boundary_face") {
        for (int bf = 0; bf != nids; ++bf) {
          int c = AmanziMesh::getBoundaryFaceInternalCell(*pres.Mesh(), bf);
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
