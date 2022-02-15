/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  The interface between state/data and the model, an EOS.
*/

#include "EOSFactory.hh"
#include "IsobaricEOSEvaluator.hh"

namespace Amanzi {
namespace AmanziEOS {

/* *******************************************************************
* Constructor takes a parameter list.
******************************************************************* */
IsobaricEOSEvaluator::IsobaricEOSEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  // Process the list for my provided field.
  std::string mode = plist_.get<std::string>("eos basis", "molar");
  if (mode == "molar") {
    mode_ = EOS_MODE_MOLAR;
  } else if (mode == "mass") {
    mode_ = EOS_MODE_MASS;
  } else if (mode == "both") {
    mode_ = EOS_MODE_BOTH;
  } else {
    AMANZI_ASSERT(0);
  }

  // my keys and tags
  Key key;
  tag_ = Tags::DEFAULT;
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    key = plist_.get<std::string>("molar density key");
    my_keys_.push_back(std::make_pair(key, tag_));
  }

  if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    key = plist_.get<std::string>("mass density key");
    my_keys_.push_back(std::make_pair(key, tag_));
  }

  // set up my dependencies
  std::string domain = Keys::getDomain(key);

  // -- temperature
  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(std::make_pair(temp_key_, tag_));

  // -- pressure
  pres_key_ = plist_.get<std::string>("pressure key", "atmospheric_pressure");

  // Construct my EOS model
  AMANZI_ASSERT(plist_.isSublist("EOS parameters"));
  EOSFactory<EOS_Density> eos_fac;
  eos_ = eos_fac.Create(plist_.sublist("EOS parameters"));
};


IsobaricEOSEvaluator::IsobaricEOSEvaluator(const IsobaricEOSEvaluator& other) :
    EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    eos_(other.eos_),
    mode_(other.mode_),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_) {};


Teuchos::RCP<Evaluator> IsobaricEOSEvaluator::Clone() const {
  return Teuchos::rcp(new IsobaricEOSEvaluator(*this));
}


void IsobaricEOSEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  auto temp = S.GetPtr<CompositeVector>(temp_key_, tag_);
  double pres = S.Get<double>(pres_key_);

  int index = 0; // index to the results list
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    auto result = results[index];
    for (auto comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp));
      Epetra_MultiVector& result_v = *(result->ViewComponent(*comp));

      int count = result->size(*comp);
      for (int id = 0; id != count; ++id) {
        result_v[0][id] = eos_->MolarDensity(temp_v[0][id], pres);
      }
    }
    index++;
  }

  if (mode_ == EOS_MODE_BOTH) {
    // calculate density from MolarDensity and molar mass.
    double M = eos_->MolarMass();
    results[1]->Update(M, *(results[0]), 0.0);
  } else if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    auto result = results[index];
    for (auto comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

      int count = result->size(*comp);
      for (int id = 0; id != count; ++id) {
        result_v[0][id] = eos_->Density(temp_v[0][id], pres);
      }
    }
  }
}


void IsobaricEOSEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results) 
{
  auto temp = S.GetPtr<CompositeVector>(temp_key_, tag_);
  double pres = S.Get<double>(pres_key_, tag_);

  if (wrt_key == temp_key_) {

    int index = 0; // index to the results list
    if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
      // evaluate DMolarDensityDT()
      auto result = results[index];
      for (auto comp = result->begin(); comp != result->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
        Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

        int count = result->size(*comp);
        for (int id=0; id!=count; ++id) {
          result_v[0][id] = eos_->DMolarDensityDT(temp_v[0][id], pres);
        }
      }
      index++;
    }

    if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
      auto result = results[index];
      for (auto comp = result->begin(); comp != result->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
        Epetra_MultiVector& result_v = *result->ViewComponent(*comp);

        int count = result->size(*comp);
        for (int id = 0; id != count; ++id) {
          result_v[0][id] = eos_->DDensityDT(temp_v[0][id], pres);
        }
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace AmanziEOS
}  // namespace Amanzi
