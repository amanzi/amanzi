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
    SecondaryVariablesFieldEvaluator(plist)
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

  // my keys
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    a_key_ = plist_.get<std::string>("molar density key");
    my_keys_.push_back(a_key_);
  }

  if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    a_key_ = plist_.get<std::string>("mass density key");
    my_keys_.push_back(a_key_);
  }

  // Set up my dependencies.
  std::size_t end = a_key_.find_first_of("_");
  std::string domain_name = a_key_.substr(0,end);
  if (domain_name == std::string("density") ||
      domain_name == std::string("molar") ||
      domain_name == std::string("mass")) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // -- temperature
  temp_key_ = plist_.get<std::string>("temperature key",
          domain_name+std::string("temperature"));
  dependencies_.insert(temp_key_);

  // -- pressure
  pres_key_ = plist_.get<std::string>("pressure key", "atmospheric_pressure");

  // Construct my EOS model
  AMANZI_ASSERT(plist_.isSublist("EOS parameters"));
  EOSFactory eos_fac;
  eos_ = eos_fac.CreateEOS(plist_.sublist("EOS parameters"));
};


IsobaricEOSEvaluator::IsobaricEOSEvaluator(const IsobaricEOSEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    eos_(other.eos_),
    mode_(other.mode_),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_) {};


Teuchos::RCP<FieldEvaluator> IsobaricEOSEvaluator::Clone() const {
  return Teuchos::rcp(new IsobaricEOSEvaluator(*this));
}


void IsobaricEOSEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const double> pres = S->GetScalarData(pres_key_);

  int index = 0; // index to the results list
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    // evaluate MolarDensity()
    Teuchos::Ptr<CompositeVector> result = results[index];
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
      Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

      int count = result->size(*comp);
      for (int id=0; id!=count; ++id) {
        result_v[0][id] = eos_->MolarDensity(temp_v[0][id], *pres);
      }
    }
    index++;
  }

  if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass()) {
    // calculate MassDensity from MolarDensity and molar mass.
    double M = eos_->MolarMass();
    results[1]->Update(M, *(results[0]), 0.0);
  } else if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    // evaluate MassDensity()
    Teuchos::Ptr<CompositeVector> result = results[index];
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
      Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

      int count = result->size(*comp);
      for (int id=0; id!=count; ++id) {
        result_v[0][id] = eos_->MassDensity(temp_v[0][id], *pres);
      }
    }
  }
}


void IsobaricEOSEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const double> pres = S->GetScalarData(pres_key_);

  if (wrt_key == temp_key_) {

    int index = 0; // index to the results list
    if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
      // evaluate DMolarDensityDT()
      Teuchos::Ptr<CompositeVector> result = results[index];
      for (CompositeVector::name_iterator comp=result->begin();
           comp!=result->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

        int count = result->size(*comp);
        for (int id=0; id!=count; ++id) {
          result_v[0][id] = eos_->DMolarDensityDT(temp_v[0][id], *pres);
        }
      }
      index++;
    }

    if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass()) {
      // calculate DMassDensityDT from DMolarDensityDT and molar mass.
      double M = eos_->MolarMass();
      results[1]->Update(M, *results[0], 0.0);
    } else if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
      // evaluate DMassDensityDT()
      Teuchos::Ptr<CompositeVector> result = results[index];
      for (CompositeVector::name_iterator comp=result->begin();
           comp!=result->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

        int count = result->size(*comp);
        for (int id=0; id!=count; ++id) {
          result_v[0][id] = eos_->DMassDensityDT(temp_v[0][id], *pres);
        }
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace AmanziEOS
}  // namespace Amanzi
