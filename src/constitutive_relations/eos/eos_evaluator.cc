/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_factory.hh"
#include "eos_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,EOSEvaluator> EOSEvaluator::factory_("eos");


EOSEvaluator::EOSEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist) {

  // Process the list for my provided field.
  std::string mode = plist_.get<std::string>("EOS basis", "molar");
  if (mode == "molar") {
    mode_ = EOS_MODE_MOLAR;
  } else if (mode == "mass") {
    mode_ = EOS_MODE_MASS;
  } else if (mode == "both") {
    mode_ = EOS_MODE_BOTH;
  } else {
    ASSERT(0);
  }

  // my keys
  Key key;
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    key = plist_.get<std::string>("molar density key");
    my_keys_.push_back(key);
  }

  if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    key = plist_.get<std::string>("mass density key");
    my_keys_.push_back(key);
  }
  setLinePrefix(my_keys_[0]+std::string(" evaluator"));

  // Set up my dependencies.
  std::size_t end = key.find_first_of("_");
  std::string domain_name = key.substr(0,end);
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
  pres_key_ = plist_.get<std::string>("pressure key",
          domain_name+std::string("effective_pressure"));
  dependencies_.insert(pres_key_);

  // -- logging
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME)) {
  Teuchos::OSTab tab = getOSTab();
    for (KeySet::const_iterator dep=dependencies_.begin();
         dep!=dependencies_.end(); ++dep) {
      *out_ << " dep: " << *dep << std::endl;
    }
  }

  // Construct my EOS model
  ASSERT(plist_.isSublist("EOS parameters"));
  EOSFactory eos_fac;
  eos_ = eos_fac.createEOS(plist_.sublist("EOS parameters"));
};


EOSEvaluator::EOSEvaluator(const EOSEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    eos_(other.eos_),
    mode_(other.mode_),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_) {}


Teuchos::RCP<FieldEvaluator> EOSEvaluator::Clone() const {
  return Teuchos::rcp(new EOSEvaluator(*this));
}


void EOSEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  Teuchos::Ptr<CompositeVector> molar_dens, mass_dens;
  if (mode_ == EOS_MODE_MOLAR) {
    molar_dens = results[0];
  } else if (mode_ == EOS_MODE_MASS) {
    mass_dens = results[0];
  } else {
    molar_dens = results[0];
    mass_dens = results[1];
  }


  if (molar_dens != Teuchos::null) {
    // evaluate MolarDensity()
    for (CompositeVector::name_iterator comp=molar_dens->begin();
         comp!=molar_dens->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
      const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
      Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

      int count = dens_v.MyLength();
      for (int id=0; id!=count; ++id) {
        dens_v[0][id] = eos_->MolarDensity(temp_v[0][id], pres_v[0][id]);
      }
    }
  }

  if (mass_dens != Teuchos::null) {
    for (CompositeVector::name_iterator comp=mass_dens->begin();
         comp!=mass_dens->end(); ++comp) {
      if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
          molar_dens->has_component(*comp)) {
        // calculate MassDensity from MolarDensity and molar mass.
        double M = eos_->MolarMass();

        mass_dens->ViewComponent(*comp,false)->Update(M,
                *molar_dens->ViewComponent(*comp,false), 0.);
      } else {
        // evaluate MassDensity() directly
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
        Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int id=0; id!=count; ++id) {
          dens_v[0][id] = eos_->MassDensity(temp_v[0][id], pres_v[0][id]);
        }
      }
    }
  }
}


void EOSEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) {

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  Teuchos::Ptr<CompositeVector> molar_dens, mass_dens;
  if (mode_ == EOS_MODE_MOLAR) {
    molar_dens = results[0];
  } else if (mode_ == EOS_MODE_MASS) {
    mass_dens = results[0];
  } else {
    molar_dens = results[0];
    mass_dens = results[1];
  }

  if (wrt_key == pres_key_) {
    if (molar_dens != Teuchos::null) {
      // evaluate MolarDensity()
      for (CompositeVector::name_iterator comp=molar_dens->begin();
           comp!=molar_dens->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
        Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int id=0; id!=count; ++id) {
          dens_v[0][id] = eos_->DMolarDensityDp(temp_v[0][id], pres_v[0][id]);
        }
      }
    }

    if (mass_dens != Teuchos::null) {
      for (CompositeVector::name_iterator comp=mass_dens->begin();
           comp!=mass_dens->end(); ++comp) {
        if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
            molar_dens->has_component(*comp)) {
          // calculate MassDensity from MolarDensity and molar mass.
          double M = eos_->MolarMass();

          mass_dens->ViewComponent(*comp,false)->Update(M,
                  *molar_dens->ViewComponent(*comp,false), 0.);
        } else {
          // evaluate MassDensity() directly
          const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
          const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
          Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

          int count = dens_v.MyLength();
          for (int id=0; id!=count; ++id) {
            dens_v[0][id] = eos_->DMassDensityDp(temp_v[0][id], pres_v[0][id]);
          }
        }
      }
    }

  } else if (wrt_key == temp_key_) {

    if (molar_dens != Teuchos::null) {
      // evaluate MolarDensity()
      for (CompositeVector::name_iterator comp=molar_dens->begin();
           comp!=molar_dens->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
        const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
        Epetra_MultiVector& dens_v = *(molar_dens->ViewComponent(*comp,false));

        int count = dens_v.MyLength();
        for (int id=0; id!=count; ++id) {
          dens_v[0][id] = eos_->DMolarDensityDT(temp_v[0][id], pres_v[0][id]);
        }
      }
    }

    if (mass_dens != Teuchos::null) {
      for (CompositeVector::name_iterator comp=mass_dens->begin();
           comp!=mass_dens->end(); ++comp) {
        if (mode_ == EOS_MODE_BOTH && eos_->IsConstantMolarMass() &&
            molar_dens->has_component(*comp)) {
          // calculate MassDensity from MolarDensity and molar mass.
          double M = eos_->MolarMass();

          mass_dens->ViewComponent(*comp,false)->Update(M,
                  *molar_dens->ViewComponent(*comp,false), 0.);
        } else {
          // evaluate MassDensity() directly
          const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
          const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp,false));
          Epetra_MultiVector& dens_v = *(mass_dens->ViewComponent(*comp,false));

          int count = dens_v.MyLength();
          for (int id=0; id!=count; ++id) {
            dens_v[0][id] = eos_->DMassDensityDT(temp_v[0][id], pres_v[0][id]);
          }
        }
      }
    }

  } else {
    ASSERT(0);
  }
}

} // namespace
} // namespace
