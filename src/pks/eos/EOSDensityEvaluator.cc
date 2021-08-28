/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOSFieldEvaluator is the interface between state/data and the 
  model, an EOS.
*/

#include "errors.hh"

#include "EOSDensityEvaluator.hh"
#include "EOSDensityFactory.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Constructor
****************************************************************** */
EOSDensityEvaluator::EOSDensityEvaluator(Teuchos::ParameterList& plist) :
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
  Key key;
  if (mode_ == EOS_MODE_MOLAR || mode_ == EOS_MODE_BOTH) {
    key = plist_.get<std::string>("molar density key");
    my_keys_.push_back(key);
  }

  if (mode_ == EOS_MODE_MASS || mode_ == EOS_MODE_BOTH) {
    key = plist_.get<std::string>("mass density key");
    my_keys_.push_back(key);
  }

  // set up my dependencies
  std::string domain = Keys::getDomain(key);
  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(temp_key_);

  pres_key_ = plist_.get<std::string>("pressure key", Keys::getKey(domain, "pressure"));
  dependencies_.insert(pres_key_);

  // logging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    for (auto dep = dependencies_.begin(); dep != dependencies_.end(); ++dep) {
      *vo_->os() << " dep: " << *dep << std::endl;
    }
  }

  // Construct my EOS model
  AMANZI_ASSERT(plist_.isSublist("EOS parameters"));
  EOSDensityFactory eos_fac;
  eos_ = eos_fac.CreateEOS(plist_.sublist("EOS parameters"));
};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
EOSDensityEvaluator::EOSDensityEvaluator(const EOSDensityEvaluator& other)
  : SecondaryVariablesFieldEvaluator(other),
    eos_(other.eos_),
    mode_(other.mode_),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_) {}


/* ******************************************************************
* Clone.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> EOSDensityEvaluator::Clone() const {
  return Teuchos::rcp(new EOSDensityEvaluator(*this));
}


/* ******************************************************************
* TBW.
****************************************************************** */
void EOSDensityEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  double p_atm = *S->GetScalarData("atmospheric_pressure");
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
    for (auto comp = molar_dens->begin(); comp != molar_dens->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
      Epetra_MultiVector& dens_v = *molar_dens->ViewComponent(*comp);

      int count = dens_v.MyLength();
      for (int i = 0; i != count; ++i) {
        double tmp = std::max<double>(pres_v[0][i], p_atm);
        dens_v[0][i] = eos_->MolarDensity(temp_v[0][i], tmp);
        if (dens_v[0][i] <= 0.0)
            Exceptions::amanzi_throw(Errors::CutTimeStep());
      }
    }
  }

  if (mass_dens != Teuchos::null) {
    for (auto comp = mass_dens->begin(); comp != mass_dens->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
      Epetra_MultiVector& dens_v = *mass_dens->ViewComponent(*comp);

      int count = dens_v.MyLength();
      for (int i = 0; i != count; ++i) {
        double tmp = std::max<double>(pres_v[0][i], p_atm);
        dens_v[0][i] = eos_->Density(temp_v[0][i], tmp);
        if (dens_v[0][i] < 0.0)
          Exceptions::amanzi_throw(Errors::CutTimeStep());
      }
    }
  }
}


/* ******************************************************************
* TBW.
****************************************************************** */
void EOSDensityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  double p_atm = *S->GetScalarData("atmospheric_pressure");
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
      for (auto comp = molar_dens->begin(); comp != molar_dens->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
        const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
        Epetra_MultiVector& dens_v = *molar_dens->ViewComponent(*comp);

        int count = dens_v.MyLength();
        for (int i = 0; i != count; ++i) {
          dens_v[0][i] = (pres_v[0][i] > p_atm) ? eos_->DMolarDensityDp(temp_v[0][i], pres_v[0][i]) : 0.0;
        }
      }
    }

    if (mass_dens != Teuchos::null) {
      for (auto comp = mass_dens->begin(); comp != mass_dens->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
        const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
        Epetra_MultiVector& dens_v = *mass_dens->ViewComponent(*comp);

        int count = dens_v.MyLength();
        for (int i = 0; i != count; ++i) {
          dens_v[0][i] = (pres_v[0][i] > p_atm) ? dens_v[0][i] = eos_->DDensityDp(temp_v[0][i], pres_v[0][i]) : 0.0;
        }
      }
    }

  } else if (wrt_key == temp_key_) {
    if (molar_dens != Teuchos::null) {
      for (auto comp = molar_dens->begin(); comp != molar_dens->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
        const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
        Epetra_MultiVector& dens_v = *molar_dens->ViewComponent(*comp);

        int count = dens_v.MyLength();
        for (int i = 0; i != count; ++i) {
          dens_v[0][i] = eos_->DMolarDensityDT(temp_v[0][i], pres_v[0][i]);
        }
      }
    }

    if (mass_dens != Teuchos::null) {
      for (auto comp = mass_dens->begin(); comp != mass_dens->end(); ++comp) {
        const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp);
        const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp);
        Epetra_MultiVector& dens_v = *mass_dens->ViewComponent(*comp);

        int count = dens_v.MyLength();
        for (int i = 0; i != count; ++i) {
          dens_v[0][i] = eos_->DDensityDT(temp_v[0][i], pres_v[0][i]);
        }
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace AmanziEOS
}  // namespace Amanzi

