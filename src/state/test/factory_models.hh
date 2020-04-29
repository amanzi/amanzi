/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Two quick physics models for use in testing.


#ifndef FACTORY_MODELS_HH_
#define FACTORY_MODELS_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziTypes.hh"
#include "evaluator/EvaluatorModel_CompositeVector.hh"
#include "evaluator/EvaluatorModelByMaterial.hh"

using namespace Amanzi;

//
// Two quick models for use in test.  Note these are actual physics models, and
// could really be in physics code.
template <class cView_type, class View_type>
class ModelCompressiblePorosity {
 public:
  static const int n_dependencies = 2;

  ModelCompressiblePorosity(Teuchos::ParameterList& plist)
      : compressibility_(plist.sublist("model parameters").get<double>("pore compressibility [Pa^-1]")),
        cutoff_(plist.sublist("model parameters").get<double>("pore compressibility inflection point [Pa]", 1000.)),
        p_atm_(plist.sublist("model parameters").get<double>("atmospheric pressure [Pa]", 101325.))
  {
    std::string result_key_ = Keys::cleanPListName(plist.name());
    std::string domain_name = Keys::getDomain(result_key_);
    base_poro_key_ = Keys::readKey(plist, domain_name, "base porosity", "base_porosity");
    pressure_key_ = Keys::readKey(plist, domain_name, "pressure", "pressure");
  }

  void SetViews(const std::vector<cView_type>& dependency_views,
                const std::vector<View_type>& result_views,
                const State& s)
  {
    AMANZI_ASSERT(dependency_views.size() == 2);
    AMANZI_ASSERT(result_views.size() == 1);

    result_ = result_views[0];
    base_poro_ = dependency_views[0];
    pressure_ = dependency_views[1];
  }

  KeyVector my_keys()
  {
    return { result_key_, };
  };

  KeyVector dependencies()
  {
    return { base_poro_key_, pressure_key_ };
  };

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    double poro = base_poro_(i);
    double p_over = pressure_(i) - p_atm_;
    if (p_over > cutoff_) {
      poro += compressibility_ * ( cutoff_ / 2. + (p_over - cutoff_));
    } else if (p_over > 0.) {
      poro += compressibility_ * (std::pow(p_over,2.) / 2. / cutoff_);
    }
    result_(i) = poro;
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  // d/d base_poro
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    result_(i) = pressure_(i) > p_atm_ ? 1.0 : 0.;
  }

  // d/dpressure
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const
  {
    double p_over = pressure_(i) - p_atm_;
    double deriv = 0.;
    if (p_over > cutoff_) {
      deriv = compressibility_;
    } else if (p_over > 0.) {
      deriv = compressibility_ * p_over / cutoff_;
    }
    result_(i) = deriv;
  }

 private:
  View_type result_;
  cView_type base_poro_;
  cView_type pressure_;

  double compressibility_;
  double cutoff_;
  double p_atm_;

  Key result_key_;
  Key base_poro_key_;
  Key pressure_key_;

};



//
// Two quick models for use in test.  Note these are actual physics models, and
// could really be in physics code.
template <class cView_type, class View_type>
class ModelWaterRetentionVanGenuchten {
 public:
  static const int n_dependencies = 1;

  ModelWaterRetentionVanGenuchten(Teuchos::ParameterList& plist)
      : alpha_(plist.sublist("model parameters").get<double>("alpha [Pa^-1]")),
        n_(plist.sublist("model parameters").get<double>("n [-]", -1.0)),
        m_(plist.sublist("model parameters").get<double>("m [-]", -1.0)),
        sr_(plist.sublist("model parameters").get<double>("residual saturation [-]"))
  {
    // reconcile n and m
    if (n_ < 0.) {
      if (m_ < 0.) {
        Errors::Message msg("Van Genuchten water retention model got neither n nor m parameters.");
        throw(msg);

      } else {
        n_ = 1.0 / (1.0 - m_);
      }
    } else {
      m_ = 1.0 - 1.0/n_;
    }

    // keys
    std::string result_key_ = Keys::cleanPListName(plist.name());
    std::string domain_name = Keys::getDomain(result_key_);
    cap_pressure_key_ = Keys::readKey(plist, domain_name, "capillary pressure", "capillary_pressure_gas_liquid");
    
  }

  void SetViews(const std::vector<cView_type>& dependency_views,
                const std::vector<View_type>& result_views,
                const State& s)
  {
    AMANZI_ASSERT(dependency_views.size() == 1);
    AMANZI_ASSERT(result_views.size() == 1);

    result_ = result_views[0];
    cap_pressure_ = dependency_views[0];
  }

  KeyVector my_keys()
  {
    return { result_key_, };
  };

  KeyVector dependencies()
  {
    return { cap_pressure_key_, };
  };

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    if (cap_pressure_(i) <= 0.) {
      result_(i) = 1.0;
    } else {
      result_(i) = std::pow(1.0 + std::pow(alpha_*cap_pressure_(i), n_), -m_) * (1.0 - sr_) + sr_;
    }
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];

  // d/d cap pressure
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const
  {
    if (cap_pressure_(i) <= 0.) {
      result_(i) = 0.;
    } else {
      result_(i) =  -m_*n_ * std::pow(1.0 + std::pow(alpha_*cap_pressure_(i), n_), -m_-1.0) * std::pow(alpha_*cap_pressure_(i), n_-1) * alpha_ * (1.0 - sr_);
    }      
  }

 private:
  View_type result_;
  cView_type cap_pressure_;

  double alpha_;
  double n_, m_;
  double sr_;

  Key result_key_;
  Key cap_pressure_key_;
};



#endif
