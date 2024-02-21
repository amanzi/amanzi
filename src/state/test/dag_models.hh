/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//!

#ifndef DAG_MODELS_HH_
#define DAG_MODELS_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "dbc.hh"
#include "AmanziTypes.hh"
#include "StateDefs.hh"

/*
  We will build the following dependencies tree:
    A -> {B, C, E, H}
    C -> {D, G}
    E -> {D, F}
    H -> F
    D -> G
    F -> G

  Primary fields are B=2 and G=3. The equations are
    A = 2*B + C*E*H = 6484
    C = 2*D + G     = 15
    E = D*F         = 36
    H = 2*F         = 12
    D = 2*G         = 6
    F = 2*G         = 6

  Derivatives are
    dA/dB = 2
    dA/dG = 8640

  WARNING: derivative of secondary field wrt to secondary field is
  not well defined. The code may throw an exception since
  intermediate derivatives are not saved.
*/

using namespace Amanzi;

namespace Amanzi {
class State;
}

template <class cView_type, class View_type>
class AModel {
 public:
  static const int n_dependencies = 4;
  static const std::string eval_type;

  AModel(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : alpha_(plist->sublist("model parameters").get<double>("alpha"))
  {}

  void setViews(const std::vector<cView_type>& dependency_views,
                const std::vector<View_type>& result_views,
                const State& s)
  {
    AMANZI_ASSERT(dependency_views.size() == 4);
    AMANZI_ASSERT(result_views.size() == 1);

    A_ = result_views[0];
    B_ = dependency_views[0];
    C_ = dependency_views[1];
    E_ = dependency_views[2];
    H_ = dependency_views[3];
  }

  KeyVector getMyKeys()
  {
    // NOTE, a real Model would parse the parameter list to get these
    return {
      "A",
    };
  };

  KeyVector getDependencies()
  {
    // NOTE, a real Model would parse the parameter list to get these
    return { "B", "C", "E", "H" };
  };

  // the model
  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    A_(i, 0) = alpha_ * B_(i, 0) + C_(i, 0) * E_(i, 0) * H_(i, 0);
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to getDependencies().  Deriv<I> must correspond to the derivative
  // with respect to getDependencies()(i,0);

  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { A_(i, 0) = alpha_; }
  // d/dC
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const
  {
    A_(i, 0) = E_(i, 0) * H_(i, 0);
  }
  // d/dE
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const int i) const
  {
    A_(i, 0) = C_(i, 0) * H_(i, 0);
  }
  // d/dH
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<3>, const int i) const
  {
    A_(i, 0) = C_(i, 0) * E_(i, 0);
  }

 private:
  View_type A_;
  cView_type B_;
  cView_type C_;
  cView_type E_;
  cView_type H_;
  double alpha_;
};


template <class cView_type, class View_type>
class CModel {
 public:
  static const int n_dependencies = 2;
  static const std::string eval_type;

  CModel(const Teuchos::RCP<Teuchos::ParameterList>& plist) {}

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(deps.size() == 2);
    AMANZI_ASSERT(res.size() == 1);

    C_ = res[0];
    D_ = deps[0];
    G_ = deps[1];
  }

  KeyVector getMyKeys()
  {
    return {
      "C",
    };
  };

  KeyVector getDependencies() { return { "D", "G" }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const { C_(i, 0) = 2 * D_(i, 0) + G_(i, 0); }

  // d/dD
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { C_(i, 0) = 2.0; }
  // d/dG
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { C_(i, 0) = 1.0; }


 private:
  View_type C_;
  cView_type D_;
  cView_type G_;
};


template <class cView_type, class View_type, bool UseAccessor = false>
class DModel_ {
 public:
  static const int n_dependencies = 1;
  static const std::string eval_type;
  using Accessor_View_type = Kokkos::View<const int*, typename View_type::device_type>;

  DModel_(const Teuchos::RCP<Teuchos::ParameterList>& plist) {}

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(deps.size() == 1);
    AMANZI_ASSERT(res.size() == 1);

    D_ = res[0];
    G_ = deps[0];
  }

  void setAccessor(const Accessor_View_type& ids) { ids_ = ids; }

  KeyVector getMyKeys()
  {
    return {
      "D",
    };
  };

  KeyVector getDependencies()
  {
    return {
      "G",
    };
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const
  {
    int j = UseAccessor ? ids_(i) : i;
    D_(i, 0) = 2 * G_(j, 0);
  }

  // derivatives
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { D_(i, 0) = 2.0; }

 private:
  View_type D_;
  cView_type G_;

  Kokkos::View<const int*, typename View_type::device_type> ids_;
};


template <class cView_type, class View_type>
using DModel = DModel_<cView_type, View_type, false>;
template <class cView_type, class View_type>
using DModelAccessor = DModel_<cView_type, View_type, true>;


template <class cView_type, class View_type>
class EModel {
 public:
  static const int n_dependencies = 2;
  static const std::string eval_type;

  EModel(const Teuchos::RCP<Teuchos::ParameterList>& plist) {}

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(deps.size() == 2);
    AMANZI_ASSERT(res.size() == 1);

    E_ = res[0];
    D_ = deps[0];
    F_ = deps[1];
  }

  KeyVector getMyKeys()
  {
    return {
      "E",
    };
  };

  KeyVector getDependencies() { return { "D", "F" }; }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const { E_(i, 0) = D_(i, 0) * F_(i, 0); }

  // d/dD
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { E_(i, 0) = F_(i, 0); }
  // d/dF
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const { E_(i, 0) = D_(i, 0); }


 private:
  View_type E_;
  cView_type D_;
  cView_type F_;
};


template <class cView_type, class View_type>
class FModel {
 public:
  static const int n_dependencies = 1;
  static const std::string eval_type;

  FModel(const Teuchos::RCP<Teuchos::ParameterList>& plist) {}

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(deps.size() == 1);
    AMANZI_ASSERT(res.size() == 1);

    F_ = res[0];
    G_ = deps[0];
  }

  KeyVector getMyKeys()
  {
    return {
      "F",
    };
  };

  KeyVector getDependencies()
  {
    return {
      "G",
    };
  }

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const { F_(i, 0) = 2 * G_(i, 0); }

  // d/dG
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { F_(i, 0) = 2.0; }


 private:
  View_type F_;
  cView_type G_;
};


template <class cView_type, class View_type>
class HModel {
 public:
  static const int n_dependencies = 1;
  static const std::string eval_type;

  HModel(const Teuchos::RCP<Teuchos::ParameterList>& plist) {}

  void
  setViews(const std::vector<cView_type>& deps, const std::vector<View_type>& res, const State& s)
  {
    AMANZI_ASSERT(deps.size() == 1);
    AMANZI_ASSERT(res.size() == 1);

    H_ = res[0];
    F_ = deps[0];
  }

  KeyVector getMyKeys()
  {
    return {
      "H",
    };
  };

  KeyVector getDependencies()
  {
    return {
      "F",
    };
  }


  KOKKOS_INLINE_FUNCTION void operator()(const int i) const { H_(i, 0) = 2 * F_(i, 0); }

  // d/dF
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const { H_(i, 0) = 2.0; }


 private:
  View_type H_;
  cView_type F_;
};


template <class cView_type, class View_type>
const std::string AModel<cView_type, View_type>::eval_type("A");

template <class cView_type, class View_type>
const std::string CModel<cView_type, View_type>::eval_type("C");

template <class cView_type, class View_type, bool UseAccessor>
const std::string DModel_<cView_type, View_type, UseAccessor>::eval_type("D");

template <class cView_type, class View_type>
const std::string EModel<cView_type, View_type>::eval_type("E");

template <class cView_type, class View_type>
const std::string FModel<cView_type, View_type>::eval_type("F");

template <class cView_type, class View_type>
const std::string HModel<cView_type, View_type>::eval_type("H");


#endif
