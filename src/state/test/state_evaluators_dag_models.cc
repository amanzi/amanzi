/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/


//! <MISSING_ONELINE_DOCSTRING>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "State.hh"
#include "evaluator/EvaluatorSecondaryMonotype.hh"
#include "evaluator/EvaluatorPrimary.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziMesh;

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

// MOVE THESE TO UTILS
//
// Tag type for derivatives
template<int> struct Deriv {};

//
// END MOVE THESE TO UTILS



template<class cView_type, class View_type>
class AModel {
 public:
  static const int n_dependencies = 4;
  
  AModel(Teuchos::ParameterList& plist) :
      alpha_(plist.get<double>("alpha", 2.0)) {}

  void SetViews(const std::vector<cView_type>& deps,
                const std::vector<View_type>& res)
  {
    AMANZI_ASSERT(deps.size() == 4);
    AMANZI_ASSERT(res.size() == 1);
    
    A_ = res[0];
    B_ = deps[0];
    C_ = deps[1];
    E_ = deps[2];
    H_ = deps[3];
  }

  KeyVector my_keys() {
    return { "A", };
  };

  KeyPairVector dependencies() {
    return { {"B", ""}, {"C", ""}, {"E", ""}, {"H", ""} };
  };
  
  KOKKOS_INLINE_FUNCTION void operator() (const int i) const {
    A_[i] = alpha_ * B_[i] + C_[i] * E_[i] * H_[i];
  }

  // generic match, throws
  template<class Deriv_type>
  void operator()(Deriv_type, const typename View_type::size_type i) const {
    AMANZI_ASSERT(false);
  }
  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const typename View_type::size_type i) const {
    A_[i] = alpha_;
  }
  // d/dC
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const typename View_type::size_type i) const {
    A_[i] = E_[i] * H_[i];
  }
  // d/dE
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const typename View_type::size_type i) const {
    A_[i] = C_[i] * H_[i];
  }
  // d/dH
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<3>, const typename View_type::size_type i) const {
    A_[i] = C_[i] * E_[i];
  }

  
private:
  View_type A_;
  cView_type B_;
  cView_type C_;
  cView_type E_;
  cView_type H_;
  double alpha_;
};


template<class cView_type, class View_type>
class CModel {
 public:
  static const int n_dependencies = 1;

  CModel(Teuchos::ParameterList& plist) {}

  void SetViews(const std::vector<cView_type >& deps,
                const std::vector<View_type >& res)
  {
    AMANZI_ASSERT(deps.size() == 2);
    AMANZI_ASSERT(res.size() == 1);
    
    C_ = res[0];
    D_ = deps[0];
    G_ = deps[1];
  }

  KeyVector my_keys() {
    return { "C", };
  };

  KeyPairVector dependencies() {
    return { {"D", ""}, {"G", ""} };
  };
  
  KOKKOS_INLINE_FUNCTION void operator() (const int i) const {
    C_[i] = 2 * D_[i] + G_[i];
  }

  // generic match, throws
  template<class Deriv_type>
  void operator()(Deriv_type, const typename View_type::size_type i)  const {
    AMANZI_ASSERT(false);
  }

  // d/dD
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const typename View_type::size_type i)  const {
    C_[i] = 2.0;
  }
  // d/dG
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const typename View_type::size_type i)  const {
    C_[i] = 1.0;
  }

  
private:
  View_type C_;
  cView_type D_;
  cView_type G_;
};


template<class cView_type, class View_type>
class DModel {
 public:
  static const int n_dependencies = 1;

  DModel(Teuchos::ParameterList& plist) {}

  void SetViews(const std::vector<cView_type >& deps,
                const std::vector<View_type >& res)
  {
    AMANZI_ASSERT(deps.size() == 1);
    AMANZI_ASSERT(res.size() == 1);
    
    D_ = res[0];
    G_ = deps[0];
  }

  KeyVector my_keys() {
    return { "D", };
  };

  KeyPairVector dependencies() {
    return { {"G", ""} };
  };

  KOKKOS_INLINE_FUNCTION void operator()(const typename View_type::size_type i) const {
    D_[i] = 2 * G_[i];
  }

  // generic match, throws
  template<class Deriv_type>
  void operator()(Deriv_type, const typename View_type::size_type i)  const {
    AMANZI_ASSERT(false);
  }

  // derivatives
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const typename View_type::size_type i)  const {
    D_[i] = 2.0;
  }
  
private:
  View_type D_;
  cView_type G_;
};


template<class cView_type, class View_type>
class EModel {
 public:
  static const int n_dependencies = 2;

  EModel(Teuchos::ParameterList& plist) {}

  void SetViews(const std::vector<cView_type >& deps,
                const std::vector<View_type >& res)
  {
    AMANZI_ASSERT(deps.size() == 2);
    AMANZI_ASSERT(res.size() == 1);
    
    E_ = res[0];
    D_ = deps[0];
    F_ = deps[1];
  }

  KeyVector my_keys() {
    return { "E", };
  };

  KeyPairVector dependencies() {
    return { {"D", ""}, {"F", ""} };
  };
  
  KOKKOS_INLINE_FUNCTION void operator() (const int i) const {
    E_[i] = D_[i] * F_[i];
  }

  // generic match, throws
  template<class Deriv_type>
  void operator()(Deriv_type, const typename View_type::size_type i) const {
    AMANZI_ASSERT(false);
  }
  // d/dD
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const typename View_type::size_type i) const {
    E_[i] = F_[i];
  }
  // d/dF
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const typename View_type::size_type i) const {
    E_[i] = D_[i];
  }

  
private:
  View_type E_;
  cView_type D_;
  cView_type F_;
};


template<class cView_type, class View_type>
class FModel {
 public:
  static const int n_dependencies = 1;
  
  FModel(Teuchos::ParameterList& plist) {}

  void SetViews(const std::vector<cView_type >& deps,
                const std::vector<View_type >& res)
  {
    AMANZI_ASSERT(deps.size() == 1);
    AMANZI_ASSERT(res.size() == 1);
    
    F_ = res[0];
    G_ = deps[0];
  }

  KeyVector my_keys() {
    return { "F", };
  };

  KeyPairVector dependencies() {
    return { {"G", ""} };
  };
  
  KOKKOS_INLINE_FUNCTION void operator() (const int i) const {
    F_[i] = 2 * G_[i];
  }

  // generic match, throws
  template<class Deriv_type>
  void operator()(Deriv_type, const typename View_type::size_type i) const {
    AMANZI_ASSERT(false);
  }
  // d/dG
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const typename View_type::size_type i) const {
    F_[i] = 2.0;
  }
  

private:
  View_type F_;
  cView_type G_;
};


template<class cView_type, class View_type>
class HModel {
 public:
  static const int n_dependencies = 1;

  HModel(Teuchos::ParameterList& plist) {}

  void SetViews(const std::vector<cView_type >& deps,
                const std::vector<View_type >& res)
  {
    AMANZI_ASSERT(deps.size() == 1);
    AMANZI_ASSERT(res.size() == 1);
    
    H_ = res[0];
    F_ = deps[0];
  }

  KeyVector my_keys() {
    return { "H", };
  };
  
  KeyPairVector dependencies() {
    return { {"F", ""} };
  };

  
  KOKKOS_INLINE_FUNCTION void operator() (const int i) const {
    H_[i] = 2 * F_[i];
  }

  // generic match, throws
  template<class Deriv_type>
  void operator()(Deriv_type, const typename View_type::size_type i) const {
    AMANZI_ASSERT(false);
  }
  // d/dF
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const typename View_type::size_type i) const {
    H_[i] = 2.0;
  }

  
private:
  View_type H_;
  cView_type F_;
};



template<template<class,class> class Model, class Device_type=AmanziDefaultDevice>
class EvaluatorModel_CompositeVector : public EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace> {
 public:
  using Model_type = Model<cVectorView_type<Device_type>, VectorView_type<Device_type>>;
  using View_type = VectorView_type<Device_type>;
  using cView_type = cVectorView_type<Device_type>;
  
  EvaluatorModel_CompositeVector(Teuchos::ParameterList& plist)
      : EvaluatorSecondaryMonotype<CompositeVector,CompositeVectorSpace>(plist),
        model_(Teuchos::rcp(new Model_type(plist))),
        name_(plist.name())
  {
    dependencies_ = model_->dependencies();
  }

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EvaluatorModel_CompositeVector<Model,Device_type>(*this));
  }

  virtual void
  Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override
  {
    AMANZI_ASSERT(results.size() > 0);
    for (const auto& comp : *results[0]) {
      std::vector<cView_type> dependency_views;
      for (const auto& dep : dependencies_) {
        const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
        if (vec.getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorModel: expects only dependencies with one DoF.");
          throw(msg);
        }        
        dependency_views.emplace_back(vec.ViewComponent(comp, 0, false));
      }

      // FIXME -- no defaultdevice!
      std::vector<View_type> result_views;
      for (auto result : results) {
        if (result->getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorModel: expects only results with one DoF.");
          throw(msg);
        }        
        result_views.emplace_back(result->ViewComponent(comp, 0, false));
      }

      // set up the model and range and then dispatch
      model_->SetViews(dependency_views, result_views);
      Kokkos::RangePolicy<typename Device_type::execution_space> range(0, result_views[0].extent(0));
      Kokkos::parallel_for(name_, range, *model_);
    }
  }

  virtual void
  EvaluatePartialDerivative_(const State& S, const Key& wrt_key,
                             const Key& wrt_tag,
                             const std::vector<CompositeVector*>& results) override {
    AMANZI_ASSERT(results.size() > 0);
    for (const auto& comp : *results[0]) {
      std::vector<cView_type> dependency_views;
      for (const auto& dep : dependencies_) {
        const auto& vec = S.Get<CompositeVector>(dep.first, dep.second);
        if (vec.getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorModel: expects only dependencies with one DoF.");
          throw(msg);
        }        
        dependency_views.emplace_back(vec.ViewComponent(comp, 0, false));
      }

      // FIXME -- no defaultdevice!
      std::vector<View_type> result_views;
      for (auto result : results) {
        if (result->getNumVectors(comp) != 1) {
          Errors::Message msg("EvaluatorModel: expects only results with one DoF.");
          throw(msg);
        }        
        result_views.emplace_back(result->ViewComponent(comp, 0, false));
      }

      // set up the model and range and then dispatch
      model_->SetViews(dependency_views, result_views);

      // which derivative?  This is a bit awkward to ensure compile-time
      // constants.  Can this be fixed?
      auto wrt = std::make_pair(wrt_key, wrt_tag);

      
      
      if (dependencies_.size() > 0 && wrt == dependencies_[0]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<0>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 1 && wrt == dependencies_[1]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<1>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 2 && wrt == dependencies_[2]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<2>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 3 && wrt == dependencies_[3]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<3>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 4 && wrt == dependencies_[4]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<4>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 5 && wrt == dependencies_[5]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<5>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 6 && wrt == dependencies_[6]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<6>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 7 && wrt == dependencies_[7]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<7>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 8 && wrt == dependencies_[8]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<8>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 9 && wrt == dependencies_[9]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<9>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      } else if (dependencies_.size() > 10 && wrt == dependencies_[10]) {
        // NOTE: the 0 must be compile-time.  HELP.
        Kokkos::RangePolicy<Deriv<10>, typename Device_type::execution_space> range(0, result_views[0].extent(0));
        Kokkos::parallel_for(name_, range, *model_);

      }      
    }
  }

  Teuchos::RCP<Model_type> model_;
  std::string name_;
};


class make_state {
 public:
  make_state()
  {
    // create a mesh
    auto comm = Amanzi::getDefaultComm();
    MeshFactory meshfac(comm);
    auto mesh = meshfac.create(0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 2, 2, 2);

    // create a state
    // State S;
    S.RegisterDomainMesh(mesh);

    // Secondary fields
    fa_eval = requireSecondary<AModel>("A", { "B", "G" });
    fc_eval = requireSecondary<CModel>("C", {});
    fd_eval = requireSecondary<DModel>("D", {});
    fe_eval = requireSecondary<EModel>("E", { "G", });
    ff_eval = requireSecondary<FModel>("F", {});
    fh_eval = requireSecondary<HModel>("H", {});

    // Primary fields
    fb_eval = requirePrimary("B");
    fg_eval = requirePrimary("G");

    // Setup fields initialize
    S.Setup();
    S.GetW<CompositeVector>("B", "B").putScalar(2.0);
    S.GetRecordW("B", "B").set_initialized();
    S.GetW<CompositeVector>("G", "G").putScalar(3.0);
    S.GetRecordW("G", "G").set_initialized();
    S.Initialize();
  }

  template<template<class,class> class Model>
  Teuchos::RCP<Evaluator>
  requireSecondary(const std::string& name, const std::vector<std::string>& derivs) {
    Teuchos::ParameterList es_list(name);
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName(name);
    es_list.set("tag", "");
    S.Require<CompositeVector,CompositeVectorSpace>(name, "", name)
        .SetMesh(S.GetMesh())
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    for (const auto& deriv : derivs)
      S.RequireDerivative<CompositeVector,CompositeVectorSpace>(name, "", deriv, "");
    
    //  S.RequireDerivative<double>("fa", "", "fb", "");
    //  S.RequireDerivative<double>("fa", "", "fg", "");
    auto f_eval = Teuchos::rcp(new EvaluatorModel_CompositeVector<Model>(es_list));
    S.SetEvaluator(name, f_eval);
    return f_eval;
  }


  Teuchos::RCP<EvaluatorPrimary_>
  requirePrimary(const std::string& name) {
    Teuchos::ParameterList es_list(name);
    es_list.sublist("verbose object")
        .set<std::string>("verbosity level", "extreme");
    es_list.setName(name);
    es_list.set("tag", "");
    S.Require<CompositeVector,CompositeVectorSpace>(name, "", name)
        .SetMesh(S.GetMesh())
        ->SetComponent("cell", AmanziMesh::CELL, 1);
    
    //  S.RequireDerivative<double>("fa", "", "fb", "");
    //  S.RequireDerivative<double>("fa", "", "fg", "");
    auto f_eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector,CompositeVectorSpace>(es_list));
    S.SetEvaluator(name, f_eval);
    return f_eval;
  }

  void check_close(double val, const std::string& name) {
    auto cvv = S.Get<CompositeVector>(name, "").ViewComponent<AmanziDefaultHost>("cell", 0, false);
    CHECK_CLOSE(val, cvv(0), 1.e-10);
  }

  void check_close_deriv(double val, const std::string& name, const std::string& wrt) {
    auto cvv = S.GetDerivative<CompositeVector>(name, "", wrt, "").ViewComponent<AmanziDefaultHost>("cell", 0, false);
    CHECK_CLOSE(val, cvv(0), 1.e-10);
  }
  
  
 public:
  State S;
  Teuchos::RCP<Evaluator> fa_eval;
  Teuchos::RCP<Evaluator> fc_eval;
  Teuchos::RCP<Evaluator> fd_eval;
  Teuchos::RCP<Evaluator> fe_eval;
  Teuchos::RCP<Evaluator> ff_eval;
  Teuchos::RCP<Evaluator> fh_eval;
  Teuchos::RCP<EvaluatorPrimary_> fb_eval, fg_eval;
};

SUITE(DAG)
{
  TEST_FIXTURE(make_state, DAG_TWO_FIELDS)
  {
    // check initialized properly
    check_close(2.0,"B");
    check_close(3.0,"G");

    // calculate field A
    std::cout << "Calculate field A:" << std::endl;
    bool changed = fa_eval->Update(S, "main");
    check_close(6484.0,"A");
    CHECK(changed);

    // check intermediate steps got updated too
    check_close(6.0, "D");

    // calculate dA/dB
    std::cout << "Calculate derivative of field A wrt field B:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "A", "B", "");
    check_close_deriv(2.0, "A", "B");
    CHECK(changed);

    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "A", "G", "");
    check_close_deriv(8640.0, "A", "G");
    CHECK(changed);

    // calculate dE/dG:
    std::cout << "Calculate derivative of field E wrt field G:" << std::endl;
    changed = fe_eval->UpdateDerivative(S, "E", "G", "");
    check_close_deriv(24.0, "E", "G");
    CHECK(changed);

    // Now we repeat some calculations. Since no primary fields changed,
    // the result should be the same
    // calculate dA/dG
    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    changed = fa_eval->UpdateDerivative(S, "A", "G", "");
    check_close_deriv(8640.0, "A", "G");
    CHECK(!changed);

    std::cout << "Calculate derivative of field A wrt field G:" << std::endl;
    fb_eval->SetChanged();
    changed = fa_eval->UpdateDerivative(S, "A", "G", "");
    check_close_deriv(8640.0, "A", "G");
    CHECK(changed);
  }
}
