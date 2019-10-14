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


//
// Tag type for derivatives
template<int> struct Deriv {};



template<class cView_type, class View_type>
class AModel {
 public:

  static const int n_dependencies = 4;
  static const std::string name;
  
  AModel(Teuchos::ParameterList& plist) :
      alpha_(plist.get<double>("alpha", 2.0)) {}

  void SetViews(const std::vector<cView_type>& dependency_views,
                const std::vector<View_type>& result_views)
  {
    AMANZI_ASSERT(dependency_views.size() == 4);
    AMANZI_ASSERT(result_views.size() == 1);
    
    A_ = result_views[0];
    B_ = dependency_views[0];
    C_ = dependency_views[1];
    E_ = dependency_views[2];
    H_ = dependency_views[3];
  }

  KeyVector my_keys() {
    // NOTE, a real Model would parse the parameter list to get these
    return { "A", };
  };

  KeyPairVector dependencies() {
    // NOTE, a real Model would parse the parameter list to get these
    return { {"B", ""}, {"C", ""}, {"E", ""}, {"H", ""} };
  };

  // the model
  KOKKOS_INLINE_FUNCTION void operator() (const int i) const {
    A_[i] = alpha_ * B_[i] + C_[i] * E_[i] * H_[i];
  }

  // derivatives
  //
  // NOTE: the order of these function tags (i.e. Deriv<0>, ...) is set by the
  // above call to dependencies().  Deriv<I> must correspond to the derivative
  // with respect to dependencies()[I];
  
  // d/dB
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    A_[i] = alpha_;
  }
  // d/dC
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const {
    A_[i] = E_[i] * H_[i];
  }
  // d/dE
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<2>, const int i) const {
    A_[i] = C_[i] * H_[i];
  }
  // d/dH
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<3>, const int i) const {
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
const std::string AModel<cView_type,View_type>::name("A");


template<class cView_type, class View_type>
class CModel {
 public:
  static const int n_dependencies = 2;

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

  // d/dD
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i)  const {
    C_[i] = 2.0;
  }
  // d/dG
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i)  const {
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

  KOKKOS_INLINE_FUNCTION void operator()(const int i) const {
    D_[i] = 2 * G_[i];
  }

  // derivatives
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i)  const {
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

  // d/dD
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    E_[i] = F_[i];
  }
  // d/dF
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<1>, const int i) const {
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

  // d/dG
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
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

  // d/dF
  KOKKOS_INLINE_FUNCTION void operator()(Deriv<0>, const int i) const {
    H_[i] = 2.0;
  }

  
private:
  View_type H_;
  cView_type F_;
};



// NOTE: This is tricky code (unfortunately).  This is a compile-time "loop"
// (through recursion and partial class specialization) to generate, at
// compile-time, the correct Tags needed to execute the correct derivative
// operation.  Each Launcher object's launch() call checks if wrt is dependency
// I, (starting at n_dependencies - 1), and calls Kokkos for_each if so, or
// recurses to I-1 if not.  When I-1 == -1, then wrt is not in dependencies,
// and it throws an error (ending the recursion).  Deriv<I> is the tag
// corresponding to dependency I, as indexed in the Model's list of
// dependencies.
template<int I, class Model_type, class Device_type>
class Launcher {
 public:
  Launcher(const std::string& name,
           const KeyPair& wrt,
           const KeyPairVector& dependencies,
           const Model_type& model,
           const int extent) :
      name_(name),
      wrt_(wrt),
      dependencies_(dependencies),
      model_(model),
      extent_(extent) {}

  void launch() {
    // std::cout << "wrt: " << wrt_.first << " matches " << dependencies_[I].first << " in slot " << I << "?";
    if (wrt_ == dependencies_[I]) {
      // std::cout << "  YUP, go!" << std::endl;
      Kokkos::RangePolicy<Deriv<I>, typename Device_type::execution_space> range(0, extent_);
      Kokkos::parallel_for(name_, range, model_);
    } else {
      // std::cout << "  nope, recurse!" << std::endl;
      Launcher<I-1,Model_type,Device_type> launcher(name_, wrt_, dependencies_, model_, extent_);
      launcher.launch();
    }              
  }

 private:
  const std::string& name_;
  const KeyPair& wrt_;
  const KeyPairVector& dependencies_;
  const Model_type& model_;
  const int extent_;  
};

template<class Model_type, class Device_type>
class Launcher<-1, Model_type, Device_type> {
 public:
  Launcher(const std::string& name,
           const KeyPair& wrt,
           const KeyPairVector& dependencies,
           const Model_type& model,
           const int extent)
      : name_(name),
        wrt_(wrt) {}

  void launch() {
    Errors::Message msg;
    msg << "EvaluatorModel (" << name_ << "): requested derivative with respect to ( \""
        << wrt_.first << "\",\"" << wrt_.second << "\" ) is not provided.";
    throw(msg);
  }

 private:
  const std::string& name_;
  const KeyPair& wrt_;
  
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
      auto wrt = std::make_pair(wrt_key, wrt_tag);

      Launcher<Model_type::n_dependencies-1, Model_type, Device_type> launcher(name_, wrt, dependencies_, *model_, result_views[0].extent(0));
      launcher.launch();
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
