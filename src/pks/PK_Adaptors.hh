/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
*/

//! PK Adaptors map between the mixin classes and the virtual interface.

/* Developer's note:

These adaptors shouldn't need to be touched.  Specific PKs can implement some,
or all, of this interface, using or not using mixins as they desire.  They
should then register this Adaptor class with the PK_Factory.

Example:

class MyPK : public MyMixin< ... < PK_Default > > {
 public:
  void Setup() { ... }

 private:
  static RegisteredPKFactory<PK_Adaptor<MyPK> > reg_;
};

*/


#ifndef AMANZI_PK_ADAPTOR_HH_
#define AMANZI_PK_ADAPTOR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "PK.hh"

namespace Amanzi {

class Debugger;
class TreeVector;

template<class Base_t>
class PK_Adaptor : public PK,
                   public Base_t {
 public:
  using Base_t::Base_t;

  virtual void ConstructChildren() override final {
    Base_t::ConstructChildren(); }
  virtual void Setup(const TreeVector& soln) override final {
    Base_t::Setup(soln); }
  virtual void Initialize() override final {
    Base_t::Initialize(); }
  virtual bool AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                           const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    return Base_t::AdvanceStep(tag_old, soln_old, tag_new, soln_new); }
  virtual bool ValidStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                         const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    return Base_t::ValidStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void CommitStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                          const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    Base_t::CommitStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void FailStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                        const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    Base_t::FailStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void CalculateDiagnostics(const Key& tag) override final {
    Base_t::CalculateDiagnostics(tag); }
  virtual void ChangedSolutionPK(const Key& tag) override final {
    Base_t::ChangedSolutionPK(tag); }
  virtual std::string name() override final {
    return Base_t::name(); }  
  virtual double get_dt() override final {
    return Base_t::get_dt(); }
  virtual Teuchos::Ptr<Debugger> debugger() override final {
    return Base_t::debugger(); }
  virtual void StateToSolution(TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::StateToSolution(soln, tag, suffix); }
  virtual void SolutionToState(TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::SolutionToState(soln, tag, suffix); }
  virtual void SolutionToState(const TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::SolutionToState(soln, tag, suffix); }
      
  
};


template<class Base_t>
class PK_Explicit_Adaptor : public PK_Explicit<>,
                            public Base_t {
 public:

  PK_Explicit_Adaptor(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& solution) :
      Base_t(pk_tree, global_plist, S, solution) {}

  // the PK interface
  virtual void ConstructChildren() override final {
    Base_t::ConstructChildren(); }
  virtual void Setup(const TreeVector& soln) override final {
    Base_t::Setup(soln); }
  virtual void Initialize() override final {
    Base_t::Initialize(); }
  virtual bool AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                           const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    return Base_t::AdvanceStep(tag_old, soln_old, tag_new, soln_new); }
  virtual bool ValidStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                         const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    return Base_t::ValidStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void CommitStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                          const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    Base_t::CommitStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void FailStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                        const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    Base_t::FailStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void CalculateDiagnostics(const Key& tag) override final {
    Base_t::CalculateDiagnostics(tag); }
  virtual void ChangedSolutionPK(const Key& tag) override final {
    Base_t::ChangedSolutionPK(tag); }
  virtual std::string name() override final {
    return Base_t::name(); }  
  virtual double get_dt() override final {
    return Base_t::get_dt(); }
  virtual Teuchos::Ptr<Debugger> debugger() override final {
    return Base_t::debugger(); }
  virtual void StateToSolution(TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::StateToSolution(soln, tag, suffix); }
  virtual void SolutionToState(TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::SolutionToState(soln, tag, suffix); }
  virtual void SolutionToState(const TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::SolutionToState(soln, tag, suffix); }

  // the Explicit_TI::fnBase interface
  virtual void Functional(double t, const TreeVector& u, TreeVector& f) override final {
    Base_t::Functional(t,u,f); }
};


template<class Base_t>
class PK_BDF_Adaptor : public PK_BDF<>,
                            public Base_t {
 public:

  PK_BDF_Adaptor(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution) :
      Base_t(pk_tree, global_plist, S, solution) {}

  // the PK interface
  virtual void ConstructChildren() override final {
    Base_t::ConstructChildren(); }
  virtual void Setup(const TreeVector& soln) override final {
    Base_t::Setup(soln); }
  virtual void Initialize() override final {
    Base_t::Initialize(); }
  virtual bool AdvanceStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                           const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    return Base_t::AdvanceStep(tag_old, soln_old, tag_new, soln_new); }
  virtual bool ValidStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                         const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    return Base_t::ValidStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void CommitStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                          const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    Base_t::CommitStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void FailStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                        const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) override final {
    Base_t::FailStep(tag_old, soln_old, tag_new, soln_new); }
  virtual void CalculateDiagnostics(const Key& tag) override final {
    Base_t::CalculateDiagnostics(tag); }
  virtual void ChangedSolutionPK(const Key& tag) override final {
    Base_t::ChangedSolutionPK(tag); }
  virtual std::string name() override final {
    return Base_t::name(); }  
  virtual double get_dt() override final {
    return Base_t::get_dt(); }
  virtual Teuchos::Ptr<Debugger> debugger() override final {
    return Base_t::debugger(); }

  // the BDFfnBase interface
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) override final {
    Base_t::Functional(t_old, t_new, u_old, u_new, f);
  }

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override final {
    return Base_t::ApplyPreconditioner(u, Pu);
  }

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override final {
    return Base_t::ErrorNorm(u, du);
  }

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override final {
    Base_t::UpdatePreconditioner(t, up, h);
  }

  // check the admissibility of a solution
  // override final with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override final {
    return Base_t::IsAdmissible(up);
  }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override final {
    return Base_t::ModifyPredictor(h, u0, u);
  }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult 
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) override final {
    return Base_t::ModifyCorrection(h, res, u, du);
  }

  // update the continuation parameter
  virtual void UpdateContinuationParameter(double lambda) override final {
    Base_t::UpdateContinuationParameter(lambda);
  }

  // calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() override final {
    Base_t::ChangedSolution();
  }

  // experimental routine -- returns the number of linear iterations.
  virtual int ReportStatistics() override final {
    return Base_t::ReportStatistics();
  }

  virtual void StateToSolution(TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::StateToSolution(soln, tag, suffix);
  }
  virtual void SolutionToState(TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::SolutionToState(soln, tag, suffix);
  }
  virtual void SolutionToState(const TreeVector& soln, const Key& tag, const Key& suffix) override final {
    Base_t::SolutionToState(soln, tag, suffix);
  }
  

};

  




}  // namespace Amanzi

#endif
