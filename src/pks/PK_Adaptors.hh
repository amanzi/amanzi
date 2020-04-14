/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
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

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {

class Debugger;

template <class Base_t>
class PK_Adaptor : public PK, public Base_t {
 public:
  using Base_t::Base_t;

  virtual void Setup() override final { Base_t::Setup(); }
  virtual void Initialize() override final { Base_t::Initialize(); }
  virtual bool
  AdvanceStep(const Key& tag_old, const Key& tag_new) override final
  {
    return Base_t::AdvanceStep(tag_old, tag_new);
  }
  virtual bool ValidStep(const Key& tag_old, const Key& tag_new) override final
  {
    return Base_t::ValidStep(tag_old, tag_new);
  }
  virtual void CommitStep(const Key& tag_old, const Key& tag_new) override final
  {
    Base_t::CommitStep(tag_old, tag_new);
  }
  virtual void FailStep(const Key& tag_old, const Key& tag_new) override final
  {
    Base_t::FailStep(tag_old, tag_new);
  }
  virtual void CalculateDiagnostics(const Key& tag) override final
  {
    Base_t::CalculateDiagnostics(tag);
  }
  virtual void ChangedSolutionPK(const Key& tag) override final
  {
    Base_t::ChangedSolutionPK(tag);
  }
  virtual std::string name() override final { return Base_t::name(); }
  virtual double get_dt() override final { return Base_t::get_dt(); }
  virtual Teuchos::Ptr<Debugger> debugger() override final
  {
    return Base_t::debugger();
  }

 protected:
  virtual void ConstructChildren() override final
  {
    Base_t::ConstructChildren();
  }
  virtual Teuchos::RCP<TreeVectorSpace> SolutionSpace() override final
  {
    return Base_t::SolutionSpace();
  }
  virtual void StateToSolution(TreeVector& soln, const Key& tag,
                               const Key& suffix) override final
  {
    Base_t::StateToSolution(soln, tag, suffix);
  }
  virtual void SolutionToState(const Key& tag, const Key& suffix) override final
  {
    Base_t::SolutionToState(tag, suffix);
  }
  virtual void
  StateToState(const Key& tag_from, const Key& tag_to) override final
  {
    Base_t::StateToState(tag_from, tag_to);
  }

 private:
  // factory registration
  static RegisteredPKFactory<PK_Adaptor<Base_t>> reg_;

};

template <class Base_t>
class PK_Explicit_Adaptor : public PK_Explicit<>, public Base_t {
 public:
  using Base_t::Base_t;

  // the PK interface
  virtual void Setup() override final { Base_t::Setup(); }
  virtual void Initialize() override final { Base_t::Initialize(); }
  virtual bool
  AdvanceStep(const Key& tag_old, const Key& tag_new) override final
  {
    return Base_t::AdvanceStep(tag_old, tag_new);
  }
  virtual bool ValidStep(const Key& tag_old, const Key& tag_new) override final
  {
    return Base_t::ValidStep(tag_old, tag_new);
  }
  virtual void CommitStep(const Key& tag_old, const Key& tag_new) override final
  {
    Base_t::CommitStep(tag_old, tag_new);
  }
  virtual void FailStep(const Key& tag_old, const Key& tag_new) override final
  {
    Base_t::FailStep(tag_old, tag_new);
  }
  virtual void CalculateDiagnostics(const Key& tag) override final
  {
    Base_t::CalculateDiagnostics(tag);
  }
  virtual void ChangedSolutionPK(const Key& tag) override final
  {
    Base_t::ChangedSolutionPK(tag);
  }
  virtual std::string name() override final { return Base_t::name(); }
  virtual double get_dt() override final { return Base_t::get_dt(); }
  virtual Teuchos::Ptr<Debugger> debugger() override final
  {
    return Base_t::debugger();
  }

  // the Explicit_TI::fnBase interface
  virtual void FunctionalTimeDerivative(double t, const TreeVector& u,
                                        TreeVector& f) override final
  {
    Base_t::FunctionalTimeDerivative(t, u, f);
  }

  virtual void ConstructChildren() override final
  {
    Base_t::ConstructChildren();
  }
  virtual Teuchos::RCP<TreeVectorSpace> SolutionSpace() override final
  {
    return Base_t::SolutionSpace();
  }
  virtual void StateToSolution(TreeVector& soln, const Key& tag,
                               const Key& suffix) override final
  {
    Base_t::StateToSolution(soln, tag, suffix);
  }
  virtual void SolutionToState(const Key& tag, const Key& suffix) override final
  {
    Base_t::SolutionToState(tag, suffix);
  }
  virtual void
  StateToState(const Key& tag_from, const Key& tag_to) override final
  {
    Base_t::StateToState(tag_from, tag_to);
  }

 private:
  // factory registration
  static RegisteredPKFactory<PK_Explicit_Adaptor<Base_t>> reg_;
};

template <class Base_t>
class PK_Implicit_Adaptor : public PK_Implicit<TreeVector>, public Base_t {
 public:
  using Base_t::Base_t;

  // the PK interface
  virtual void Setup() override final { Base_t::Setup(); }
  virtual void Initialize() override final { Base_t::Initialize(); }
  virtual bool
  AdvanceStep(const Key& tag_old, const Key& tag_new) override final
  {
    return Base_t::AdvanceStep(tag_old, tag_new);
  }
  virtual bool ValidStep(const Key& tag_old, const Key& tag_new) override final
  {
    return Base_t::ValidStep(tag_old, tag_new);
  }
  virtual void CommitStep(const Key& tag_old, const Key& tag_new) override final
  {
    Base_t::CommitStep(tag_old, tag_new);
  }
  virtual void FailStep(const Key& tag_old, const Key& tag_new) override final
  {
    Base_t::FailStep(tag_old, tag_new);
  }
  virtual void CalculateDiagnostics(const Key& tag) override final
  {
    Base_t::CalculateDiagnostics(tag);
  }
  virtual void ChangedSolutionPK(const Key& tag) override final
  {
    Base_t::ChangedSolutionPK(tag);
  }
  virtual std::string name() override final { return Base_t::name(); }
  virtual double get_dt() override final { return Base_t::get_dt(); }
  virtual Teuchos::Ptr<Debugger> debugger() override final
  {
    return Base_t::debugger();
  }

  // the BDFfnBase interface
  // computes the non-linear functional f = f(t,u,udot)
  virtual void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new,
                     Teuchos::RCP<TreeVector> f) override final
  {
    Base_t::FunctionalResidual(t_old, t_new, u_old, u_new, f);
  }

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> Pu) override final
  {
    return Base_t::ApplyPreconditioner(u, Pu);
  }

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) override final
  {
    return Base_t::ErrorNorm(u, du);
  }

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up,
                                    double h) override final
  {
    Base_t::UpdatePreconditioner(t, up, h);
  }

  // check the admissibility of a solution
  // override final with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override final
  {
    return Base_t::IsAdmissible(up);
  }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override final
  {
    return Base_t::ModifyPredictor(h, u0, u);
  }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override final
  {
    return Base_t::ModifyCorrection(h, res, u, du);
  }

  // update the continuation parameter
  virtual void UpdateContinuationParameter(double lambda) override final
  {
    Base_t::UpdateContinuationParameter(lambda);
  }

  // calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() override final { Base_t::ChangedSolution(); }

  virtual void ConstructChildren() override final
  {
    Base_t::ConstructChildren();
  }
  virtual Teuchos::RCP<TreeVectorSpace> SolutionSpace() override final
  {
    return Base_t::SolutionSpace();
  }
  virtual void StateToSolution(TreeVector& soln, const Key& tag,
                               const Key& suffix) override final
  {
    Base_t::StateToSolution(soln, tag, suffix);
  }
  virtual void SolutionToState(const Key& tag, const Key& suffix) override final
  {
    Base_t::SolutionToState(tag, suffix);
  }
  virtual void
  StateToState(const Key& tag_from, const Key& tag_to) override final
  {
    Base_t::StateToState(tag_from, tag_to);
  }
 private:
  // factory registration
  static RegisteredPKFactory<PK_Implicit_Adaptor<Base_t>> reg_;
};

template <class Base_t>
class PK_ImplicitExplicit_Adaptor : public PK_ImplicitExplicit<TreeVector>,
                                    public Base_t {
 public:
  using Base_t::Base_t;

  // the PK interface
  virtual void Setup() override final { Base_t::Setup(); }
  virtual void Initialize() override final { Base_t::Initialize(); }
  virtual bool
  AdvanceStep(const Key& tag_old, const Key& tag_new) override final
  {
    return Base_t::AdvanceStep(tag_old, tag_new);
  }
  virtual bool ValidStep(const Key& tag_old, const Key& tag_new) override final
  {
    return Base_t::ValidStep(tag_old, tag_new);
  }
  virtual void CommitStep(const Key& tag_old, const Key& tag_new) override final
  {
    Base_t::CommitStep(tag_old, tag_new);
  }
  virtual void FailStep(const Key& tag_old, const Key& tag_new) override final
  {
    Base_t::FailStep(tag_old, tag_new);
  }
  virtual void CalculateDiagnostics(const Key& tag) override final
  {
    Base_t::CalculateDiagnostics(tag);
  }
  virtual void ChangedSolutionPK(const Key& tag) override final
  {
    Base_t::ChangedSolutionPK(tag);
  }
  virtual std::string name() override final { return Base_t::name(); }
  virtual double get_dt() override final { return Base_t::get_dt(); }
  virtual Teuchos::Ptr<Debugger> debugger() override final
  {
    return Base_t::debugger();
  }

  // the BDFfnBase interface
  // computes the non-linear functional f = f(t,u,udot)
  virtual void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new,
                     Teuchos::RCP<TreeVector> f) override final
  {
    Base_t::FunctionalResidual(t_old, t_new, u_old, u_new, f);
  }

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> Pu) override final
  {
    return Base_t::ApplyPreconditioner(u, Pu);
  }

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) override final
  {
    return Base_t::ErrorNorm(u, du);
  }

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up,
                                    double h) override final
  {
    Base_t::UpdatePreconditioner(t, up, h);
  }

  // check the admissibility of a solution
  // override final with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override final
  {
    return Base_t::IsAdmissible(up);
  }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override final
  {
    return Base_t::ModifyPredictor(h, u0, u);
  }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override final
  {
    return Base_t::ModifyCorrection(h, res, u, du);
  }

  // update the continuation parameter
  virtual void UpdateContinuationParameter(double lambda) override final
  {
    Base_t::UpdateContinuationParameter(lambda);
  }

  // calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() override final { Base_t::ChangedSolution(); }

  // the Explicit_TI::fnBase interface
  virtual void FunctionalTimeDerivative(double t, const TreeVector& u,
                                        TreeVector& f) override final
  {
    Base_t::FunctionalTimeDerivative(t, u, f);
  }

  virtual void ConstructChildren() override final
  {
    Base_t::ConstructChildren();
  }
  virtual Teuchos::RCP<TreeVectorSpace> SolutionSpace() override final
  {
    return Base_t::SolutionSpace();
  }
  virtual void StateToSolution(TreeVector& soln, const Key& tag,
                               const Key& suffix) override final
  {
    Base_t::StateToSolution(soln, tag, suffix);
  }
  virtual void SolutionToState(const Key& tag, const Key& suffix) override final
  {
    Base_t::SolutionToState(tag, suffix);
  }
  virtual void
  StateToState(const Key& tag_from, const Key& tag_to) override final
  {
    Base_t::StateToState(tag_from, tag_to);
  }

 private:
  // factory registration
  static RegisteredPKFactory<PK_ImplicitExplicit_Adaptor<Base_t>> reg_;

};

} // namespace Amanzi

#endif
