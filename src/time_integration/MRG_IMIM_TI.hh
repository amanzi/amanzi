
/*
This class implements Implicit-Implicit Multirate Time Integration Methods:

TODO: expand documentation
<List of methods>

 - User Defined


NOTE: Can only support TreeVectors

*/

#ifndef AMANZI_MRG_IMIM_TI_HH_
#define AMANZI_MRG_IMIM_TI_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include <typeinfo>
#include <functional>
#include "errors.hh"
#include "dbc.hh"

#include "VerboseObject.hh"
#include "Solver.hh"
#include "SolverFactory.hh"
#include "SolverDefs.hh"
#include "SolverFnBase.hh"


#include "PartitionFnBase.hh"
#include "MRG_IMIM_FnBase.hh"
#include "MRG_IMIM_SolverFnBase_Fast.hh"
#include "MRG_IMIM_SolverFnBase_Full.hh"

namespace Amanzi
{

  /**
   * @brief Method enums
   * TODO: Add Names and an order 3 method
   */
  enum IMIM_method_t
  {
    MrGark_IMIM_2,
    MrGark_IMIM_temp3,
    MrGark_IMIM_user_defiend
  };

  /**
   * @brief Class for the IMIM MR GARK Time Integration method
   * 
   * Note: Vector is need to support easy concatination of vectors for full solve
   * 
   * @tparam Vector 
   * @tparam VectorSpace 
   */
  template <class Vector, class VectorSpace>
  class MRG_IMIM_TI
  {
  protected:
    Teuchos::RCP<MRG_IMIM_FnBase<Vector>> fn_;
    
    IMIM_method_t method_;

    //FIXME: This should give a Linking error
    Teuchos::RCP<AmanziSolvers::Solver<Vector,VectorSpace> > solverfull_;
    Teuchos::RCP<MRG_IMIM_SolverFnBase_Full<Vector> > solver_fn_full_;
    Teuchos::ParameterList plist_full_;

    Teuchos::RCP<AmanziSolvers::Solver<Vector,VectorSpace> > solverfast_;
    Teuchos::RCP<MRG_IMIM_SolverFnBase_Fast<Vector> > solver_fn_fast_;
    Teuchos::ParameterList plist_fast_;

    //FIXME: Set these up
    Teuchos::RCP<VerboseObject> vo_fast_;
    Teuchos::RCP<AmanziSolvers::ResidualDebugger> db_fast_;

    Teuchos::RCP<VerboseObject> vo_full_;
    Teuchos::RCP<AmanziSolvers::ResidualDebugger> db_full_;

    int order_;
    int stage_;
    std::vector<double> a_fastfast_, b_fast_, c_fast_;
    std::vector<double> a_slowslow_, b_slow_, c_slow_;
    // Both functions for coupling will take M and lambda respectively
    std::function<void(double, double, std::vector<double>&)> a_slowfast_;
    std::function<void(double, double, std::vector<double>&)> a_fastslow_;

    // Allocate the coupling matrices based of M and update when a new M is given
    std::vector<double> a_slowfast_local_;
    std::vector<double> a_fastslow_local_;

    //Scalings for the coupled nonlinear solver
    std::vector<double> scaling_non = {0.0,0.0,0.0,0.0};

    std::vector<Teuchos::RCP<Vector>> k_fast_;
    std::vector<Teuchos::RCP<Vector>> k_slow_;
    std::vector<Teuchos::RCP<Vector>> k_slowfast_;

    //Full Vectors for the nonlinear solver
    //Data used through the algorithim
    Teuchos::RCP<Vector> y_n_full_;
    Teuchos::RCP<Vector> y_n_full_old_;
    Teuchos::RCP<Vector> y_exp_full_;
    Teuchos::RCP<Vector> y_lambda_;

    //Each points to the respective subvector of the fulls
    //Only Views
    Teuchos::RCP<Vector> y_exp_f_;
    Teuchos::RCP<Vector> y_exp_s_;
    Teuchos::RCP<Vector> y_ns;
    Teuchos::RCP<Vector> y_nf;
    Teuchos::RCP<Vector> y_ns_old;
    Teuchos::RCP<Vector> y_nf_old;


    void InitMethod(const IMIM_method_t method);
    void InitMemory(const Teuchos::RCP<Vector> initvector_fast, const Teuchos::RCP<Vector> initvector_full);
    void InitCoefficentMemory();
    void InitSolvers(Teuchos::ParameterList& plist_fast,Teuchos::ParameterList& plist_full , const Teuchos::RCP<Vector> initvector_fast, const Teuchos::RCP<Vector> initvector_full);
    int SolveNonlinearSystemFast(double t_n1, double t_n0, double scaling, Teuchos::RCP<Vector> y_old, const Teuchos::RCP<Vector>& y_new);
    int SolveNonlinearSystemFull(double t_n1, double t_n0, std::vector<double> scalings, Teuchos::RCP<Vector> y_old, const Teuchos::RCP<Vector>& y_new);
    
  public:
    MRG_IMIM_TI(MRG_IMIM_FnBase<Vector> &fn,
               const IMIM_method_t method,
               Teuchos::ParameterList& plist_full,
               Teuchos::ParameterList& plist_fast,
                const Teuchos::RCP<Vector> initvector_fast,
                const Teuchos::RCP<Vector> initvector_full);

    MRG_IMIM_TI(MRG_IMIM_FnBase<Vector> &fn,
               const int order,
               const int stage,
               const std::vector<double> a_fastfast,
               const std::vector<double> b_fast,
               const std::vector<double> c_fast,
               const std::vector<double> a_slowslow,
               const std::vector<double> b_slow,
               const std::vector<double> c_slow,
               const std::function<void(double, double, std::vector<double>&)> a_slowfast_,
               const std::function<void(double, double, std::vector<double>&)> a_fastslow_,
               Teuchos::ParameterList& plist_full,
               Teuchos::ParameterList& plist_fast,
                const Teuchos::RCP<Vector> initvector_fast,
                const Teuchos::RCP<Vector> initvector_full);

    int TimeStep(const double t, const double h, const int M, const Teuchos::RCP<Vector> y, Teuchos::RCP<Vector> y_new);

    int order() { return order_; };
    int stage() { return stage_; };
  };

  template <class Vector, class VectorSpace>
  void MRG_IMIM_TI<Vector, VectorSpace>::InitMethod(const IMIM_method_t method)
  {

    switch (method)
    {
    case MrGark_IMIM_2:
      // Order 2 method
      {
      double alpha = 1.0 - 1.0/sqrt(2.0);

      a_fastfast_ = {alpha, 0, 1.0 - alpha, alpha};
      b_fast_= {1.0 - alpha, alpha};
      c_fast_ = {alpha, 1.0};

      a_slowslow_ = {alpha, 0, 1.0 - alpha, alpha};
      b_slow_= {1.0 - alpha, alpha};
      c_slow_ = {alpha, 1.0};

      // The conditional coupling coefficents
      a_fastslow_ = [alpha] ( double M, double lambda, std::vector<double>& a){
        a[0] = (lambda - 1.0 + alpha) / M;
        if (lambda == M)
        {
          a[2] = 1.0 - alpha;
          a[3] = alpha;
        }
        else
        {
          a[2] = lambda / M;
          a[3] = 0.0;
        }
      };

      a_slowfast_ = [alpha] ( double M, double lambda, std::vector<double>& a){
        a[0] = (lambda == 1.0) ? M * alpha : 0.0;
        a[2] = (1.0 - alpha);
        a[3] = alpha;
      };

      //TEST: These are testing coefficents. TODO: Remove bad ones.

      // // The conditional coupling coefficents
      // a_fastslow_ = [alpha] ( double M, double lambda, std::vector<double>& a){
      //   a[0] = alpha * (2.0*lambda - 1.0) / M;
      //   a[1] = (lambda - 1.0)*(1.0 - 2.0*alpha)/M;
      //   a[2] = (1.0 - 3.0*alpha + 2.0 * lambda * alpha) / M;
      //   a[3] = (3.0*alpha - 1.0 + (1.0 - 2.0 * alpha)*lambda)/M;
      // };

      // // The conditional coupling coefficents
      // a_fastslow_ = [alpha, this] ( double M, double lambda, std::vector<double>& a){
      //   if (lambda == 1.0)
      //   {
      //     a = a_fastfast_;
      //   }
      //   else
      //   {
      //     double denominator = M * (alpha - 1);
      //     a[0] = -alpha*((M-2.0)*alpha + 3.0) + (2.0*alpha-1.0)*lambda + 1.0 / denominator;
      //     a[1] = alpha*((M - 1.0)* alpha - lambda + 1.0) / denominator;
      //     a[2] = (-M * alpha * alpha + 2.0 * lambda * alpha - lambda) / denominator;
      //     a[3] = alpha * (M * alpha - lambda) / denominator;
      //   }
      // };

      // a_slowfast_ = [alpha, this] ( double M, double lambda, std::vector<double>& a){
      //   if (lambda == 1.0)
      //   {
      //     a = a_fastfast_;
      //   }
      //   else
      //   {
      //     a[0] = 0;
      //     a[1] = 0;
      //     a[2] = 0;
      //     a[3] = 0;
      //   }        
      // };
      }

      
      stage_ = 2;
      order_ = 2;


      break;

    case MrGark_IMIM_temp3:
      // Order 3 method
      stage_ = 3;
      order_ = 3;

      // TODO: get coefficents
      break;
    default:
      // TODO: warn user of improper method
      stage_ = c_fast_.size();
      order_ = 0;

      break;
    }
  }

  /**
   * @brief Intalize memory for coefficents
   * 
   * Memory for local storage of coefficents
   * 
   * @tparam Vector 
   * @tparam VectorSpace 
   */
template <class Vector, class VectorSpace>
  void MRG_IMIM_TI<Vector, VectorSpace>::InitCoefficentMemory()
  {
    a_slowfast_local_.resize(stage_ * stage_);
    a_fastslow_local_.resize(stage_ * stage_);
  }

/**
 * @brief Intalize memory for the internal stage evaluations 
 * 
 * @tparam Vector 
 * @tparam VectorSpace 
 * @param initvector 
 */
  template <class Vector, class VectorSpace>
  void MRG_IMIM_TI<Vector, VectorSpace>::InitMemory(const Teuchos::RCP<Vector> initvector_fast, const Teuchos::RCP<Vector> initvector_full)
  {
    k_slow_.resize(stage_);
    k_fast_.resize(stage_);
    k_slowfast_.resize(stage_);
    for (int i = 0; i < stage_; ++i)
    {
      k_slow_[i] = Teuchos::rcp(new Vector(*initvector_fast));
      k_fast_[i] = Teuchos::rcp(new Vector(*initvector_fast));
      k_slowfast_[i] = Teuchos::rcp(new Vector(*initvector_fast));
    }

    //nonlinear Solver Vectors
    //Data Vectors
    //Each have two entries
    y_exp_full_ = Teuchos::rcp(new Vector(*initvector_full));
    y_n_full_ = Teuchos::rcp(new Vector(*initvector_full));
    y_n_full_old_ = Teuchos::rcp(new Vector(*initvector_full));

    y_lambda_ = Teuchos::rcp(new Vector(*initvector_fast));

    //Reference Pointers to view the sub tree subvectors
    y_exp_f_ = y_exp_full_->SubVector(0);
    y_exp_s_ = y_exp_full_->SubVector(1);

    y_nf = y_n_full_->SubVector(0);
    y_ns = y_n_full_->SubVector(1);
    
    y_nf_old = y_n_full_old_->SubVector(0);
    y_ns_old = y_n_full_old_->SubVector(1);
    
  }
  
  template <class Vector, class VectorSpace>
  void MRG_IMIM_TI<Vector, VectorSpace>::InitSolvers(Teuchos::ParameterList& plist_fast,Teuchos::ParameterList& plist_full , const Teuchos::RCP<Vector> initvector_fast, const Teuchos::RCP<Vector> initvector_full)
  {
    
    AmanziSolvers::SolverFactory<Vector,VectorSpace> factory;

    // update the verbose options
    vo_fast_ = Teuchos::rcp(new VerboseObject(initvector_fast->Comm(), "TI::MRG_IMIM_Fast", plist_fast_));
    db_fast_ = Teuchos::rcp(new AmanziSolvers::ResidualDebugger(plist_fast_.sublist("residual debugger")));

    // Set up the nonlinear solver
    // -- initialized the SolverFnBase interface
    solver_fn_fast_ = Teuchos::rcp(new MRG_IMIM_SolverFnBase_Fast<Vector>(plist_fast_, fn_));

    solverfast_ = factory.Create(plist_fast_);
    solverfast_->set_db(db_fast_);
    solverfast_->Init(solver_fn_fast_, initvector_fast->Map());

    
    AmanziSolvers::SolverFactory<Vector,VectorSpace> factory_full;

    // update the verbose options
    vo_full_ = Teuchos::rcp(new VerboseObject(initvector_full->Comm(), "TI::MRG_IMIM_Full", plist_full_));
    db_full_ = Teuchos::rcp(new AmanziSolvers::ResidualDebugger(plist_full_.sublist("residual debugger")));

    // Set up the nonlinear solver
    // -- initialized the SolverFnBase interface
    solver_fn_full_ = Teuchos::rcp(new MRG_IMIM_SolverFnBase_Full<Vector>(plist_full_, fn_));

    solverfull_ = factory_full.Create(plist_full_);
    solverfull_->set_db(db_full_);
    solverfull_->Init(solver_fn_full_, initvector_full->Map());
  }
    

  template <class Vector, class VectorSpace>
  MRG_IMIM_TI<Vector, VectorSpace>::MRG_IMIM_TI(MRG_IMIM_FnBase<Vector> &fn,
                                              const IMIM_method_t method,
                                              Teuchos::ParameterList& plist_full,
                                              Teuchos::ParameterList& plist_fast,
                                              const Teuchos::RCP<Vector> initvector_fast,
                                              const Teuchos::RCP<Vector> initvector_full) : plist_fast_(plist_fast), plist_full_(plist_full)
  {
    
    fn_ = Teuchos::rcpFromRef(fn);

    MRG_IMIM_TI<Vector, VectorSpace>::InitMethod(method);
    MRG_IMIM_TI<Vector, VectorSpace>::InitMemory(initvector_fast, initvector_full);
    MRG_IMIM_TI<Vector, VectorSpace>::InitCoefficentMemory();
    MRG_IMIM_TI<Vector, VectorSpace>::InitSolvers(plist_fast_, plist_full_, initvector_fast, initvector_full);
  }

  template <class Vector, class VectorSpace>
  MRG_IMIM_TI<Vector, VectorSpace>::MRG_IMIM_TI(MRG_IMIM_FnBase<Vector> &fn,
                                              const int order,
                                              const int stage,
                                              const std::vector<double> a_fastfast,
                                              const std::vector<double> b_fast,
                                              const std::vector<double> c_fast,
                                              const std::vector<double> a_slowslow,
                                              const std::vector<double> b_slow,
                                              const std::vector<double> c_slow,
                                              const std::function<void(double, double, std::vector<double>&)> a_slowfast,
                                              const std::function<void(double, double, std::vector<double>&)> a_fastslow,
                                              Teuchos::ParameterList& plist_full,
                                              Teuchos::ParameterList& plist_fast,
                                              const Teuchos::RCP<Vector> initvector_fast,
                                              const Teuchos::RCP<Vector> initvector_full) : order_(order), stage_(stage), method_(MrGark_IMIM_user_defiend),
                                                                          a_slowfast_(a_slowfast), a_fastslow_(a_fastslow), plist_fast_(plist_fast), plist_full_(plist_full)
  {
    
    fn_ = Teuchos::rcpFromRef(fn);

    b_fast_.resize(b_fast.size());
    c_fast_.resize(c_fast.size());
    b_slow_.resize(b_slow.size());
    c_slow_.resize(c_slow.size());

    b_fast_ = b_fast;
    c_fast_ = c_fast;
    b_slow_ = b_slow;
    c_slow_ = c_slow;

    a_fastfast_.resize(a_fastfast.size());
    a_slowslow_.resize(a_slowslow.size());
    a_fastfast_ = a_fastfast;
    a_slowslow_ = a_slowslow;

    stage_ = c_fast_.size();

    MRG_IMIM_TI<Vector, VectorSpace>::InitMemory(initvector_fast, initvector_full);
    MRG_IMIM_TI<Vector, VectorSpace>::InitCoefficentMemory();
    MRG_IMIM_TI<Vector, VectorSpace>::InitSolvers(plist_fast_, plist_full_, initvector_fast, initvector_full);
  
  }


  template <class Vector, class VectorSpace>
  int MRG_IMIM_TI<Vector, VectorSpace>::TimeStep(const double t, const double h, const int M, const Teuchos::RCP<Vector> y, Teuchos::RCP<Vector> y_new)
  {

    // Could simplify with templates and commands as many operations repeat

    double sum_time = t;
    double sum_time_old = t;
    double M_cast = static_cast<double>(M);

    *y_lambda_ = *y;
    *y_ns_old = *y;
    *y_nf_old = *y;

    int q_slow = 0;
    double h_fast = h / static_cast<double>(M);

    for (int i = 0; i < M; ++i)
    {
      // intialize the coupling coefficents
      a_fastslow_(M_cast, static_cast<double>(i + 1), a_fastslow_local_);
      a_slowfast_(M_cast, static_cast<double>(i + 1), a_slowfast_local_ );

      for (int j = 0; j < stage_; ++j)
      {
        sum_time = t + (static_cast<double>(i) + c_fast_[j]) * h_fast;

        y_exp_f_->Scale(0.0);

        // Add coupled slow stages
        for (int k = 0; k < q_slow; ++k)
        {
          int index = j * stage_ + k;
          if (a_fastslow_local_[index] != 0)
          {
            y_exp_f_->Update(a_fastslow_local_[index], *k_slow_[k], 1.0);
          }
        }

        // Add fast stages
        for (int k = 0; k < j; ++k)
        {
          int index = j * stage_ + k;
          if (a_fastfast_[index] != 0)
          {
            y_exp_f_->Update(a_fastfast_[index], *k_fast_[k], 1.0);
          }
        }

        //Assume there are no slow only stages
        if (q_slow < stage_ && a_fastslow_local_[j*stage_ + q_slow] != 0 && a_slowfast_local_[q_slow * stage_ + j])
        {
          y_exp_s_->Scale(0.0);

          for (int k = 0; k < q_slow; k++)
          {
            int index = q_slow * stage_ + k;
            if (a_slowslow_[index] != 0)
            {
              y_exp_s_->Update(a_slowslow_[index], *k_slow_[k], 1.0);
            }
          }

          for (int k = 0; k < j; k++)
          {
            int index = q_slow * stage_ + k;
            if (a_slowfast_local_[index] != 0)
            {
              y_exp_s_->Update(a_slowfast_local_[index], *k_fast_[k], 1.0);
            }
          }
          if (i > 0)
          {
            y_exp_s_->Update(1.0, *k_slowfast_[q_slow], 1.0);
          }

          scaling_non[0] = h_fast * a_fastfast_[j * (stage_ + 1)];
          scaling_non[1] = h * a_fastslow_local_[j * stage_ + q_slow];
          scaling_non[2] = h_fast * a_slowfast_local_[q_slow * stage_ + j];
          scaling_non[3] = h * a_slowslow_[q_slow * (stage_ + 1)];

          int flag_full = SolveNonlinearSystemFull(sum_time, sum_time_old, scaling_non, y_n_full_old_, y_n_full_);

          if (flag_full)
          {
            //TODO: add flag based operations
          }

          fn_->ModifySolutionSlow(sum_time, y_ns);
          fn_->FunctionalTimeDerivativeSlow(sum_time, y_ns, k_slow_[q_slow]);
          k_slow_[q_slow]->Scale(h);

          q_slow++;
        }
        else
        {
          int flag_fast = SolveNonlinearSystemFast(sum_time, sum_time_old, h_fast*a_fastfast_[j * (stage_ + 1)], y_nf_old,  y_nf);

          if (flag_fast)
          {
            //TODO: add flag based operations
          }
          
        }

        fn_->ModifySolutionFast(sum_time, y_nf);
        fn_->FunctionalTimeDerivativeFast(sum_time, y_nf, k_fast_[j]);
        k_fast_[j]->Scale(h_fast);

        sum_time_old = sum_time;

      }

      //Accumalate y_lambda for next fast macro step
      for (int j = 0; j < stage_; j++)
      {
        if (b_fast_[j] != 0)
        {
          y_lambda_->Update(b_fast_[j], *k_fast_[j], 1.0);
        }
        for (int k = q_slow; k < stage_; k++)
        {
          int index = j*stage_ + k;
          if (a_slowfast_local_[index] != 0)
          {
            k_slowfast_[j]->Update(a_slowfast_local_[index], *k_fast_[k], (i > 0) ? 1.0 : 0.0);
          }
        }
      }

      *y_nf_old = *y_lambda_;
    }

    //Assume that all slows are coupled with one fast so all stages are solved unlike EXIM

    *y_new = *y_lambda_;
    for (int i = 0; i < stage_; i++)
    {
      if (b_slow_[i] != 0)
      {
        y_new->Update(b_slow_[i], *k_slow_[i], 1.0);
      }
    }

    return 0;
  }


  
  //TODO: Documentation
  template <class Vector, class VectorSpace>
  int MRG_IMIM_TI<Vector, VectorSpace>::SolveNonlinearSystemFast(double t_n1, double t_n0, double scaling, Teuchos::RCP<Vector> y_old, const Teuchos::RCP<Vector>& y_new)
  {
    db_fast_->StartIteration<VectorSpace>(t_n0, 0, 1, y_new->Map());

    *y_new = *y_old;
    
    solver_fn_fast_->SetTimes(t_n0, t_n1);
    solver_fn_fast_->SetPreviousTimeSolution(y_old);
    solver_fn_fast_->SetExplicitTerms(y_exp_f_);
    solver_fn_fast_->SetScaling(scaling);

    // Solve the nonlinear system.
    int ierr, code, itr;
    try {
      ierr = solverfast_->Solve(y_new);
      itr = solverfast_->num_itrs();
      code = solverfast_->returned_code();
    } catch (const Errors::CutTimeStep& e) {
      ierr = 1;
      itr = -1;  // This should not be summed up into the global counter.
      code = AmanziSolvers::SOLVER_INTERNAL_EXCEPTION;
      if (vo_fast_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_fast_->os() << e.what() << std::endl;
      }
    }
    
    if (ierr == 0) {
      if (vo_fast_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_fast_->os() << "success: " << itr << " nonlinear Fast itrs" 
                  << " error=" << solverfast_->residual() << std::endl;
      }
    } else {
      if (vo_fast_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_fast_->os() << vo_fast_->color("red") << "step failed with error code " << code << vo_fast_->reset() << std::endl;
      }
    } 


    return ierr;
  }

  //TODO: Documentation
  template <class Vector, class VectorSpace>
  int MRG_IMIM_TI<Vector, VectorSpace>::SolveNonlinearSystemFull(double t_n1, double t_n0, std::vector<double> scalings, Teuchos::RCP<Vector> y_old, const Teuchos::RCP<Vector>& y_new)
  {
    db_full_->StartIteration<VectorSpace>(t_n0, 0, 1, y_new->Map());

    *y_new = *y_old;
    
    solver_fn_full_->SetTimes(t_n0, t_n1);
    solver_fn_full_->SetPreviousTimeSolution(y_old);
    solver_fn_full_->SetExplicitTerms(y_exp_full_);
    solver_fn_full_->SetScaling(scalings);

      // Solve the nonlinear system.
    int ierr, code, itr;
    try {
      ierr = solverfull_->Solve(y_new);
      itr = solverfull_->num_itrs();
      code = solverfull_->returned_code();
    } catch (const Errors::CutTimeStep& e) {
      ierr = 1;
      itr = -1;  // This should not be summed up into the global counter.
      code = AmanziSolvers::SOLVER_INTERNAL_EXCEPTION;
      if (vo_full_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_full_->os() << e.what() << std::endl;
      }
    }
    
    if (ierr == 0) {
      if (vo_full_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_full_->os() << "success: " << itr << " nonlinear itrs" 
                  << " error=" << solverfull_->residual() << std::endl;
      }
    } else {
      if (vo_full_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_full_->os() << vo_full_->color("red") << "step failed with error code " << code << vo_full_->reset() << std::endl;
      }
    } 


    return ierr;
  }

}

#endif