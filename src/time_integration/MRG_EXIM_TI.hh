
/*
This class implements Implicit-Explicit Multirate Time Integration Methods:

TODO: expand documentation
<List of methods>

 - User Defined

 



*/

#ifndef AMANZI_MRG_EXIM_TI_HH_
#define AMANZI_MRG_EXIM_TI_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Teuchos_RCP.hpp"
#include <functional>
#include "errors.hh"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "dbc.hh"

#include "PartitionFnBase.hh"
#include "MRG_EXIM_FnBase.hh"
#include "MRG_EXIM_SolverFnBase.hh"
#include "VerboseObject.hh"


#include "Solver.hh"
#include "SolverFactory.hh"
#include "SolverDefs.hh"

namespace Amanzi
{

  /**
   * @brief Method enums
   * TODO: Add appropiate names and order 3 coefficents
   */
  enum method_t
  {
    MrGark_EXIM_2,
    temp3,
    user_defiend
  };

  template <class Vector, class VectorSpace>
  class MRG_EXIM_TI
  {
  private:
    //RHS
    Teuchos::RCP<MRG_EXIM_FnBase<Vector>> fn_;

    method_t method_;

    // Solver Parameters
    Teuchos::RCP<AmanziSolvers::Solver<Vector,VectorSpace> > solverslow_;
    Teuchos::RCP<MRG_EXIM_SolverFnBase<Vector> > solver_fn_slow_;
    Teuchos::ParameterList plist_;

    //FIXME: Set these up
    Teuchos::RCP<VerboseObject> vo_;
    Teuchos::RCP<AmanziSolvers::ResidualDebugger> db_;

    //Method parameters
    int order_;
    int stage_;
    std::vector<double> a_fastfast_, b_fast_, c_fast_;
    std::vector<double> a_slowslow_, b_slow_, c_slow_;

    // Both functions for coupling will take M and lambda respectively
    std::function<void(double, double, std::vector<double>&)> a_slowfast_;
    std::function<void(double, double, std::vector<double>&)> a_fastslow_;

    // Allocate the coupling matrices based of M and update when a new M is given
    int prev_M_ = 0;
    std::vector<double> a_fastfast_local_, b_fast_local_, c_fast_local_;
    std::vector<double> a_slowfast_local_;
    std::vector<double> a_fastslow_local_;

    //Internal Stage Derivatives
    std::vector<Teuchos::RCP<Vector>> k_fast_;
    std::vector<Teuchos::RCP<Vector>> k_slow_;

    //Extra Storage to prevent reallocation
    Teuchos::RCP<Vector> ya_slowfast_;
    Teuchos::RCP<Vector> y_lambda_;
    Teuchos::RCP<Vector> y_ns;
    Teuchos::RCP<Vector> y_ns_old;
    Teuchos::RCP<Vector> y_exp_f_;
    Teuchos::RCP<Vector> y_exp_s_;

    void InitMethod(const method_t method);
    void InitMemory(const Teuchos::RCP<Vector> initvector);
    void InitCoefficentMemory();
    void InitSolvers(Teuchos::ParameterList& plist, const Teuchos::RCP<Vector> initvector);
    int SolveNonlinearSystem(double t_n1, double t_n0, double h, double scaling, Teuchos::RCP<Vector> y_old, const Teuchos::RCP<Vector>& y_new);

  public:

    MRG_EXIM_TI(MRG_EXIM_FnBase<Vector> &fn,
                const method_t method,
                Teuchos::ParameterList& plist,
                const Teuchos::RCP<Vector> initvector);

    MRG_EXIM_TI(MRG_EXIM_FnBase<Vector> &fn,
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
               Teuchos::ParameterList& plist,
               const Teuchos::RCP<Vector> initvector);

    /**
     * @brief Time step procedure
     * 
     * Will advance the solution by h amount but will advance the fast solution over more steps
     * 
     * @param t current time
     * @param h slow step
     * @param M ratio: M = h_s/h_f
     * @param y Current solution
     * @param y_new Computed solution
     * @return int The flag on any problems
     *          0: The time step was a success
     *          ow: TODO: Need to define some flags
     */
    int TimeStep(const double t, const double h, const int M, const Teuchos::RCP<Vector> y, Teuchos::RCP<Vector> y_new);

    int order() { return order_; };
    int stage() { return stage_; };
  };

  template <class Vector, class VectorSpace>
  void MRG_EXIM_TI<Vector, VectorSpace>::InitMethod(const method_t method)
  {

    switch (method)
    {
    case MrGark_EXIM_2:
      {
      // Order 2 method
      double invsqrt_2 = 1.0/sqrt(2.0);

      a_fastfast_ = {0, 0, 2.0 / 3.0, 0};
      b_fast_= {1.0/4.0, 3.0/4.0};
      c_fast_ = {0, 2.0/3.0};

      a_slowslow_ = {1 - invsqrt_2, 0, invsqrt_2, 1 - invsqrt_2};
      b_slow_= {invsqrt_2, 1 - invsqrt_2};
      c_slow_ = {1 - invsqrt_2, 1};

      // The conditional coupling coefficents
      a_fastslow_ = [=] ( double M, double lambda, std::vector<double>& a){
        a[0] = (lambda - 1) / M;

        a[2] = (3.0*lambda - 1) / ( 3.0 * M);
      };

      a_slowfast_ = [invsqrt_2] ( double M, double lambda, std::vector<double>& a){
        a[0] = (lambda == 1) ? (M - M * invsqrt_2) : 0;
        a[2] = 1.0/4.0;
        a[3] = 3.0/4.0;
      };
      }
      
      stage_ = 2;
      order_ = 2;


      break;

    case temp3:
      // Order 3 method
      stage_ = 3;
      order_ = 3;

      // TODO: input coefficents
      break;
    default:
      // TODO:
      
      stage_ = -1;
      order_ = -1;

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
  void MRG_EXIM_TI<Vector, VectorSpace>::InitCoefficentMemory() {
    a_fastfast_local_.resize(stage_ * stage_);
    b_fast_local_.resize(stage_);
    c_fast_local_.resize(stage_);
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
  void MRG_EXIM_TI<Vector, VectorSpace>::InitMemory(const Teuchos::RCP<Vector> initvector)
  {
    k_slow_.resize(stage_);
    k_fast_.resize(stage_);
    for (int i = 0; i < stage_; ++i)
    {
      k_slow_[i] = Teuchos::rcp(new Vector(*initvector));
      k_fast_[i] = Teuchos::rcp(new Vector(*initvector));
    }

    ya_slowfast_ = Teuchos::rcp(new Vector(*initvector));

    y_lambda_ = Teuchos::rcp(new Vector(*initvector));
    y_exp_f_ = Teuchos::rcp(new Vector(*initvector));
    y_exp_s_ = Teuchos::rcp(new Vector(*initvector));
    y_ns = Teuchos::rcp(new Vector(*initvector));
    y_ns_old = Teuchos::rcp(new Vector(*initvector));
  }

  /**
   * @brief Initalize the nonlinear solvers
   * 
   * @tparam Vector 
   * @tparam VectorSpace the map
   * @param plist the parameters for the nonlinear solver
   * @param initvector vector for intial sizeing
   */
  template<class Vector, class VectorSpace>
  void MRG_EXIM_TI<Vector, VectorSpace>::InitSolvers(Teuchos::ParameterList& plist, const Teuchos::RCP<Vector> initvector){
    
    // update the verbose options
    vo_ = Teuchos::rcp(new VerboseObject(initvector->Comm(), "TI::MRG_EXIM", plist));
    db_ = Teuchos::rcp(new AmanziSolvers::ResidualDebugger(plist.sublist("residual debugger")));

    // Set up the nonlinear solver
    // -- initialized the SolverFnBase interface
    solver_fn_slow_ = Teuchos::rcp(new MRG_EXIM_SolverFnBase<Vector>(plist_, fn_));

    //FIXME: Linkage Errors from lines below

    AmanziSolvers::SolverFactory<Vector,VectorSpace> factory;
    solverslow_ =  factory.Create(plist_);
    solverslow_->set_db(db_);
    solverslow_->Init(solver_fn_slow_, initvector->Map());
  }



  template <class Vector, class VectorSpace>
  MRG_EXIM_TI<Vector, VectorSpace>::MRG_EXIM_TI(MRG_EXIM_FnBase<Vector> &fn,
                                              const method_t method,
                                              Teuchos::ParameterList& plist,
                                              const Teuchos::RCP<Vector> initvector) : plist_(plist)
  {
    fn_ = Teuchos::rcpFromRef(fn);

    MRG_EXIM_TI<Vector,VectorSpace>::InitMethod(method);
    MRG_EXIM_TI<Vector,VectorSpace>::InitMemory(initvector);
    MRG_EXIM_TI<Vector, VectorSpace>::InitCoefficentMemory();
    MRG_EXIM_TI<Vector, VectorSpace>::InitSolvers( plist, initvector);

  }

  template <class Vector, class VectorSpace>
  MRG_EXIM_TI<Vector, VectorSpace>::MRG_EXIM_TI(MRG_EXIM_FnBase<Vector> &fn,
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
                                              Teuchos::ParameterList& plist,
                                              const Teuchos::RCP<Vector> initvector) : order_(order), stage_(stage), method_(user_defiend),
                                                                          a_slowfast_(a_slowfast), a_fastslow_(a_fastslow), plist_(plist)
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

    MRG_EXIM_TI<Vector, VectorSpace>::InitMemory(initvector);
    MRG_EXIM_TI<Vector, VectorSpace>::InitCoefficentMemory();
    MRG_EXIM_TI<Vector, VectorSpace>::InitSolvers( plist, initvector);
  }


  /**
   * @brief Step f in time
   * 
   * @tparam Vector 
   * @tparam VectorSpace 
   * @param t 
   * @param h the slow step
   * @param M the ratio between the slow and fast steps
   * @param y 
   * @param y_new output
   */
  template <class Vector, class VectorSpace>
  int MRG_EXIM_TI<Vector, VectorSpace>::TimeStep(const double t, const double h, const int M, const Teuchos::RCP<Vector> y, Teuchos::RCP<Vector> y_new)
  {

    double sum_timef = t;

    *y_lambda_ = *y;
    int q_slow = 0;
    double h_fast = h / static_cast<double>(M);

    double slow_time_old = t;
    double slow_time = t;
    
    double M_cast = static_cast<double>(M);

    // Intalize coefficents for scaled M only when it is different
    if (M != prev_M_)
    {
      a_fastfast_local_ = a_fastfast_;
      for (int i = 0; i < a_fastfast_local_.size(); i++)
      {
        a_fastfast_local_[i] /= M_cast;
      }
      b_fast_local_ = b_fast_;
      c_fast_local_ = c_fast_;
      for (int i = 0; i < stage_; i++)
      {
        b_fast_local_[i] /= M_cast;
        c_fast_local_[i] /= M_cast;
      }

      prev_M_ = M;
    }

    for (int i = 0; i < M; ++i)
    {
      // intialize the coupling coefficents
      a_fastslow_(M_cast, static_cast<double>(i + 1), a_fastslow_local_);
      a_slowfast_(M_cast, static_cast<double>(i + 1), a_slowfast_local_);


      for (int j = 0; j < stage_; ++j)
      {
        sum_timef = t + ((static_cast<double>(i) + c_fast_[j])/M_cast) * h_fast;

        //Check if a slow stage needs to be computed
        if (q_slow < stage_ && a_fastslow_local_[j*stage_ + q_slow] != 0)
        {
          slow_time = t + h*c_slow_[q_slow];
          y_exp_s_->Update(1.0, *y, 0.0);

          //Accumulate slow stages
          for (int k = 0; k < q_slow; k++)
          {
            int index = q_slow * stage_ + k;
            if (a_slowslow_[index] != 0)
            {
              y_exp_s_->Update(a_slowslow_[index], *k_slow_[k], 1.0);
            }
          }

          //Accumulate fast coupled stages
          for (int k = 0; k < j; k++)
          {
            int index = j * stage_ + k;
            if (a_slowfast_local_[index] != 0)
            {
              y_exp_s_->Update(a_slowfast_local_[index], *k_fast_[k], 1.0);
            }
          }
          y_exp_s_->Update(1.0, *ya_slowfast_, 1.0);

          //Solve nonlinear system
          //FIXME: Linkage Errors 
          int flag = SolveNonlinearSystem(slow_time, slow_time_old, h, a_slowslow_[q_slow * (stage_ + 1)], y_ns_old, y_ns);

          if (flag)
          {
            
            //TODO: add flag based operations
          }
          

          fn_->ModifySolutionSlow(sum_timef, y_ns);
          fn_->FunctionalTimeDerivativeSlow(sum_timef, y_ns, k_slow_[j]);
          k_slow_[j]->Scale(h);

          q_slow++;
          *y_ns_old = *y_ns;
        }

        
        y_exp_f_->Update(1.0, *y_lambda_, 0.0);

        // Add coupled slow stages
        for (int k = 0; k < q_slow; ++k)
        {
          int index = q_slow * stage_ + k;
          if (a_fastslow_local_[index] != 0)
          {
            y_exp_f_->Update(a_fastslow_local_[index], *k_slow_[k], 1.0);
          }
        }

        // Add fast stages
        for (int k = 0; k < j; ++k)
        {
          int index = j * stage_ + k;
          if (a_fastfast_local_[index] != 0)
          {
            y_exp_f_->Update(a_fastfast_local_[ index ], *k_fast_[k], 1.0);
          }
        }

        fn_->ModifySolutionFast(sum_timef, y_exp_f_);
        fn_->FunctionalTimeDerivativeFast(sum_timef, y_exp_f_, k_fast_[j]);
        k_fast_[j]->Scale(h_fast);
      }

      for (int j = 0; j < stage_; j++)
      {
        if (b_fast_local_[j] != 0)
        {
          y_lambda_->Update(b_fast_local_[j], *k_fast_[j], 1.0);
        }
        for (int k = 0; k < stage_; k++)
        {
          int index = j*stage_ + k;
          if (a_slowfast_local_[index] != 0)
          {
            ya_slowfast_->Update(a_slowfast_local_[index], *k_slow_[k], 1.0);
          }
        }
      }
    }

    *y_new = *y_lambda_;
    for (int i = 0; i < stage_; i++)
    {
      if (b_slow_[i] != 0)
      {
        y_new->Update(b_slow_[i], *k_slow_[i], 1.0);
      }
    }

    //Time step was succeful
    return 0;
  }

  /**
   * @brief Solve the nonlinear system F(u) = 0 generated from implicit stages
   * 
   * TODO: Switch return type to an enum for updating results of failure.
   * 
   * FIXME: Add the verbose debugger displays
   * 
   * @tparam Vector vector type
   * @tparam VectorSpace map to vector
   * @param t_n1 
   * @param t_n0 
   * @param h 
   * @param scaling 
   * @param y_old 
   * @return int the a_slowslowociated error code TODO: Add error codes (Could be an enum instead)
   */
  template <class Vector, class VectorSpace>
  int MRG_EXIM_TI<Vector, VectorSpace>::SolveNonlinearSystem(double t_n1, double t_n0, double h, double scaling, Teuchos::RCP<Vector> y_old, const Teuchos::RCP<Vector>& y_new){
    solver_fn_slow_->SetTimes(t_n0, t_n1);
    solver_fn_slow_->SetPreviousTimeSolution(y_old);
    solver_fn_slow_->SetExplicitTerms(y_exp_s_);

     // Solve the nonlinear system.
    int ierr, code, itr;
    try {
      ierr = solverslow_->Solve(y_new);
      itr = solverslow_->num_itrs();
      code = solverslow_->returned_code();
    } catch (const Errors::CutTimeStep& e) {
      ierr = 1;
      itr = -1;  // This should not be summed up into the global counter.
      code = AmanziSolvers::SOLVER_INTERNAL_EXCEPTION;
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << e.what() << std::endl;
      }
    }


    return ierr;
  }

}

#endif