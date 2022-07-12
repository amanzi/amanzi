
/*
This class implements Implicit-Implicit Multirate Time Integration Methods:

TODO: expand documentation
<List of methods>

 - User Defined


*/

#ifndef AMANZI_MRG_IMIM_TI_HH_
#define AMANZI_MRG_IMIM_TI_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include <functional>

#include "PartitionFnBase.hh"
#include "MRG_IMIM_FnBase.hh"
#include "MRG_IMIM_SolverFnBase.hh"
#include "VerboseObject.hh"
#include "../solvers/Solver.hh"
#include "../solvers/SolverFactory.hh"

namespace Amanzi
{

  /**
   * @brief Method enums
   * TODO: Add Names and an order 3 method
   */
  enum method_t
  {
    temp2,
    temp3,
    user_defiend
  };

  template <class Vector, class VectorSpace>
  class MRG_IMIM_TI
  {
  private:
    Teuchos::RCP<MRG_IMIM_FnBase<Vector>> fn_;

    Teuchos::RCP<AmanziSolvers::Solver<Vector,VectorSpace> > solverfull_;
    Teuchos::RCP<BDF1_SolverFnBase<Vector> > solver_fnfull_;

    Teuchos::RCP<AmanziSolvers::Solver<Vector,VectorSpace> > solverfast_;
    Teuchos::RCP<BDF1_SolverFnBase<Vector> > solver_fnfast_;

    int order_;
    int stage_;
    std::vector<double> a_fastfast_, b_fast_, c_fast_;
    std::vector<double> a_slowslow_, b_slow_, c_slow_;
    // Both functions for coupling will take M and lambda respectively
    std::function<void(double, double, std::vector<double>&)> a_slowfast_;
    std::function<void(double, double, std::vector<double>&)> a_fastslow_;

    // Allocate the coupling matrices based of M and update when a new M is given
    // TODO:
    int prev_M_ = 0;
    std::vector<double> a_fastfast_local_, b_fast_local_, c_fast_local_;
    std::vector<double> a_slowfast_local_;
    std::vector<double> a_fastslow_local_;

    std::vector<Teuchos::RCP<Vector>> k_fast_;
    std::vector<Teuchos::RCP<Vector>> k_slow_;
    Teuchos::RCP<Vector> ya_slowfast_;
    Teuchos::RCP<Vector> y_lambda_;
    Teuchos::RCP<Vector> y_exp_f_;
    Teuchos::RCP<Vector> y_exp_s_;
    Teuchos::RCP<Vector> y_ns;
    Teuchos::RCP<Vector> y_ns_old;
    Teuchos::RCP<Vector> y_nf;
    Teuchos::RCP<Vector> y_nf_old;

    void InitMethod(const method_t method);
    void InitMemory(const Teuchos::RCP<Vector> initvector);
    void InitCoefficentMemory();

    //TODO: Set up Construction of both Nonlinear Solvers
    void InitSolvers();

  public:
    MRG_IMIM_TI(MRG_IMIM_FnBase<Vector> &fn,
               const method_t method,
               Teuchos::ParameterList& plist,
               const Teuchos::RCP<Vector> initvector);

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
               Teuchos::ParameterList& plist,
               const Teuchos::RCP<Vector> initvector);

    void TimeStep(const double t, const double h, const int M, const Teuchos::RCP<Vector> y, Teuchos::RCP<Vector> y_new);

    int order() { return order_; };
    int stage() { return stage_; };
  };

  template <class Vector, class VectorSpace>
  void MRG_IMIM_TI<Vector, VectorSpace>::InitMethod(const method_t method)
  {

    switch (method)
    {
    case temp2:
      // Order 2 method
      {
      double lam = 1 - 1/sqrt(2);

      a_fastfast_ = {lam, 0, 1 - lam, lam};
      b_fast_= {1 - lam, lam};
      c_fast_ = {lam, 1};

      a_slowslow_ = {lam, 0, 1 - lam, lam};
      b_slow_= {1 - lam, lam};
      c_slow_ = {lam, 1};

      // The conditional coupling coefficents
      a_fastslow_ = [lam] ( double M, double lambda, std::vector<double>& a){
        a[0] = (lambda - 1 + lam) / M;
        if (lambda == M)
        {
          a[2] = 1 - lam;
          a[3] = lam;
        }
        else
        {
          a[2] = lambda / M;
          a[3] = 0;
        }
      };

      a_slowfast_ = [lam] ( double M, double lambda, std::vector<double>& a){
        a[0] = (M == 1) ? lam : 0;
        a[2] = (1 - lam) / M;
        a[3] = lam / M;
      };
      }

      
      stage_ = 2;
      order_ = 2;


      break;

    case temp3:
      // Order 3 method
      stage_ = 3;
      order_ = 3;

      // TODO: get coefficents
      break;
    default:
      // TODO: warn user of improper method
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
  void MRG_IMIM_TI<Vector, VectorSpace>::InitCoefficentMemory()
  {
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
  void MRG_IMIM_TI<Vector, VectorSpace>::InitMemory(const Teuchos::RCP<Vector> initvector)
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
    y_nf = Teuchos::rcp(new Vector(*initvector));
    y_nf_old = Teuchos::rcp(new Vector(*initvector));
  }

  template <class Vector, class VectorSpace>
  MRG_IMIM_TI<Vector, VectorSpace>::MRG_IMIM_TI(MRG_IMIM_FnBase<Vector> &fn,
                                              const method_t method,
                                              Teuchos::ParameterList& plist,
                                              const Teuchos::RCP<Vector> initvector)
  {
    
    fn_ = Teuchos::rcpFromRef(fn);

    MRG_IMIM_TI<Vector>::InitMethod(method);
    MRG_IMIM_TI<Vector>::InitMemory(initvector);
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
                                              const std::function<void(double, double, std::vector<double>&)> a_slowfast_,
                                              const std::function<void(double, double, std::vector<double>&)> a_fastslow_,
                                              Teuchos::ParameterList& plist,
                                              const Teuchos::RCP<Vector> initvector) : order_(order), stage_(stage), method_(user_defined),
                                                                          a_slowfast_(a_slowfast), a_fastslow_(a_fastslow)
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

    MRG_IMIM_TI<Vector>::InitMemory(initvector);
  }


  template <class Vector, class VectorSpace>
  void MRG_IMIM_TI<Vector, VectorSpace>::TimeStep(const double t, const double h, const int M, const Teuchos::RCP<Vector> y, Teuchos::RCP<Vector> y_new)
  {

    // Could simplify with templates and commands as many operations repeat

    double sum_timef = t;
    double M_cast = static_cast<double>(M);

    *y_lambda_ = *y;
    int q_slow = 0;
    double h_fast = h / static_cast<double>(M);

    // Intalize coefficents for scaled M
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

    for (int i = 0; i < M; ++i)
    {
      // intialize the coupling coefficents
      a_fastslow_(M_cast, static_cast<double>(i + 1), a_fastslow_local_);
      a_slowfast_(M_cast, static_cast<double>(i + 1), a_slowfast_local_ );

      for (int j = 0; j < stage_; ++j)
      {
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
            y_exp_f_->Update(a_fastfast_local_[index], *k_fast_[k], 1.0);
          }
        }

        if (q_slow <= stage_ && a_fastslow_local_[j*stage_ + q_slow] != 0 && a_slowfast_local_[q_slow * stage + j])
        {
          y_exp_s->Update(1.0, *y, 0.0);

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
            int index = j * stage_ + k;
            if (a_slowfast_local_[index] != 0)
            {
              y_exp_s_->Update(a_slowfast_local_[index], *k_fast_[k], 1.0);
            }
          }
          y_exp_s_->Update(1.0, *ya_slowfast_, 1.0);

          /*
            TODO:
            Will solve nonlinear system for the coupled problem

            Will return yf and ys as solutions
          */

          fns_->ModifySolutionSlow(sum_timef, y_ns);
          fns_->FunctionalTimeDerivativeSlow(sum_timef, y_ns, *k_slow_[j]);
          k_slow_[j]->Scale(h);

          q_slow++;
        }
        else
        {
          /*
            TODO:

            Will solve nonlinear system for only the fast system

            Will return yf
          */
        }

        fnf_->ModifySolutionFast(sum_timef, y_nf);
        fnf_->FunctionalTimeDerivativeFast(sum_timef, y_nf, *k_fast_[j]);
        k_fast_[j]->Scale(h_fast);
      }

      for (int j = 0; j < stage_; j++)
      {
        if (b_fast_local_[j] != 0)
        {
          y_lambda_.Update(b_fast_local_[j], *k_fast_[j], 1.0);
        }
        for (int k = 0; k < stage_; k++)
        {
          int index = j*stage_ + k;
          if (a_slowfast_local_[index] != 0)
          {
            ka_slowfast_.Update(a_slowfast_local_[index], *k_slow_[k], 1.0);
          }
        }
      }
    }

    y_new = y_lambda_;
    for (int i = 0; i < stage_; i++)
    {
      if (b_slow_[i] != 0)
      {
        y_new->Update(b_slow_[i], *k_slow_[i], 1.0);
      }
    }
  }

}

#endif