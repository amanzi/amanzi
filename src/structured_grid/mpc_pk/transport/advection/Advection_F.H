#ifndef _Advection_F_H_
#define _Advection_F_H_

#include <REAL.H>

#if defined(BL_FORT_USE_UPPERCASE)
# define BDS_EDGE_STATES    BDS_EDGE_STATES
# define BDS_SET_NGROWHYP   BDS_SET_NGROWHYP
# define BDS_SET_NGROWFORCE BDS_SET_NGROWFORCE
# define BDS_SET_NWORK      BDS_SET_NWORK
# define ADV_FDIV           ADV_FDIV
# define ADV_UPDATE         ADV_UPDATE
# define BDS_EST_EIGEN      BDS_EST_EIGEN
#elif defined(BL_FORT_USE_LOWERCASE)
# define BDS_EDGE_STATES    bds_edge_states
# define BDS_SET_NGROWHYP   bds_set_ngrowhyp
# define BDS_SET_NGROWFORCE bds_set_ngrowforce
# define BDS_SET_NWORK      bds_set_nwork
# define ADV_FDIV           adv_fdiv
# define ADV_UPDATE         adv_update
# define BDS_EST_EIGEN      bds_est_eigen
#elif defined(BL_FORT_USE_UNDERSCORE)
# define BDS_EDGE_STATES    bds_edge_states_
# define BDS_SET_NGROWHYP   bds_set_ngrowhyp_
# define BDS_SET_NGROWFORCE bds_set_ngrowforce_
# define BDS_SET_NWORK      bds_set_nwork_
# define ADV_FDIV           adv_fdiv_
# define ADV_UPDATE         adv_update_
# define BDS_EST_EIGEN      bds_est_eigen_
#endif

#include <ArrayLim.H>

extern "C" {
  void BDS_SET_NGROWHYP(int* nGrowHyp);
  void BDS_SET_NGROWFORCE(int* nGrowForce);
  void BDS_SET_NWORK(int* nWork);
  void BDS_EDGE_STATES(const Real* s, const int* ng_s, const Real* ci, const int* ng_ci,
                       Real* sedgex, Real* sedgey, const int* ng_se,
                       const Real* umac, const Real* vmac, const int* ng_u, const Real* force, const int* ng_f,
                       Real* slope, const int* ng_c, const int* nw, const Real* dx, const Real* dt,
                       const int* num_comp, const int* is_cons,
                       const int* lo, const int* hi, const int* bc);

  void ADV_FDIV(Real* aofs, const int* ng_a, Real* edgex, Real* edgey, const int* ng_se,
                const Real* u_dat, const Real* v_dat, const int* ng_u,
                const Real* ax_dat, const Real* ay_dat, const Real* vol_dat,
                const int* num_comp, const int* lo, const int* hi);

  void ADV_UPDATE(const Real* sold, Real* snew, const int* ng_s,
                  const Real* capo, const Real* capn, const int* ng_cap,
                  const Real* aofs, const int* ng_a, const Real* phi, const int* ng_phi,
                  const Real* dt, const int* num_comp, const int* lo, const int* hi);

  void BDS_EST_EIGEN(const Real* s, const int* ng_s, const Real* umac, const Real* vmac, const int* ng_u,
                     const Real* phi, const int * ng_phi, const int* lo, const int* hi, Real* eigenval);

}
#endif
