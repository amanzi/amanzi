#ifndef _driver_support_H
#define _driver_support_H

#ifndef WIN32
#include <unistd.h>
#endif

#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <ParmParse.H>

#include <TensorDiffusionBndry.H>
#include <TensorOp.H>
#include <LO_BCTYPES.H>

namespace Amanzi {
  namespace AmanziTransport {

    const int MCLinOp_grow = 1; // FIXME: Need this to be a class variable
    const int Geom_Grow = 1;
    enum SolveMode {PREDICTOR, CORRECTOR, ONEPASS};

    void
    loadBndryData (MCInterpBndryData&  bd,
                   MultiFab&           S_fine, int sComp_S_fine,
                   MultiFab*           S_crse, int sComp_S_crse,
                   const Array<BCRec>& bc,
                   const Geometry&     geom,
                   int                 ratio,
                   int                 nComp);

    void
    diffuse_tracer(Real                   t_old,
                   Real                   t_new,
                   Real                   be_cn_theta,
                   const MultiFab&        S_old,     int sComp_S_old,
                   MultiFab&              S_new,     int sComp_S_new,
                   MultiFab*              W_old,     int sComp_W_old,
                   MultiFab*              W_new,     int sComp_W_new,
                   MultiFab*              W_half,    int sComp_W_half,
                   int                    W_flag,
                   const MCInterpBndryData& bd_old,  int sComp_bd_old,
                   const MCInterpBndryData& bd_new,  int sComp_bd_new,
                   MultiFab* const*       fluxn,
                   MultiFab* const*       fluxnp1,   int dComp_flux,
                   MultiFab*              delta_rhs, int sComp_rhs,
                   const MultiFab*        alpha_in,  int sComp_alpha_in,
                   const MultiFab* const* betan,     int sComp_betan,
                   const MultiFab* const* betanp1,   int sComp_betanp1,
                   const MultiFab* const* beta1n,    int sComp_beta1n,
                   const MultiFab* const* beta1np1,  int sComp_beta1np1,
                   int                    nComp,
                   const SolveMode&       solve_mode,
                   int                    max_order,
                   bool                   add_old_time_divFlux);

  } /* AmanziTransport */
} /* Amanzi */

#endif
