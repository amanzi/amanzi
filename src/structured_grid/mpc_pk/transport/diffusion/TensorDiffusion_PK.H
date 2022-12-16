/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
#include <ViscBndry.H>
#include <ABecHelper.H>
#include <TensorOp.H>
#include <LO_BCTYPES.H>

namespace Amanzi {
  namespace AmanziTransport {

    const int MCLinOp_grow = 1; // FIXME: Need this to be a class variable
    const int Geom_Grow = 1;

    TensorOp*
    getOp(Real                     a,
          Real                     b,
          const MCInterpBndryData& bd,       int sComp_bd, int nComp_bd,
          MultiFab*                W,        int sComp_W,  int nComp_W,
          MultiFab*                W_half,   int sComp_W_half,
          int                      W_flag,
          const MultiFab* const    beta[BL_SPACEDIM],     int sComp_beta, int nComp_beta,
          const MultiFab* const    beta1[BL_SPACEDIM],    int sComp_beta1, int nComp_beta1,
          const MultiFab&          volume,
          MultiFab*   area,
          const MultiFab*          alpha_in, int sComp_alpha_in);

    ABecHelper*
    getOp(Real                     a,
          Real                     b,
          const ViscBndry&         bd,       int sComp_bd,
          MultiFab*                W,        int sComp_W, int nComp_W,
          MultiFab*                W_half,   int sComp_W_half,
          int                      W_flag,
          const MultiFab* const  beta[BL_SPACEDIM],     int sComp_beta, int nComp_beta,
          const MultiFab&          volume,
          MultiFab*  area,
          const MultiFab*          alpha_in, int sComp_alpha_in,
          int                      nComp);

  } /* AmanziTransport */
} /* Amanzi */

#endif
