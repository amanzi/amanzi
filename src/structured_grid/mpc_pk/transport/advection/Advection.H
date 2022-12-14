/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _Advection_H_
#define _Advection_H_

#include <FArrayBox.H>
#include <Geometry.H>

class Advection
{
public:
  static int nGrowHyp();
  static int nGrowForce();
  static int nWork();
  static void FluxDivergence(const FArrayBox& CoCC, FArrayBox& CnCC,       int Ccomp,    int Cng,
                             const FArrayBox& SoCC, const FArrayBox& SnCC, int Scomp,    int Sng,
                             const PArray<FArrayBox>& UEC,                 int Ucomp,    int Ung,
                             FArrayBox& FDivCC,                            int FDivcomp, int FDivng,
                             PArray<FArrayBox>& FluxEC,                    int Flcomp,   int Flng,
                             const FArrayBox& FcCC,                        int Fccomp,   int Fcng,
                             const FArrayBox& FsCC,                        int FScomp,   int Fsng,
                             const FArrayBox& VolCC, const PArray<FArrayBox>& AEC, const FArrayBox& PhiCC,
                             const Box& vbox, const Real* dx, Real dt, const int* bc, int nc,
                             bool is_conservative);

  static void AdvUpdate(const FArrayBox& CoCC, FArrayBox& CnCC,       int Ccomp,    int Cng,
                        const FArrayBox& SoCC, const FArrayBox& SnCC, int Scomp,    int Sng,
                        FArrayBox& FDivCC,                            int FDivcomp, int FDivng,
                        PArray<FArrayBox>& FluxEC,                    int Flcomp,   int Flng,
                        const FArrayBox& PhiCC,                                     int Phing,
                        const Box& vbox, Real dt, int nc);

  static void EstimateMaxEigenvalues(const FArrayBox&        SatCC, int Scomp, int Satng,
                                     D_DECL(const FArrayBox& Ux,
                                            const FArrayBox& Uy,
                                            const FArrayBox& Uz),     int Ucomp, int Ung,
                                     const FArrayBox&        PhiCC,            int Phing,
                                     const Box&              vbox,
                                     Real*                   eigmax);

};

#endif
