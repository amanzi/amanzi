/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _DiffDomRelSrc_H
#define _DiffDomRelSrc_H

#include <RegionData.H>

/*
  Diffusion-Dominated Release model for AmanziS

  This is a special form of a RegionData object, used to evaluate a source term of the form:

  Q(time) = I/L . Sqrt(Deff / (Pi . time))

  Where:
     I     total inventory
     L     mixing length of source
     Deff  effective diffusion coefficient

  such that, at time = t, a total of M(t) = Int(Q . dt)_0^t = 2 (I/L) Sqrt(t . Deff / Pi)
  will enter the domain, per unit volume.

  There are (at least) two potential problems that arise when a volumetric source term takes
  this (nonlinear) form:

  (1) If a first-order time integration is used, a pointwise evaluation of Q will lead to an
  integral of the source term != M, where the discrepancy is a function of the time stepsize.
  Thus, if used with temporal subcycling, we must formulate to guarantee consistent sources.

  (2) The form is not defined as time -> 0, and thus will not interact properly with a
  convergence analysis.

  Here, we instead return Qav, the average value of Q over the interval (t1:t2):

  Qav(t1,t2) = Qav2 = ( M(t2) - M(t1) ) / (t2 - t1), provided t2 > t1 > tau

  Here, tau is an arbitrary short timescale, such that if t2 < tau,

  Qav = Qav1 = ( M(tau) - M(0) ) / tau

  The final result returned is a weighted average of Qav1 and Qav2, based on t1 and t2
  relative to tau:

        Qav(t1,t2) = eta.Qav1 + (1 - eta).Qav2( max(tau, t1), t2)
        eta = max(0, (tau - t1)/(t2 - t1)

 */
class DiffDomRelSrc
  : public RegionData
{
public:
  DiffDomRelSrc(const std::string&          label,
                const Array<const Region*>& regions,
                const std::string&          typeStr,
                Real                        mixingLength,
                Real                        Deff,
                Real                        totalInventory,
                Real                        startTime,
                Real                        endTime,
                Real                        timeScale);

  virtual ~DiffDomRelSrc() {}
  virtual Array<Real> operator() (Real time = 0) const;
  virtual Array<Real> operator() (Real t1, Real t2) const;
  virtual void apply(FArrayBox&  fab,
                     const Real* dx,
                     int         scomp,
                     int         ncomp,
                     Real        time = 0) const;
  virtual void apply(FArrayBox&  fab,
                     const Real* dx,
                     int         scomp,
                     int         ncomp,
                     Real        t1,
                     Real        t2) const;
protected:
  Real Func(Real t1, Real t2) const;

  Real mMixingLength, mDeff, mTotalInventory, mStartTime, mEndTime, mTimeScale, mQav1;
};

#endif
