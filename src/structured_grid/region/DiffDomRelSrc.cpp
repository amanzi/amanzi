/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <DiffDomRelSrc.H>

static Real Pi = 2 * std::asin(1.0);

DiffDomRelSrc::DiffDomRelSrc(const std::string&          label,
                             const Array<const Region*>& regions,
                             const std::string&          typeStr,
                             Real                        mixingLength,
                             Real                        Deff,
                             Real                        totalInventory,
                             Real                        startTime,
                             Real                        endTime,
                             Real                        timeScale)
  : RegionData(label, regions, typeStr, 0),
    mMixingLength(mixingLength), mDeff(Deff),
    mTotalInventory(totalInventory), mStartTime(startTime),
    mEndTime(endTime), mTimeScale(timeScale)
{
  BL_ASSERT(mMixingLength > 0);
  BL_ASSERT(mDeff > 0);
  BL_ASSERT(mTimeScale > 0);
  mQav1 = (2 * mTotalInventory / mMixingLength) * std::sqrt(Deff / (mTimeScale * Pi));
}

void
DiffDomRelSrc::apply(FArrayBox&  fab,
                     const Real* dx,
                     int         scomp,
                     int         ncomp,
                     Real        time) const
{
  BoxLib::Abort("Must evaluate DiffDomRelSrc over time interval");
};

Array<Real>
DiffDomRelSrc::operator() (Real time) const
{
  BoxLib::Abort("Must evaluate DiffDomRelSrc over time interval");
};

Real
DiffDomRelSrc::Func(Real t1, Real t2) const
{
  BL_ASSERT(t2 > t1);
  Real M1 = 2 * (mTotalInventory / mMixingLength) * std::sqrt(mDeff * t1 / Pi);
  Real M2 = 2 * (mTotalInventory / mMixingLength) * std::sqrt(mDeff * t2 / Pi);
  return (M2 - M1) / (t2 - t1);
}

Array<Real>
DiffDomRelSrc::operator() (Real t1, Real t2) const
{
  BL_ASSERT(t2 > t1);
  Real t1s = t1 - mStartTime;
  Real t2s = t2 - mStartTime;
  Real eta = std::max(0., (std::min(mTimeScale,t2s) - t1s) / (t2 - t1)); // Fraction source is in intitial phase (t<mTimeScale)
  Real eta2 = std::max(0., t2 - std::max(t1,mEndTime)) / (t2 - t1);      // Fraction source is off
  Real Qav2 = (eta == 1 || eta2 == 1  ?  0 : Func( std::max(mTimeScale, t1s), std::min(t2,mEndTime)));
  Real ret = eta * mQav1  +  (1. - eta - eta2) * Qav2;
  return Array<Real>(1,ret);
}

void
DiffDomRelSrc::apply(FArrayBox&  fab,
                     const Real* dx,
                     int         scomp,
                     int         ncomp,
                     Real        t1,
                     Real        t2) const
{
  Array<Real> val = (*this)(t1,t2);
  const Box& box = fab.box();
  FArrayBox mask(box,1); mask.setVal(-1);

  for (int j=0; j<regions.size(); ++j) {
    regions[j]->setVal(mask,1,0,dx,0);
  }

  for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
    if (mask(iv,0) > 0) {
      for (int n=0; n<ncomp; ++n) {
        fab(iv,scomp+n) = val[n];
      }
    }
  }
}
