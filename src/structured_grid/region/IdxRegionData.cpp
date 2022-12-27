/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <IdxRegionData.H>

IdxRegionData::IdxRegionData(const std::string&          label,
                             const Array<const Region*>& regions,
                             const std::string&          typeStr,
                             const IdxRDEval&            eval)
  : mLabel(label), mType(typeStr)
{
  mEvaluator = eval.clone();
  SetRegions(regions);
}

IdxRegionData::IdxRegionData(const std::string&          label,
                             const Array<const Region*>& regions,
                             const std::string&          typeStr,
                             Real                        val)
  : mLabel(label), mType(typeStr)
{
  mEvaluator = new IdxRDEval(val);
  SetRegions(regions);
}

IdxRegionData::IdxRegionData(const std::string&          label,
			     const Array<const Region*>& regions,
			     const std::string&          typeStr,
			     const Array<Real>&          vals,
			     const Array<Real>&          times,
			     const Array<std::string>&   forms)
  : mLabel(label), mType(typeStr)
{
  mEvaluator = new IdxRDEval(vals,times,forms);
  SetRegions(regions);
}

IdxRegionData::~IdxRegionData()
{
  delete mEvaluator;
}

void
IdxRegionData::SetRegions(const Array<const Region*>& regions)
{
  int nregions=regions.size();
  mRegions.resize(nregions);
  for (int i=0; i<nregions; ++i) {
    mRegions[i] = regions[i];
  }
}

void
IdxRegionData::apply(FArrayBox&       fab,
                     const IArrayBox& idx,
                     const Real*      dx,
                     int              dcomp,
                     Real             time) const
{
  const Box& box = fab.box();
  FArrayBox mask(box,1); mask.setVal(-1);
  for (int j=0; j<mRegions.size(); ++j) {
    mRegions[j]->setVal(mask,1,0,dx,0);
  }

  for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
    if (mask(iv,0) > 0) {
      const std::vector<Real>& val = (*mEvaluator)(idx(iv,0),time);
      for (int i=0; i<val.size(); ++i) {
        fab(iv,dcomp+i) = val[i];
      }
    }
  }
}

Array<Real>
IdxRegionData::IdxRDEval::Time() const
{
  return mFunc.x_;
};

Array<Real>
IdxRegionData::Time() const
{
  return mEvaluator->Time();
};

IdxRegionData::IdxRDEval::~IdxRDEval() {}

IdxRegionData::IdxRDEval *
IdxRegionData::IdxRDEval::clone () const
{
  return new IdxRDEval(*this);
}

IdxRegionData::IdxRDEval::IdxRDEval() {}

IdxRegionData::IdxRDEval::IdxRDEval(const IdxRDEval& rhs)
{
  mFunc = rhs.mFunc;
}

IdxRegionData::IdxRDEval::IdxRDEval(const Array<Real>&        vals,
				    const Array<Real>&        times,
				    const Array<std::string>& forms)
  : mFunc(times,vals,forms)
{
}

IdxRegionData::IdxRDEval::IdxRDEval(Real val)
  : mFunc(Array<Real>(1,0),Array<Real>(1,val),Array<std::string>(0))
{
}

std::vector<Real>
IdxRegionData::IdxRDEval::operator()(int i, Real t) const
{
  mRetData.resize(NComp());
  for (int n=0; n<NComp(); ++n) {
    mRetData[n] = mFunc(t);
  }
  return mRetData; // Copied
}
