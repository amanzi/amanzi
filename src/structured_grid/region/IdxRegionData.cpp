#include <IdxRegionData.H>

IdxRegionData::IdxRegionData(const std::string&          label,
                             const Array<const Region*>& regions,
                             const std::string&          typeStr,
                             IdxRDEval*                  eval)
  : mLabel(label), mType(typeStr), mEvaluator(eval)
{
  BL_ASSERT(mEvaluator != 0);
  SetRegions(regions);
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
  BL_ASSERT(mEvaluator != 0);

  const Box& box = fab.box();
  FArrayBox mask(box,1); mask.setVal(-1);
  for (int j=0; j<mRegions.size(); ++j) { 
    mRegions[j]->setVal(mask,1,0,dx,0);
  }

  for (IntVect iv=box.smallEnd(); iv<=box.bigEnd(); box.next(iv)) {
    if (mask(iv,0) > 0) {
      Real val = (*mEvaluator)(idx(iv,0),time);
      fab(iv,dcomp) = val;
    }
  }
}

