/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <RegionData.H>

RegionData::RegionData(const std::string&    label,
                       const Array<const Region*>& regions,
                       const std::string&    typeStr,
                       const Array<Real>&    vals)
    : label(label), type(typeStr), vals(vals), nComp(vals.size())
{
    setRegions(regions);
}

RegionData::RegionData(const std::string&    label,
                       const Array<const Region*>& regions,
                       const std::string&    typeStr,
                       Real                  val)
    : label(label), type(typeStr), nComp(1), vals(Array<Real>(nComp,val))
{
    setRegions(regions);
}


void
RegionData::setRegions(const Array<const Region*>& regions_)
{
    int nregions=regions_.size();

    // Get a copy of the pointers to regions in a structure that wont
    //   remove them when it leaves scope
    regions.resize(nregions);
    for (int i=0; i<nregions; ++i) {
      regions[i] = regions_[i];
    }
}

void
RegionData::apply(FArrayBox&  fab,
                  const Real* dx,
                  int         scomp,
                  int         ncomp,
                  Real        t) const
{
    Array<Real> val = (*this)(t);
    const Box& box = fab.box();
    FArrayBox mask(box,1); mask.setVal(-1);
    for (int j=0; j<regions.size(); ++j)
    {
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

ArrayRegionData::ArrayRegionData(const std::string&                label,
                                 const Array<Array<Real> >&        x,
                                 const Array<Array<Real> >&        y,
                                 const Array<Array<std::string> >& form,
                                 const Array<const Region*>&             regions,
                                 const std::string&                typeStr)
    : RegionData(label,regions,typeStr,Array<Real>(x.size()))
{
    BL_ASSERT(nComp==y.size());
    BL_ASSERT(nComp==form.size());
    funcs.resize(nComp);
    for (int i=0; i<nComp; ++i) {
        funcs[i] = TabularFunction(x[i],y[i],form[i]);
    }
}

ArrayRegionData::ArrayRegionData(const std::string&       label,
                                 const Array<Real>&       x,
                                 const Array<Real>&       y,
                                 const Array<std::string> form,
                                 const Array<const Region*>&    regions,
                                 const std::string&       typeStr,
                                 int                      nComp)
    : RegionData(label,regions,typeStr,Array<Real>(nComp))
{
    funcs.resize(1);
    funcs[0] = TabularFunction(x,y,form);
}

Array<Real>
ArrayRegionData::operator() (Real time) const
{
    Array<Real> newVals(nComp);
    for (int i=0; i<nComp; ++i) {
        if (funcs.size()==1) {
            newVals[i] = funcs[0](time);
        }
        else {
            newVals[i] = funcs[i](time);
        }
    }
    return newVals;
};

Array<Real>
ArrayRegionData::time() const
{
  Array<Real> newVals;
  if (funcs.size()>0)
      newVals = funcs[0].x_;
  return newVals;
};
