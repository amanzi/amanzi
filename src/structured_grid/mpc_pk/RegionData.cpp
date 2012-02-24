#include <RegionData.H>
#include <PorousMedia.H>
#include <POROUSMEDIA_F.H>


RegionData::RegionData(const std::string&    label,
                       const PArray<Region>& regions,
                       const std::string&    typeStr,
                       const Array<Real>&    vals)
    : label(label), type(typeStr), vals(vals), nComp(vals.size())
{
    setRegions(regions);
}

RegionData::RegionData(const std::string&    label,
                       const PArray<Region>& regions,
                       const std::string&    typeStr,
                       Real                  val)
    : label(label), type(typeStr), nComp(1), vals(Array<Real>(nComp,val))
{
    setRegions(regions);
}


void
RegionData::setRegions(const PArray<Region>& regions_)
{
    regions.clear();
    int nregions=regions_.size();

    // Get a copy of the pointers to regions in a structure that wont 
    //   remove them when it leaves scope
    regions.resize(nregions,PArrayNoManage);
    for (int i=0; i<nregions; ++i)
    {
        Region& r = const_cast<Region&>(regions_[i]);
        regions.set(i,&(r));
    }
}

void
RegionData::apply(FArrayBox&  fab, 
                  const Real* dx, 
                  int         scomp,
                  int         ncomp,
                  Real        time) const
{
    Array<Real> val = (*this)(time);
    for (int j=0; j<regions.size(); ++j)
    { 
        regions[j].setVal(fab,val,dx,0,0,val.size());
    }
}

ArrayRegionData::ArrayRegionData(const std::string&                label,
                                 const Array<Array<Real> >&        x,
                                 const Array<Array<Real> >&        y,
                                 const Array<Array<std::string> >& form,
                                 const PArray<Region>&             regions,
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
                                 const PArray<Region>&    regions,
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

