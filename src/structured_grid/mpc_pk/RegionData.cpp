#include <RegionData.H>
#include <PorousMedia.H>
#include <POROUSMEDIA_F.H>


RegionData::RegionData(const PArray<Region>& regions,
                       const std::string&    typeStr,
                       const Array<Real>&    vals)
    : type(typeStr), vals(vals), nComp(vals.size())
{
    setRegions(regions);
}

RegionData::RegionData(const PArray<Region>& regions,
                       const std::string&    typeStr,
                       Real                  val)
    : type(typeStr), nComp(1), vals(Array<Real>(val,nComp))
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

ArrayRegionData::ArrayRegionData(const Array<Array<Real> >&        x,
                                 const Array<Array<Real> >&        y,
                                 const Array<Array<std::string> >& form,
                                 const PArray<Region>&             regions,
                                 const std::string&                typeStr)
    : RegionData(regions,typeStr,Array<Real>(x.size()))
{
    BL_ASSERT(nComp==y.size());
    BL_ASSERT(nComp==form.size());
    funcs.resize(nComp);
    for (int i=0; i<nComp; ++i) {
        funcs[i] = TabularFunction(x[i],y[i],form[i]);
    }
}

ArrayRegionData::ArrayRegionData(const Array<Real>&       x,
                                 const Array<Real>&       y,
                                 const Array<std::string> form,
                                 const PArray<Region>&    regions,
                                 const std::string&       typeStr,
                                 int                      nComp)
    : RegionData(regions,typeStr,Array<Real>(nComp))
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
FluxToArrayBC::operator() (Real time) const
{
    // Here, a flux is converted to an array of rho.saturation values
    // based on the material properties of the rock at the boundary
    
    Real gravity = PorousMedia::getGravity();
    const Array<Real>& density = PorousMedia::Density(); // Assumes 1 component per phase
    int ncomps = density.size();
    BL_ASSERT(ncomps>0 && ncomps<=2);
    
    
    Real lkappa = rock.permeability[0];
    Real gstar;
    if (density.size() > 1)
        gstar = -lkappa*(density[0]-density[1])*gravity;
    else
        gstar = -lkappa*(density[0])*gravity;
    
    // Compute saturation given Aqueous flow rate
    int lkrtype = rock.krType;
    Real lkrcoef = rock.krParam[0];
    Real lsatres = rock.krParam[1];            
    int nc = 1;
    Real vtot = 0.; // Zero total velocity
    Real sol;
    const Array<Real>& visc = PorousMedia::Viscosity();
    
    Array<Real> rhoSat(ncomps);
    Real vel = func(time);
    FORT_FIND_INV_FLUX(&sol, &vel, &nc, &vtot,&gstar,visc.dataPtr(),&ncomps,&lkrtype,&lkrcoef);
    
    rhoSat[0] = density[0]*(sol*(1.0-lsatres)+lsatres);
    if (ncomps > 1) {
        rhoSat[1] = density[1]*(1.0-rhoSat[0]/density[0]);
    }
    
    return rhoSat;
}
