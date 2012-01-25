#ifndef BCDATA_H_
#define BCDATA_H_

#include <TabularFunction.H>
#include <Rock.H>

class RegionData
{
public:
    RegionData(const PArray<Region>& regions,
               const std::string&    typeStr,
               const Array<Real>&    vals);

    RegionData(const PArray<Region>& regions,
               const std::string&    typeStr,
               Real                  val);

    const PArray<Region>& Regions() const {return regions;}
    std::string Type() const {return type;}
    virtual Array<Real> operator() (Real time=0) const {return vals;}

protected:
    void setRegions(const PArray<Region>& regions_);

    int nComp;
    std::string type;
    PArray<Region> regions;
    Array<Real> vals;
};

class ArrayRegionData
    : public RegionData
{
public:
    ArrayRegionData(const Array<Array<Real> >&        x,
                    const Array<Array<Real> >&        y,
                    const Array<Array<std::string> >& form,
                    const PArray<Region>&             regions,
                    const std::string&                typeStr);

    ArrayRegionData(const Array<Real>&       x,
                    const Array<Real>&       y,
                    const Array<std::string> form,
                    const PArray<Region>&    regions,
                    const std::string&       typeStr,
                    int                      nComp);
    
    virtual Array<Real> operator() (Real time=0) const;

protected:
    Array<TabularFunction> funcs;
};



class FluxToArrayBC
    : public RegionData
{
public:
    FluxToArrayBC(const Array<Real>&        x,
                  const Array<Real>&        y,
                  const Array<std::string>& form,
                  const PArray<Region>&     regions,
                  const std::string&        typeStr,
                  int                       nComp,
                  const Rock&               rock)
        : RegionData(regions,typeStr,Array<Real>(nComp)), func(x,y,form), rock(rock) {}
    
    virtual Array<Real> operator() (Real time) const;

protected:
    TabularFunction func;
    const Rock& rock;
};




#endif
