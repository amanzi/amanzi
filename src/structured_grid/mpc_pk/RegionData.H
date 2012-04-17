#ifndef BCDATA_H_
#define BCDATA_H_

#include <TabularFunction.H>
#include <Rock.H>

class RegionData
{
public:
    RegionData(const std::string&    label,
               const PArray<Region>& regions,
               const std::string&    typeStr,
               const Array<Real>&    vals);

    RegionData(const std::string&    label,
               const PArray<Region>& regions,
               const std::string&    typeStr,
               Real                  val);

    const PArray<Region>& Regions() const {return regions;}
    std::string Type() const {return type;}
    std::string Label() const {return label;}
    virtual Array<Real> operator() (Real time=0) const {return vals;}
    virtual Array<Real> time() const {return vals;}

    virtual void apply(FArrayBox&  fab, 
                       const Real* dx, 
                       int         scomp,
                       int         ncomp,
                       Real        time = 0) const;
    

protected:
    void setRegions(const PArray<Region>& regions_);

    int nComp;
    std::string type, label;
    PArray<Region> regions;
    Array<Real> vals;
};

class ArrayRegionData
    : public RegionData
{
public:
    ArrayRegionData(const std::string&                label,
                    const Array<Array<Real> >&        x,
                    const Array<Array<Real> >&        y,
                    const Array<Array<std::string> >& form,
                    const PArray<Region>&             regions,
                    const std::string&                typeStr);

    ArrayRegionData(const std::string&       label,
                    const Array<Real>&       x,
                    const Array<Real>&       y,
                    const Array<std::string> form,
                    const PArray<Region>&    regions,
                    const std::string&       typeStr,
                    int                      nComp);
    
    virtual Array<Real> time() const;

    virtual Array<Real> operator() (Real time=0) const;

protected:
    Array<TabularFunction> funcs;
};


struct ArrayTransform
{
    virtual ArrayTransform* clone() const {return new ArrayTransform(*this);}
    virtual Array<Real> transform(Real time) const {return Array<Real>(1,time);}
};


class Transform_S_AR_For_BC
    : public RegionData
{
public:
    Transform_S_AR_For_BC(const std::string&        label,
                          const Array<Real>&        x,
                          const Array<Real>&        y,
                          const Array<std::string>& form,
                          const PArray<Region>&     regions,
                          const std::string&        typeStr,
                          int                       nComp,
                          ArrayTransform&           tfunc)
        : RegionData(label,regions,typeStr,Array<Real>(nComp)), func(x,y,form), tfunc(tfunc.clone()) {}

    ~Transform_S_AR_For_BC() {delete tfunc;}
    virtual Array<Real> operator() (Real time) const {
        Real t1 = func(time);
        const Array<Real> t2 = tfunc->transform(t1);
        return t2;
    }
    Real tabularFunction(Real time) const  {return func(time);}

protected:
    TabularFunction func;
    ArrayTransform* tfunc;
};




#endif
