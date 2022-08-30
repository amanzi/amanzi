#ifndef BCDATA_H_
#define BCDATA_H_

#include <PArray.H>
#include <Region.H>
#include <TabularFunction.H>

class RegionData
{
public:
    RegionData(const std::string&          label,
               const Array<const Region*>& regions,
               const std::string&          typeStr,
               const Array<Real>&          vals);

    RegionData(const std::string&          label,
               const Array<const Region*>& regions,
               const std::string&          typeStr,
               Real                        val);
    virtual ~RegionData() {};

    const Array<const Region*>& Regions() const {return regions;}
    std::string Type() const {return type;}
    std::string Label() const {return label;}
    virtual Array<Real> operator() (Real time=0) const {return vals;}
    virtual Array<Real> time() const {return Array<Real>(vals.size(),0);}

    virtual void apply(FArrayBox&  fab, 
                       const Real* dx, 
                       int         scomp,
                       int         ncomp,
                       Real        time = 0) const;
    

protected:
    void setRegions(const Array<const Region*>& regions_);

    int nComp;
    std::string type, label;
    Array<const Region*> regions;
    Array<Real> vals;
};

class ArrayRegionData : public RegionData {
public:
    ArrayRegionData(const std::string&                label,
                    const Array<Array<Real> >&        x,
                    const Array<Array<Real> >&        y,
                    const Array<Array<std::string> >& form,
                    const Array<const Region*>&       regions,
                    const std::string&                typeStr);

    ArrayRegionData(const std::string&          label,
                    const Array<Real>&          x,
                    const Array<Real>&          y,
                    const Array<std::string>    form,
                    const Array<const Region*>& regions,
                    const std::string&          typeStr,
                    int                         nComp);
    virtual ~ArrayRegionData() {}
    
    virtual Array<Real> time() const;

    virtual Array<Real> operator() (Real time=0) const;

protected:
    Array<TabularFunction> funcs;
};

struct ArrayTransform
{
  virtual ~ArrayTransform() {};
  virtual ArrayTransform* clone() const { return new ArrayTransform(*this); }
  virtual Array<Real> transform(Real time) const { return Array<Real>(1,time); }
};


class Transform_S_AR_For_BC : public RegionData {
 public:
  Transform_S_AR_For_BC(const std::string&          label,
                        const Array<Real>&          x,
                        const Array<Real>&          y,
                        const Array<std::string>&   form,
                        const Array<const Region*>& regions,
                        const std::string&          typeStr,
                        int                         nComp,
                        ArrayTransform&             tfunc)
      : RegionData(label,regions,typeStr,Array<Real>(nComp)), func(x,y,form), tfunc(tfunc.clone()) {}

  ~Transform_S_AR_For_BC() { delete tfunc; }
  virtual Array<Real> operator() (Real time) const override {
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
