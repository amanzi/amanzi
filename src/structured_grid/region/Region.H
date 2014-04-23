#ifndef _REGION_H_
#define _REGION_H_

#include <string>
#include <map>
#include <sstream>

#include <Array.H>
#include <FArrayBox.H>
#include <Box.H>
#include <BoxArray.H>
#include <ParallelDescriptor.H>
using std::map;

//enum PurposeT {any, rock, source, output, boundary}; 

//
// Base abstract class
//
class Region
{
public: 

  static std::map<std::string,int> create_region_map()
  {
    std::map<std::string,int> m;
    m["point"] = 0;
    m["box"]   = 1;
    m["color_function"] = 2;
    m["surface"]   = 3;

    return m;
  };

  Region(const std::string& r_name,
         const std::string& r_purpose,
         const std::string& r_type);
  virtual ~Region() {}

  // Note: No virtual templates, must manually create
  virtual void setVal(BaseFab<Real>&     fab, 
                      const Array<Real>& val, 
                      const Real*        dx, 
                      int                ng,
                      int                scomp,
                      int                ncomp) const;

  virtual void setVal(BaseFab<Real>& fab,
                      const Real&    val, 
                      const int      idx,
                      const Real*    dx, 
                      int            ng) const;
  
  virtual void setVal(BaseFab<int>&     fab, 
                      const Array<int>& val, 
                      const Real*       dx, 
                      int               ng,
                      int               scomp,
                      int               ncomp) const;

  virtual void setVal(BaseFab<int>& fab,
                      const int&    val, 
                      const int     idx,
                      const Real*   dx, 
                      int           ng) const;
  
  virtual bool inRegion (const Array<Real>& x) const = 0;
  virtual std::ostream& print (std::ostream& os) const;

  // Identifier for the region 
  std::string name;
  std::string purpose;
  std::string type;

  static Real geometry_eps;
  static Array<Real> domlo, domhi;

protected:

  template<class T>
  void setVal_DoIt(BaseFab<T>&     fab, 
                   const Array<T>& val, 
                   const Real*     dx, 
                   int             ng,
                   int             scomp,
                   int             ncomp) const;
};

std::ostream& operator<< (std::ostream& os, const Region& rhs);

class PointRegion : public Region
{
public: 
  
  PointRegion (const std::string& r_name,
               const std::string& r_purpose,
               const Array<Real>& coor)
      : Region(r_name,r_purpose,"point"), coor(coor) {} ;


  // Note: No virtual templates, must manually create
  virtual void setVal(BaseFab<Real>&     fab, 
                      const Array<Real>& val, 
                      const Real*        dx, 
                      int                ng,
                      int                scomp,
                      int                ncomp) const;

  virtual void setVal(BaseFab<Real>& fab,
                      const Real&    val, 
                      const int      idx,
                      const Real*    dx, 
                      int            ng) const;
  
  virtual void setVal(BaseFab<int>&     fab, 
                      const Array<int>& val, 
                      const Real*       dx, 
                      int               ng,
                      int               scomp,
                      int               ncomp) const;

  virtual void setVal(BaseFab<int>& fab,
                      const int&    val, 
                      const int     idx,
                      const Real*   dx, 
                      int           ng) const;
  
  virtual ~PointRegion() {}

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;

  Array<Real> coor;

protected:

  template<class T>
  void setVal_DoIt(BaseFab<T>&     fab, 
                   const Array<T>& val, 
                   const Real*     dx, 
                   int             ng,
                   int             scomp,
                   int             ncomp) const;
};

std::ostream& operator<< (std::ostream& os, const PointRegion& rhs);

class BoxRegion : public Region
{
public: 
  
  BoxRegion (const std::string& r_name,
             const std::string& r_purpose,
             const Array<Real>& lo,
             const Array<Real>& hi)
      : Region(r_name,r_purpose,"box"), lo(lo), hi(hi) {} ;
  virtual ~BoxRegion() {}

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;

  Array<Real> lo;
  Array<Real> hi;
};

std::ostream& operator<< (std::ostream& os, const BoxRegion& rhs);


#if BL_SPACEDIM == 2
class PolygonRegion : public Region
{
public: 
  enum PLANE {XYPLANE, YZPLANE, XZPLANE};

  PolygonRegion (const std::string&   r_name,
                 const std::string&   r_purpose,
                 const Array<Real>&   v1vector,
                 const Array<Real>&   v2vector
    );
  virtual ~PolygonRegion() {}

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;

  Array<Real> vector1;
  Array<Real> vector2;
};

class CircleRegion : public Region
{
public: 
  CircleRegion (const std::string&   r_name,
                const std::string&   r_purpose,
                const Array<Real>&   r_center,
                Real                 radius);
  virtual ~CircleRegion() {}

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;

  Array<Real> center;
  Real rsqrd;
};

#else

class SweptPolygonRegion : public Region
{
public: 
  enum PLANE {XYPLANE, YZPLANE, XZPLANE};

  SweptPolygonRegion (const std::string&   r_name,
                      const std::string&   r_purpose,
                      const Array<Real>&   v1vector,
                      const Array<Real>&   v2vector,
                      const PLANE&         plane,
                      const Array<Real>&   extent);
  virtual ~SweptPolygonRegion() {}

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;

  Array<Real> vector1;
  Array<Real> vector2;
  PLANE plane;
  Array<Real> extent;
};

class RotatedPolygonRegion : public Region
{
public: 
  enum PLANE {XYPLANE, YZPLANE, XZPLANE};

  RotatedPolygonRegion(const std::string& r_name,
                       const std::string& r_purpose,
                       const Array<Real>& v1vector,
                       const Array<Real>& v2vector,
                       const PLANE&       _plane,
                       const Array<Real>& _reference_pt,
                       const std::string& axis);
  virtual ~RotatedPolygonRegion() {}

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;

  int n, t1, t2;
  Array<Real> vector1;
  Array<Real> vector2;
  PLANE plane;
  Array<Real> reference_pt;
  std::string axis;
};
#endif

#if BL_SPACEDIM == 2
std::ostream& operator<< (std::ostream& os, const PolygonRegion& rhs);
#else
std::ostream& operator<< (std::ostream& os, const SweptPolygonRegion& rhs);
#endif

class ColorFunctionRegion : public BoxRegion
{
public: 
  
  ColorFunctionRegion (const std::string& r_name,
                       const std::string& r_purpose,
		       const std::string& file_name,
                       int                color_val);
  virtual ~ColorFunctionRegion() {delete m_color_map;}

  IntVect atIndex(Array<Real> x) const;
  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;

  void set_color_map();
  BaseFab<int>* m_color_map;
  Array<Real> dx;
  int m_color_val;
  std::string m_file;
};

std::ostream& operator<< (std::ostream& os, const ColorFunctionRegion& rhs);


class ComplementRegion : public Region
{
public:
  ComplementRegion (const std::string& r_name,
                    const std::string& r_purpose,
                    Region*            exclude_region);
  virtual ~ComplementRegion() {}
  virtual void SetExcludeRegion(Region* region) {exclude_region = region;}

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;
  Region* exclude_region;
};

class UnionRegion : public Region
{
public:
  UnionRegion (const std::string&    r_name,
               const std::string&    r_purpose,
               const Array<Region*>& regions);
  virtual ~UnionRegion() {}
  virtual void SetUnionRegions(const Array<Region*>& regions);

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;
  Array<Region*> regions;
};

class SubtractionRegion : public Region
{
public:
  SubtractionRegion (const std::string&    r_name,
                     const std::string&    r_purpose,
                     const Array<Region*>& regions);
  virtual ~SubtractionRegion() {}
  virtual void SetRegions(const Array<Region*>& regions);

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;
  Array<Region*> regions;
};

class CompoundRegion : public Region
{
public:
  enum Op {INIT, UNION, INTERSECT, SUBTRACT};

  CompoundRegion (const std::string& r_name,
                  const std::string& r_purpose,
                  Region&            src_region);
  virtual ~CompoundRegion() {}

  CompoundRegion& Union(Region& src_region);
  CompoundRegion& Intersect(Region& src_region);
  CompoundRegion& Subtract(Region& src_region);

  virtual bool inRegion (const Array<Real>& x) const;
  virtual std::ostream& print (std::ostream& os) const;
  Array<std::pair<Op,Region*> > regions;
};

class AllRegion : public BoxRegion
{
public:
  AllRegion ()
    : BoxRegion("All","any",domlo,domhi) {type="all";};
  virtual ~AllRegion() {}

  virtual bool inRegion (const Array<Real>& x) const;
};

class AllBCRegion : public BoxRegion
{

public: 
  
  AllBCRegion (int dir, int lo_or_hi);
  virtual ~AllBCRegion() {}

  virtual bool inRegion (const Array<Real>& x) const;

protected:
  int p_dir, p_lohi;
};

#endif /*_REGION_H_*/
