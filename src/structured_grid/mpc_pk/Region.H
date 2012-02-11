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

  Region(std::string r_name, std::string r_purpose, std::string r_type);

  void setVal(FArrayBox&         fab, 
	      const Array<Real>& val, 
	      const Real*        dx, 
	      int                ng,
	      int                scomp,
	      int                ncomp) const;

  void setVal(FArrayBox&  fab,
              const Real  val, 
              const int   idx,
	      const Real* dx, 
              int         ng) const;

  virtual bool inregion (const Array<Real>& x) const = 0;
  virtual std::ostream& operator<< (std::ostream& os) const;

  // Identifier for the region 
  std::string name;
  std::string purpose;
  std::string type;

  static Real geometry_eps;
};

std::ostream& operator<< (std::ostream& os, const Region& rhs);

class pointRegion : public Region
{

public: 
  
    pointRegion (std::string r_name, std::string r_purpose, std::string r_type, const Array<Real>& coor)
      : Region(r_name,r_purpose,r_type), coor(coor) {} ;

  virtual bool inregion (const Array<Real>& x) const;
  virtual std::ostream& operator<< (std::ostream& os) const;

  Array<Real> coor;
};

std::ostream& operator<< (std::ostream& os, const pointRegion& rhs);

class boxRegion : public Region 
{

public: 
  
    boxRegion (std::string r_name, std::string r_purpose, std::string r_type,
               const Array<Real>& lo, const Array<Real>& hi)
      : Region(r_name,r_purpose,r_type), lo(lo), hi(hi) {} ;

  virtual bool inregion (const Array<Real>& x) const;
  virtual std::ostream& operator<< (std::ostream& os) const;

  Array<Real> lo;
  Array<Real> hi;
};

std::ostream& operator<< (std::ostream& os, const boxRegion& rhs);

class colorFunctionRegion : public boxRegion 
{

public: 
  
  colorFunctionRegion (std::string r_name,
                       std::string r_purpose,
                       std::string r_type,
		       std::string file_name,
                       int         color_val);

  IntVect atIndex(Array<Real> x) const;
  virtual bool inregion (const Array<Real>& x) const;
  virtual std::ostream& operator<< (std::ostream& os) const;

  void set_color_map();
  BaseFab<int>* m_color_map;
  Array<Real> dx;
  int m_color_val;
  std::string m_file;
};

std::ostream& operator<< (std::ostream& os, const colorFunctionRegion& rhs);

#if 1
class allRegion : public boxRegion 
{

public: 
  
    allRegion (const Array<Real>& lo,
               const Array<Real>& hi)
        : boxRegion("ALL","any","all",lo,hi) {};

  virtual bool inregion (const Array<Real>& x) const;
};

class allBCRegion : public boxRegion 
{

public: 
  
    allBCRegion (int dir, int lo_or_hi,
                 const Array<Real>& lo,
                 const Array<Real>& hi);

  virtual bool inregion (const Array<Real>& x) const;

protected:
  int p_dir, p_lohi;
};
#endif

#endif /*_REGION_H_*/
