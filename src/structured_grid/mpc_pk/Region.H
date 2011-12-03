#ifndef _REGION_H_
#define _REGION_H_

#include <string>
#include <Array.H>
#include <FArrayBox.H>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <LevelBld.H>
#include <AmrLevel.H>

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

    return m;
  };

  static std::map<std::string,int> region_map;

  Region(std::string r_name, std::string r_purpose, std::string r_type);

  void setVal(FArrayBox& fab, 
	      Array<Real>& val, 
	      const Real* dx, 
	      int ng,
	      int scomp,
	      int ncomp);

  void setVal(FArrayBox& fab, const Real val, const int idx,
	      const Real* dx, int ng);

  virtual void set (Array<Real>& param)  = 0;
  virtual bool inregion (Array<Real>& x) = 0;


  // Identifier for the region 
  std::string name;
  std::string purpose;
  std::string type;
};

class pointRegion : public Region
{

public: 
  
  pointRegion (std::string r_name, std::string r_purpose, std::string r_type)
    : Region(r_name,r_purpose,r_type) {} ;

  virtual void set (Array<Real>& param);
  virtual bool inregion (Array<Real>& x);

  Real coor[BL_SPACEDIM];

};

class boxRegion : public Region 
{

public: 
  
  boxRegion (std::string r_name, std::string r_purpose, std::string r_type)
    : Region(r_name,r_purpose,r_type) {} ;

  virtual void set (Array<Real>& param);
  virtual bool inregion (Array<Real>& x);

  Real vertex_lo[BL_SPACEDIM];
  Real vertex_hi[BL_SPACEDIM];

};

class allRegion : public Region 
{

public: 
  
  allRegion ()
    : Region("ALL","any","all") {};

  virtual void set (Array<Real>& param);
  virtual bool inregion (Array<Real>& x);

  Real vertex_lo[BL_SPACEDIM];
  Real vertex_hi[BL_SPACEDIM];

};

class allBCRegion : public Region 
{

public: 
  
  allBCRegion (int dir, int lo_or_hi);

  virtual void set (Array<Real>& param);
  virtual bool inregion (Array<Real>& x);

  Real vertex_lo[BL_SPACEDIM];
  Real vertex_hi[BL_SPACEDIM];

private:
  int p_dir;
  int p_lohi;

};


#endif /*_REGION_H_*/
