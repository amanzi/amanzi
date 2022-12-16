/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <RegionManager.H>

#include <ParmParse.H>

static Real GEOMETRY_EPS_DEF = 1.e-10;

RegionManager::RegionManager()
{
  ParmParse pp("geometry");

  Array<Real> problo, probhi;
  pp.getarr("prob_lo",problo,0,BL_SPACEDIM);
  pp.getarr("prob_hi",probhi,0,BL_SPACEDIM);

  Init(problo,probhi);
}

RegionManager::RegionManager(const Array<Real>& plo,
                             const Array<Real>& phi)
{
  Init(plo,phi);
}

void
RegionManager::Init(const Array<Real>& plo,
                    const Array<Real>& phi)
{
  Region::domlo = plo;
  Region::domhi = phi;

  ParmParse pp("geometry");

  Real geometry_eps = GEOMETRY_EPS_DEF; pp.query("geometry_eps",geometry_eps);
  Region::geometry_eps = geometry_eps;

  // set up  1+2*BL_SPACEDIM default regions
  bool generate_default_regions = false; pp.query("generate_default_regions",generate_default_regions);
  int nregion_DEF = 0;
  regions.clear();
  if (generate_default_regions) {
      nregion_DEF = 1 + 2*BL_SPACEDIM;
      regions.resize(nregion_DEF);
      regions[0] = new   AllRegion();
      regions[1] = new AllBCRegion(0,0);
      regions[2] = new AllBCRegion(0,1);
      regions[3] = new AllBCRegion(1,0);
      regions[4] = new AllBCRegion(1,1);
#if BL_SPACEDIM == 3
      regions[5] = new AllBCRegion(2,0);
      regions[6] = new AllBCRegion(2,1);
#endif
  }

  int region_ctr = 0;
  for (int i=0; i<regions.size(); ++i) {
    const std::string name = regions[i]->name;
    name_to_region_idx[name] = region_ctr++;
  }

  // Get parameters for each user defined region
  int nregion = nregion_DEF;

  int nregion_user = pp.countval("regions");

  if (!generate_default_regions  && nregion_user==0) {
    BoxLib::Abort("Default regions not generated and none provided.  Perhaps omitted regions list?");
  }

  std::map<std::string,std::string> complement_exclude_regions;
  std::map<std::string,Array<std::string> > union_regions;
  std::map<std::string,Array<std::string> > subtraction_regions;

  if (nregion_user) {
    std::string r_purpose, r_type;
    Array<std::string> r_name;
    pp.getarr("regions",r_name,0,nregion_user);
    nregion += nregion_user;
    regions.resize(nregion);

    for (int j=0; j<nregion_user; ++j) {
      const std::string prefix("geometry." + r_name[j]);
      ParmParse ppr(prefix.c_str());
      ppr.get("purpose",r_purpose);
      ppr.get("type",r_type);

      if (r_type == "point") {
	Array<Real> coor;
	ppr.getarr("coordinate",coor,0,BL_SPACEDIM);
	regions[nregion_DEF+j] = new PointRegion(r_name[j],r_purpose,coor);
      }
      else if (r_type == "box" || r_type == "surface") {
	Array<Real> lo_coor,hi_coor;
	ppr.getarr("lo_coordinate",lo_coor,0,BL_SPACEDIM);
	ppr.getarr("hi_coordinate",hi_coor,0,BL_SPACEDIM);
	regions[nregion_DEF+j] = new BoxRegion(r_name[j],r_purpose,lo_coor,hi_coor);
      }
      else if (r_type == "color_function") {
	int color_value; ppr.get("color_value",color_value);
	std::string color_file; ppr.get("color_file",color_file);
	ColorFunctionRegion* cfr = new ColorFunctionRegion(r_name[j],r_purpose,color_file,color_value);
	regions[nregion_DEF+j] = cfr;
      }
      else if (r_type == "logical") {
        std::string operation; ppr.get("operation",operation);
        if (operation=="complement") {
          std::string exclude_region; ppr.get("region",exclude_region);
          complement_exclude_regions[r_name[j]] = exclude_region;
          ComplementRegion* cfr = new ComplementRegion(r_name[j],r_purpose,0); // Delay setting exclude until finished
          regions[nregion_DEF+j] = cfr;
        }
        else if (operation=="union") {
          int nr = ppr.countval("regions");
          ppr.getarr("regions",union_regions[r_name[j]],0,nr);
          UnionRegion* cfr = new UnionRegion(r_name[j],r_purpose,Array<const Region*>()); // Delay setting regions until finished
          regions[nregion_DEF+j] = cfr;
        }
        else if (operation=="subtraction") {
          int nr = ppr.countval("regions");
          ppr.getarr("regions",subtraction_regions[r_name[j]],0,nr);
          SubtractionRegion* cfr = new SubtractionRegion(r_name[j],r_purpose,Array<const Region*>()); // Delay setting regions until finished
          regions[nregion_DEF+j] = cfr;
        }
        else {
          BoxLib::Abort("Unsupported logical region operation");
        }
      }
#if BL_SPACEDIM==2
      else if (r_type == "polygon") {
        int nv = ppr.countval("v1");
        BL_ASSERT(nv>0);
        Array<Real> v1(nv), v2(nv);
        ppr.getarr("v1",v1,0,nv);
        ppr.getarr("v2",v2,0,nv);
        PolygonRegion* cfr = new PolygonRegion(r_name[j],r_purpose,v1,v2);
	regions[nregion_DEF+j] = cfr;
      }
      else if (r_type == "ellipse") {
        int nc = ppr.countval("center");
        BL_ASSERT(nc>=BL_SPACEDIM);
        Array<Real> c(nc); ppr.getarr("center",c,0,nc);

        int nr = ppr.countval("radius");
        BL_ASSERT(nr>=BL_SPACEDIM);
        Array<Real> r(nr); ppr.getarr("radius",r,0,nr);

        EllipseRegion* cfr = new EllipseRegion(r_name[j],r_purpose,c,r);
	regions[nregion_DEF+j] = cfr;
      }
#else
      else if (r_type == "swept_polygon") {
        std::string plStr; ppr.get("plane",plStr);
        BL_ASSERT(plStr == "XY" || plStr == "YZ" || plStr == "XZ");
        SweptPolygonRegion::PLANE plane = (plStr == "XY" ?
                                           SweptPolygonRegion::XYPLANE :
                                           (plStr == "YZ" ?
                                            SweptPolygonRegion::YZPLANE :
                                            SweptPolygonRegion::XZPLANE) );
        Array<Real> extent; ppr.getarr("extent",extent,0,2);
        int nv = ppr.countval("v1");
        BL_ASSERT(nv>0);
        Array<Real> v1(nv), v2(nv);
        ppr.getarr("v1",v1,0,nv);
        ppr.getarr("v2",v2,0,nv);
        SweptPolygonRegion* cfr = new SweptPolygonRegion(r_name[j],r_purpose,v1,v2,plane,extent);
	regions[nregion_DEF+j] = cfr;
      }
      else if (r_type == "rotated_polygon") {
        std::string plStr; ppr.get("plane",plStr);
        BL_ASSERT(plStr == "XY" || plStr == "YZ" || plStr == "XZ");
        RotatedPolygonRegion::PLANE plane = (plStr == "XY" ?
                                             RotatedPolygonRegion::XYPLANE :
                                             (plStr == "YZ" ?
                                              RotatedPolygonRegion::YZPLANE :
                                              RotatedPolygonRegion::XZPLANE) );
        int nv = ppr.countval("v1");
        BL_ASSERT(nv>0);
        Array<Real> v1(nv), v2(nv);
        ppr.getarr("v1",v1,0,nv);
        ppr.getarr("v2",v2,0,nv);

        Array<Real> refPt; ppr.getarr("reference_pt",refPt,0,BL_SPACEDIM);
        std::string axStr; ppr.get("axis",axStr);
        BL_ASSERT(axStr == "X" || axStr == "Y" || axStr == "Z");
        RotatedPolygonRegion* cfr = new RotatedPolygonRegion(r_name[j],r_purpose,v1,v2,plane,refPt,axStr);
	regions[nregion_DEF+j] = cfr;
      }
#endif
      else {
	std::string m = "region type not supported \"" + r_type + "\"";
	BoxLib::Abort(m.c_str());
      }
      name_to_region_idx[r_name[j]] = region_ctr++;
    }
  }

  for (std::map<std::string,std::string>::const_iterator it = complement_exclude_regions.begin();
       it != complement_exclude_regions.end(); ++it) {
    const std::string& region_name = it->first;
    const std::string& excl_region_name = it->second;
    if (name_to_region_idx.count(region_name)==0) {
      const std::string msg = "Attempted to buid complement region with undefined component region "+region_name;
      BoxLib::Abort(msg.c_str());
    }
    ComplementRegion* region = dynamic_cast<ComplementRegion*>(regions[name_to_region_idx[region_name]]);
    BL_ASSERT(region != 0);
    if (name_to_region_idx.count(excl_region_name)==0) {
      const std::string msg = "Attempted to buid complement region with undefined component region "+excl_region_name;
      BoxLib::Abort(msg.c_str());
    }
    region->SetExcludeRegion(regions[name_to_region_idx[excl_region_name]]);
  }

  for (std::map<std::string,Array<std::string> >::const_iterator it = union_regions.begin();
       it != union_regions.end(); ++it) {
    const std::string& region_name = it->first;
    const Array<std::string>& union_region_names = it->second;
    if (name_to_region_idx.count(region_name)==0) {
      const std::string msg = "Attempted to build union region with undefined component region "+region_name;
      BoxLib::Abort(msg.c_str());
    }
    UnionRegion* region = dynamic_cast<UnionRegion*>(regions[name_to_region_idx[region_name]]);
    BL_ASSERT(region != 0);
    Array<const Region*> rp(union_region_names.size());
    for (int i=0; i<rp.size(); ++i) {
      if (name_to_region_idx.count(union_region_names[i])==0) {
        const std::string msg = "Attempted to buid union region with undefined component region "+union_region_names[i];
        BoxLib::Abort(msg.c_str());
      }
      rp[i] = regions[name_to_region_idx[union_region_names[i]]];
    }
    region->SetUnionRegions(rp);
  }

  for (std::map<std::string,Array<std::string> >::const_iterator it = subtraction_regions.begin();
       it != subtraction_regions.end(); ++it) {
    const std::string& region_name = it->first;
    const Array<std::string>& subtraction_region_names = it->second;
    if (name_to_region_idx.count(region_name)==0) {
      const std::string msg = "Attempted to buid subtraction region with undefined component region "+region_name;
      BoxLib::Abort(msg.c_str());
    }
    SubtractionRegion* region = dynamic_cast<SubtractionRegion*>(regions[name_to_region_idx[region_name]]);
    BL_ASSERT(region != 0);
    Array<const Region*> rp(subtraction_region_names.size());
    for (int i=0; i<rp.size(); ++i) {
      if (name_to_region_idx.count(subtraction_region_names[i])==0) {
        const std::string msg = "Attempted to buid complement region with undefined component region "+subtraction_region_names[i];
        BoxLib::Abort(msg.c_str());
      }
      rp[i] = regions[name_to_region_idx[subtraction_region_names[i]]];
    }
    region->SetRegions(rp);
  }
}

RegionManager::~RegionManager()
{
  Clear();
}

void
RegionManager::Clear()
{
  for (int i=0; i<regions.size(); ++i) {
    delete regions[i];
  }
  regions.resize(0);
}

Array<const Region*>
RegionManager::RegionPtrArray() const
{
  Array<const Region*> ret(regions.size());
  for (int i=0; i<regions.size(); ++i) {
    ret[i] = regions[i];
  }
  return ret;
}

Array<const Region*>
RegionManager::RegionPtrArray(const Array<std::string>& region_names) const
{
  if (region_names.size()==0) {
    return RegionPtrArray();
  }

  Array<const Region*> ret(region_names.size());
  for (int i=0; i<region_names.size(); ++i) {
    const std::string& name = region_names[i];
    std::map<std::string,int>::const_iterator it = name_to_region_idx.find(name);
    if (it != name_to_region_idx.end()) {
      ret[i] = regions[it->second];
    }
    else {
      std::string m = "Named region not found: \"" + name + "\"";
      BoxLib::Error(m.c_str());
    }
  }
  return ret;
}

std::ostream& operator<<(std::ostream& os, const RegionManager& rm)
{
  os << "RegionManager: " << '\n';
  const Array<const Region*> regions = rm.RegionPtrArray();
  for (int i=0; i<regions.size(); ++i) {
    os << *regions[i] << '\n';
  }
  return os;
}
