#include <RegionManager.H>

#include <ParmParse.H>

RegionManager::RegionManager()
{
  ParmParse pp("geometry");

  Array<Real> problo, probhi;
  pp.getarr("prob_lo",problo,0,BL_SPACEDIM);
  pp.getarr("prob_hi",probhi,0,BL_SPACEDIM);
  Region::domlo = problo;
  Region::domhi = probhi;
  
  Real geometry_eps = -1; pp.get("geometry_eps",geometry_eps);
  Region::geometry_eps = geometry_eps;

  // set up  1+2*BL_SPACEDIM default regions
  bool generate_default_regions = true; pp.query("generate_default_regions",generate_default_regions);
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

  // Get parameters for each user defined region 
  int nregion = nregion_DEF;

  int nregion_user = pp.countval("regions");

  if (!generate_default_regions  && nregion_user==0) {
    BoxLib::Abort("Default regions not generated and none provided.  Perhaps omitted regions list?");
  }

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
      else {
	std::string m = "region type not supported \"" + r_type + "\"";
	BoxLib::Abort(m.c_str());
      }
    }
  }
  
  for (int i=0; i<regions.size(); ++i) {
    name_to_region_idx[regions[i]->name] = i;
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

