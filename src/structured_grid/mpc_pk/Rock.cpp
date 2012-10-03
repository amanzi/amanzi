#include <winstd.H>
#include "iostream"

#include "Rock.H"
#include "PGslib.H"

int Rock::twoexp = 1;
int Rock::max_level = 0;
Array<int>  Rock::n_cell;
Array<int>  Rock::fratio;
Array<Real> Rock::problo;
Array<Real> Rock::probhi;

std::ostream& operator<< (std::ostream& os, const Rock& rhs)
{
    rhs.operator<<(os);
    os << '\n';
}

std::map<std::string,int> Rock::rock_dist_map = Rock::create_rock_dist_map();

Rock::Rock(const std::string& name,
           Real               density,
           Real               porosity,
           int                porosity_dist_type,
           const Array<Real>& porosity_dist_param,
           const Array<Real>& permeability,
           int                permeability_dist_type,
           const Array<Real>& permeability_dist_param,
           int                krType,
           const Array<Real>& krParam,
           int                cplType,
           const Array<Real>& cplParam,
           const PArray<Region>& regions_)
    : name(name), density(density),
      porosity(porosity), porosity_dist_type(porosity_dist_type),
      porosity_dist_param(porosity_dist_param),
      permeability(permeability), permeability_dist_type(permeability_dist_type),
      permeability_dist_param(permeability_dist_param),
      krType(krType), krParam(krParam), cplType(cplType), cplParam(cplParam)
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

std::ostream&
Rock::operator<< (std::ostream& os) const
{
    os << "Rock:\n";
    os << "  name: " << name << '\n';
    os << "  regions: ";
    for (int i=0; i<regions.size(); ++i) {
        os << regions[i].name << " ";
    }
    os << '\n';
    os << "  density: " << density << '\n';

    std::string porosity_dist_name;
    for (std::map<std::string,int>::const_iterator it = rock_dist_map.begin(); it!=rock_dist_map.end(); ++it) {
        if (it->second == porosity_dist_type) {
            porosity_dist_name = it->first;
        }
    }
    os << "  porosity_dist_type: " << porosity_dist_name << " (" << porosity_dist_type << ")\n";
    if (porosity_dist_name=="uniform") {
        os << "    value: " << porosity << '\n';
    } else if (porosity_dist_name=="random"
               || porosity_dist_name=="geostatistic") {
        os << "    porosity_dist_param: ";
        for (int i=0; i<porosity_dist_param.size(); ++i) {
            os << porosity_dist_param[i] << " ";
        }
        os << '\n';
    }

    std::string permeability_dist_name;
    for (std::map<std::string,int>::const_iterator it = rock_dist_map.begin(); it!=rock_dist_map.end(); ++it) {
        if (it->second == permeability_dist_type) {
            permeability_dist_name = it->first;
        }
    }

    os << "  permeability_dist_type: " << permeability_dist_name << " (" << permeability_dist_type << ")\n";
    if (permeability_dist_name=="uniform") {
        os << "    value: ";
        for (int i=0; i<permeability.size(); ++i) {
            os << permeability[i] << " ";
        }
        os << '\n';
    } else if (permeability_dist_name=="random"
               || permeability_dist_name=="geostatistic") {
        os << "    permeability_dist_param: ";
        for (int i=0; i<permeability_dist_param.size(); ++i) {
            os << permeability_dist_param[i] << " ";
        }
        os << '\n';
    }
    os << "  krType: " << krType << '\n';
    os << "    krParam: ";
    for (int i=0; i<krParam.size(); ++i) {
        os << krParam[i] << " ";
    }
    os << '\n';
    os << "  cplType: " << cplType << '\n';
    os << "    cplParam: ";
    for (int i=0; i<cplParam.size(); ++i) {
        os << cplParam[i] << " ";
    }
    os << '\n';
}

BoxArray
Rock::ba_for_finest_data(int         max_level, 
                         Array<int>& n_cell,
                         Array<int>& fratio,
                         int         maxBaseGrid,
                         int         nGrow)
{
  //
  // Create grids at finest level.
  //

  Box bx(IntVect::TheZeroVector(),
         IntVect(D_DECL(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1)));
  bx.grow(nGrow);

  BL_ASSERT(fratio.size()>=max_level);
  twoexp = 1;
  for (int ii = 1; ii<=max_level;ii++) {
    twoexp *= fratio[ii-1];
  }
  bx.refine(twoexp);
  
  BoxArray ba(bx);
  ba.maxSize(maxBaseGrid);

  return ba;
}

void Rock::build_kmap(MultiFab&             mfdata, 
		      const std::string&    gsfile) const
{
    int nGrow = mfdata.nGrow();
    if (permeability_dist_type == rock_dist_map["uniform"]) {

        for (int i = 0; i<regions.size(); i++)
        {
            set_constant_kval(mfdata,regions[i],nGrow);
        }
    }
    else if (permeability_dist_type == rock_dist_map["random"])
    {
        for (int i = 0; i<regions.size(); i++)
        {
            PGslib::parRand(permeability,
                            permeability_dist_param[0],
                            n_cell,
                            twoexp,
                            mfdata);
	}
    }
    else if (permeability_dist_type == rock_dist_map["geostatistic"])
    {
        BL_ASSERT(gsfile!="");
        PGslib::rdpGaussianSim(permeability,
                               permeability_dist_param[0],
                               n_cell,
                               problo,
                               probhi,
                               twoexp,
                               mfdata,
                               gsfile); 
    }
    else {
        BoxLib::Abort("Unsupported permeability option");
    }
}

void Rock::build_pmap(MultiFab&          mfdata, 
		      const std::string& gsfile) const
{
    int nGrow = mfdata.nGrow();
    Array<Real> porosity_array(1,porosity);
    if (porosity_dist_type == rock_dist_map["uniform"])
    {
        for (int i= 0; i<regions.size(); i++)
        {
            set_constant_pval(mfdata,regions[i],nGrow);
        }
    }
    else if (porosity_dist_type == rock_dist_map["random"])
    {
        PGslib::parRand(porosity_array,
                        porosity_dist_param[0],
                        n_cell,
                        twoexp,
                        mfdata);
    } 
    else if (porosity_dist_type == rock_dist_map["geostatistic"])
    {
        BL_ASSERT(gsfile!="");
        PGslib::rdpGaussianSim(porosity_array,
                               porosity_dist_param[0],
                               n_cell,
                               problo,
                               probhi,
                               twoexp,
                               mfdata,
                               gsfile);
    }
    else {
        BoxLib::Abort("Unsupported permeability option");
    }
}



void Rock::set_constant_kval(MultiFab&     mfdata, 
			     const Region& region_local,
                             int           nGrow) const
{
    set_constant_val(mfdata,region_local,permeability,nGrow);
}

void Rock::set_constant_pval(MultiFab&     mfdata, 
			     const Region& region_local,
                             int           nGrow) const
{
    Array<Real> porosity_array(1,porosity);
    set_constant_val(mfdata,region_local,porosity_array,nGrow);
}

void Rock::set_constant_kval(FArrayBox&  fab, 
			     const Real* dx) const
{
  for (int ir=0; ir<regions.size(); ir++)
    regions[ir].setVal(fab,permeability,dx,0,0,permeability.size());
}

void Rock::set_constant_pval(FArrayBox&  fab, 
			     const Real* dx) const
{
  Array<Real> porosity_array(1,porosity);
  for (int ir=0; ir<regions.size(); ir++)
    regions[ir].setVal(fab,porosity_array,dx,0,0,1);
}

void Rock::set_constant_krval(FArrayBox&  fab, 
			      const Real* dx) const
{
  int nval = krParam.size()+1;
  BL_ASSERT(fab.nComp() >= nval);
  Array<Real> param_tmp(nval);
  param_tmp[0] = krType;
  for (int i=0;i<krParam.size();i++)
    param_tmp[i+1] = krParam[i];
  for (int ir=0; ir<regions.size();ir++)
    regions[ir].setVal(fab,param_tmp,dx,0,0,nval);
}

void Rock::set_constant_cplval(FArrayBox&  fab, 
			       const Real* dx) const
{
  int nval = cplParam.size()+1;
  BL_ASSERT(fab.nComp() >= nval);
  Array<Real> param_tmp(nval);
  param_tmp[0] = cplType;
  for (int i=0;i<cplParam.size();i++)
    param_tmp[i+1] = cplParam[i];
  for (int ir=0; ir<regions.size();ir++)
    regions[ir].setVal(fab,param_tmp,dx,0,0,nval); 
}

void Rock::set_constant_val(MultiFab&          mfdata, 
			    const Region&      region_local,
			    const Array<Real>& val,
                            int                nGrow) const
{
    Array<Real> dx(BL_SPACEDIM);
    int ng_twoexp = nGrow*twoexp;
    for (int i=0;i<BL_SPACEDIM; i++)  
        dx[i] = (probhi[i]-problo[i])/(n_cell[i]*twoexp);
  
    for (MFIter mfi(mfdata); mfi.isValid(); ++mfi)
        region_local.setVal(mfdata[mfi],val,dx.dataPtr(),ng_twoexp,0,val.size());
}

void Rock::set_constant_val(FArrayBox&         fab, 
			    const Region&      region_local,
			    const Array<Real>& val,
                            int                nGrow) const
{
    Array<Real> dx(BL_SPACEDIM);
    int ng_twoexp = nGrow*twoexp;
    for (int i=0;i<BL_SPACEDIM; i++)  
        dx[i] = (probhi[i]-problo[i])/(n_cell[i]*twoexp);

    region_local.setVal(fab,val,dx.dataPtr(),ng_twoexp,0,val.size());
}
