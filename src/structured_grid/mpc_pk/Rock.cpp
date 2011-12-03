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

std::map<std::string,int> Rock::rock_dist_map = Rock::create_rock_dist_map();

BoxArray
Rock::build_finest_data(int        max_level, 
			Array<int>& n_cell,
			Array<int>& fratio,
			int        maxBaseGrid)
{
  //
  // Create grids at finest level.
  //

  Box bx;
#if (BL_SPACEDIM == 3)
  bx = Box(IntVect(0,0,0),
	   IntVect(n_cell[0]-1,n_cell[1]-1,n_cell[2]-1));
#else
  bx = Box(IntVect(0,0),
	   IntVect(n_cell[0]-1,n_cell[1]-1));
#endif
  bx.grow(3);

  twoexp = 1;
  if (max_level > 0) 
    {
      twoexp = fratio[0];
      for (int ii = 1; ii<max_level;ii++)
	twoexp *= fratio[0];
    }
  bx.refine(twoexp);

  BoxArray ba(bx);
  ba.maxSize(maxBaseGrid);

  return ba;
}

void Rock::build_kmap(MultiFab&      mfdata, 
		      Array<Region*>& region_array, 
		      std::string&    gsfile)
{
  for (int i = 0; i<region.size(); i++)
    {
      if (porosity_dist_type == rock_dist_map["uniform"])
	set_constant_kval(mfdata,region_array[region[i]]);
      else if (porosity_dist_type == rock_dist_map["random"])
	{
	  PGslib::parRand(permeability,
			  permeability_dist_param[0],
			  n_cell,
			  twoexp,
			  mfdata);
	}
      else if (porosity_dist_type == rock_dist_map["geostatistic"])
	{
	  gsfile = "sgsim.inc";
	  PGslib::rdpGaussianSim(permeability,
				 permeability_dist_param[0],
				 n_cell,
				 problo,
				 probhi,
				 twoexp,
				 mfdata,
				 gsfile); 
	}
    }
}

void Rock::build_pmap(MultiFab& mfdata, 
		      Array<Region*>& region_array, 
		      std::string& gsfile)
{
  Array<Real> porosity_array(1);
  porosity_array[0] = porosity;
  for (int i= 0; i<region.size(); i++)
    {
      if (porosity_dist_type == rock_dist_map["uniform"])
	set_constant_pval(mfdata,region_array[region[i]]);
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
	  gsfile = "sgsim.inc";
	  PGslib::rdpGaussianSim(porosity_array,
				 porosity_dist_param[0],
				 n_cell,
				 problo,
				 probhi,
				 twoexp,
				 mfdata,
				 gsfile);
	}
    }
}

void Rock::set_constant_kval(MultiFab& mfdata, 
			     Region*   region_local)
{
  set_constant_val(mfdata,region_local,permeability);
}

void Rock::set_constant_pval(MultiFab& mfdata, 
			     Region*   region_local)
{
  Array<Real> porosity_array(1);
  porosity_array[0] = porosity;
  set_constant_val(mfdata,region_local,porosity_array);
}

void Rock::set_constant_krval(FArrayBox& fab, 
			      Array<Region*>& region_array,
			      const Real* dx)
{
  int nval = krParam.size()+1;
  BL_ASSERT(fab.nComp() >= nval);
  Array<Real> param_tmp(nval);
  param_tmp[0] = krType;
  for (int i=0;i<krParam.size();i++)
    param_tmp[i+1] = krParam[i];
  for (int ir=0; ir<region.size();ir++)
    region_array[region[ir]]->setVal(fab,param_tmp,dx,0,0,nval);
}
void Rock::set_constant_cplval(FArrayBox& fab, 
			       Array<Region*>& region_array,
			       const Real* dx)
{
  int nval = cplParam.size()+1;
  BL_ASSERT(fab.nComp() >= nval);
  Array<Real> param_tmp(nval);
  param_tmp[0] = cplType;
  for (int i=0;i<cplParam.size();i++)
    param_tmp[i+1] = cplParam[i];
  for (int ir=0; ir<region.size();ir++)
    region_array[region[ir]]->setVal(fab,param_tmp,dx,0,0,nval); 
}

void Rock::set_constant_val(MultiFab&    mfdata, 
			    Region*      region_local,
			    Array<Real>& val)
{
  Real dx[BL_SPACEDIM];
  int ng_twoexp = 3*twoexp;
  for (int i=0;i<BL_SPACEDIM; i++)  
    dx[i] = (probhi[i]-problo[i])/(n_cell[i]*twoexp);
  
  for (MFIter mfi(mfdata); mfi.isValid(); ++mfi)
    region_local->setVal(mfdata[mfi],val,dx,ng_twoexp,0,val.size());
}
