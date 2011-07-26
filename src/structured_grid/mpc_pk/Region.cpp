#include <winstd.H>
#include "iostream"

#include "Region.H"

Region::Region (std::string r_name, int r_purpose, int r_type)
  : name(r_name),
    purpose(r_purpose),
    type(r_type)
{
}

void
Region::setVal(FArrayBox&   fab, 
	       Array<Real>& val, 
	       const Real*  dx, 
	       int          ng,
	       int          scomp,
	       int          ncomp)
{
  int nval = val.size();
  BL_ASSERT(ncomp-scomp <= nval);

  Array<Real> x(BL_SPACEDIM);
  const int* lo = fab.loVect();
  const int* hi = fab.hiVect();


#if (BL_SPACEDIM == 2)	
  for (int iy=lo[1]; iy<hi[1]+1; iy++) 
    {
      x[1] = dx[1]*(iy+0.5);
      for (int ix=lo[0]; ix<hi[0]+1; ix++) 
	{
	  x[0] = dx[0]*(ix+0.5);
	  if (inregion(x))
	    {
	      for (int n=scomp; n<ncomp;n++)
		fab(IntVect(ix,iy),n) = val[n-scomp];
	    }
	}
    }
#else
  for (int iz=lo[2]; iz<hi[2]+1; iz++) 
    {
      x[2] = dx[2]*(iz+0.5);
      for (int iy=lo[1]; iy<hi[1]+1; iy++) 
	{
	  x[1] = dx[1]*(iy+0.5);
	  for (int ix=lo[0]; ix<hi[0]+1; ix++) 
	    {
	      x[0] = dx[0]*(ix+0.5);
	      if (inregion(x))
		for (int n=scomp; n<ncomp;n++)
		  fab(IntVect(ix,iy,iz),n) = val[n-scomp];
		
	    }
	}
    }
#endif
}

void
Region::setVal(FArrayBox&  fab, 
	       const Real  val, 
	       const int   idx, 
	       const Real* dx, 
	       int         ng)
{
  Array<Real> x(BL_SPACEDIM);

  const int* lo = fab.loVect();
  const int* hi = fab.hiVect();

#if (BL_SPACEDIM == 2)	
  for (int iy=lo[1]; iy<hi[1]+1; iy++) 
    {
      x[1] = dx[1]*(iy+0.5);
      for (int ix=lo[0]; ix<hi[0]+1; ix++) 
	{
	  x[0] = dx[0]*(ix+0.5);
	  if (inregion(x))
	    fab(IntVect(ix,iy),idx) = val;
	}
    }
#else
  for (int iz=lo[2]; iz<hi[2]+1; iz++) 
    {
      x[2] = dx[2]*(iz+0.5);
      for (int iy=lo[1]; iy<hi[1]+1; iy++) 
	{
	  x[1] = dx[1]*(iy+0.5);
	  for (int ix=lo[0]; ix<hi[0]+1; ix++) 
	    {
	      x[0] = dx[0]*(ix+0.5);
	      if (inregion(x))
		fab(IntVect(ix,iy,iz),idx) = val;
		
	    }
	}
    }
#endif
}

void
boxRegion::set (Array<Real>& param)
{

  for (int i=0; i<BL_SPACEDIM; i++)
    {
      vertex_lo[i] = param[i];
      vertex_hi[i] = param[BL_SPACEDIM+i];
    }
}

bool 
boxRegion::inregion(Array<Real>& x)
{
  bool inflag = false;
#if (BL_SPACEDIM == 2)
  if (x[0]>=vertex_lo[0] && x[0]<=vertex_hi[0] &&
      x[1]>=vertex_lo[1] && x[1]<=vertex_hi[1])
    inflag = true;
#else
  if (x[0]>=vertex_lo[0] && x[0]<=vertex_hi[0] &&
      x[1]>=vertex_lo[1] && x[1]<=vertex_hi[1] &&
      x[2]>=vertex_lo[2] && x[2]<=vertex_hi[2])
    inflag = true;
#endif
  return inflag;
}

void
allRegion::set (Array<Real>& param)
{

  for (int i=0; i<BL_SPACEDIM; i++)
    {
      vertex_lo[i] = param[i];
      vertex_hi[i] = param[BL_SPACEDIM+i];
    }
}

bool 
allRegion::inregion(Array<Real>& x)
{
  return true;
}


allBCRegion::allBCRegion (int dir, int lo_or_hi)
  : Region("BC",4,101)
{
  p_dir  = dir;
  p_lohi = lo_or_hi;
  if (dir == 0 && lo_or_hi == 0)
    name = "XLOBC";
  else if (dir == 0 && lo_or_hi == 1)
    name = "XHIBC";
  else if (dir == 1 && lo_or_hi == 0)
    name = "YLOBC";
  else if (dir == 1 && lo_or_hi == 1)
    name = "YHIBC";
#if BL_SPACEDIM == 3
  else if (dir == 2 && lo_or_hi == 0)
    name = "ZLOBC";
  else if (dir == 2 && lo_or_hi == 1)
    name = "ZHIBC";
#endif  
}

void
allBCRegion::set (Array<Real>& param)
{
  for (int i=0; i<BL_SPACEDIM; i++)
    {
      vertex_lo[i] = param[i];
      vertex_hi[i] = param[BL_SPACEDIM+i];
    }

  if (p_lohi == 0)
    vertex_hi[p_dir] = vertex_lo[p_dir];
  else 
    vertex_lo[p_dir] = vertex_hi[p_dir];
}

bool 
allBCRegion::inregion(Array<Real>& x)
{
  bool inflag = false;

  if ((p_lohi == 0 && x[p_dir] <= vertex_lo[p_dir]) ||
      (p_lohi == 1 && x[p_dir] >= vertex_hi[p_dir])) 
    inflag = true;

  return inflag;
}
