#include <winstd.H>
#include "iostream"

#include "Region.H"

std::map<std::string,int> Region::region_map = Region::create_region_map();

Region::Region (std::string r_name, 
		std::string r_purpose, 
		std::string r_type)
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
                for (int n=scomp; n<ncomp;n++) {
                    fab(IntVect(ix,iy),n) = val[n-scomp];
                }
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
	  if (inregion(x)) {
	    fab(IntVect(ix,iy),idx) = val;
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
	      if (inregion(x)) {
		fab(IntVect(ix,iy,iz),idx) = val;
              }		
	    }
	}
    }
#endif
}

void 
pointRegion::set (Array<Real>& param)
{
  for (int i = 0; i<BL_SPACEDIM; i++)
    coor[i] = param[i];
}

bool 
pointRegion::inregion (Array<Real>& x)
{
  bool inflag = false;
#if (BL_SPACEDIM == 2)
  if (x[0]>=coor[0] && x[1]>=coor[1])
    inflag = true;
#else
  if (x[0]>=coor[0] && x[1]>=coor[1] && x[2]>=coor[2] )
    inflag = true;
#endif
  return inflag;
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

colorFunctionRegion::colorFunctionRegion (std::string r_name,
                                          std::string r_purpose,
                                          std::string r_type,
					  std::string file_name,
                                          int color_val)
  : boxRegion(r_name,r_purpose,r_type), 
    dx(BL_SPACEDIM),
    m_color_val(color_val)
{
  set_color_map(file_name);
}

void
colorFunctionRegion::set_color_map(const std::string& file_name)
{
  Array<char> fileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(file_name, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream is(fileCharPtrString, std::istringstream::in);
  
  int datatype;
  is >> datatype;
  if (datatype != 0) {
    BoxLib::Abort("color function: require DATATYPE == 0");
  }

  std::string gridtype;
  is >> gridtype;
  if (!is.good()) {
    BoxLib::Abort("color function: error reading gridtype");
  }
  std::string ok_grid_type = D_PICK("1DCoRectMesh","2DCoRectMesh","3DCoRectMesh");
  if (! (gridtype == ok_grid_type) ) {
    std::string err("color function: invalid GRIDTYPE");
    err += " " + gridtype;
    BoxLib::Abort(err.c_str());
  }

  Array<int> length(BL_SPACEDIM);
  for (int k = 0; k < length.size(); ++k) {
    is >> length[k];
    if (!is.good()) {
      BoxLib::Abort("color function: error reading NXNYNZ record");
    }
    if (length[k] < 1) {
      BoxLib::Abort("color function: invalid NXNYNZ record");
    }
  }

  for (int k = 0; k < BL_SPACEDIM; ++k) {
    is >> vertex_lo[k];
    if (!is.good()) {
      BoxLib::Abort("color function: error reading CORNERLO record");
    }
  }
  for (int k = 0; k < BL_SPACEDIM; ++k) {
    is >> vertex_hi[k];
    if (!is.good()) {
      BoxLib::Abort("color function: error reading CORNERHI record");
    }
  }
  for (int k = 0; k < BL_SPACEDIM; ++k) {
    dx[k] = (vertex_hi[k] - vertex_lo[k]) / length[k];
    if (dx[k]<=0) {
        std::cout << "color function: CORNERLO  (" << vertex_lo[k]
                  <<") > CORNERHI (" << vertex_hi[k] << ")";
        BoxLib::Abort();
    }
  }

  int dataloc;
  is >> dataloc;
  IndexType iType;
  if (!is.good()) {
    BoxLib::Abort("color function: error reading DATALOC record");
  }
  if (dataloc != 0) {
    BoxLib::Abort("color function: invalid DATALOC value");
  } else {
    // Cell centered in all coordinates
    iType = IndexType(IntVect::TheUnitVector()*dataloc);
  }

  int ncomp;
  is >> ncomp;
  if (!is.good()) {
    BoxLib::Abort("color function: error reading DATACOL record");
  }
  if (ncomp < 1) {
    BoxLib::Abort("color function: invalid DATACOL value");
  }
  if (ncomp > 1) {
    BoxLib::Abort("color function: supports only a single component file");
  }

  IntVect be = IntVect(length.dataPtr()) - IntVect::TheUnitVector();
  Box map_box(IntVect::TheZeroVector(),be,iType);  
  m_color_map = new BaseFab<int>(map_box,ncomp);

  if (!map_box.numPtsOK()) {
    BoxLib::Abort("color function: input box too big");
  }
  long Npts = map_box.numPts();
  Array<int*> dptr(ncomp);
  for (int n=0; n<ncomp; ++n) {
    dptr[n] = m_color_map->dataPtr();
  }

  for (int i = 0; i < Npts; ++i) {
    for (int n=0; n<ncomp; ++n) {
      is >> dptr[n][i];
      if (!is.good()) {
	BoxLib::Abort("color function: error reading DATAVAL records");
      }
    }
  }
}

IntVect
colorFunctionRegion::atIndex(Array<Real> x) const
{
    IntVect idx;
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
      idx[d] = (int)( (x[d] - vertex_lo[d])/dx[d] - 0.5);
    }
    return idx;
}


bool
colorFunctionRegion::inregion (Array<Real>& x)
{
    if (! boxRegion::inregion(x)) {
        return false;
    }
    IntVect idx = atIndex(x);
    if ( !(m_color_map->box().contains(idx)) ) { // in hacked up boundary region
        // Push index into "domain" and check there instead
        const Box& domain = m_color_map->box();
        for (int d=0; d<BL_SPACEDIM; ++d)
        {
            idx[d] = std::max(domain.smallEnd()[d],std::min(domain.bigEnd()[d],idx[d]));
        }
    }
    BL_ASSERT(m_color_map->box().contains(idx));
    return (*m_color_map)(idx,0) == m_color_val;
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
  : Region("BC","bc","bc")
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
