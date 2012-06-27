#include <winstd.H>
#include "iostream"

#include "Region.H"

Real Region::geometry_eps = -1;
Array<Real> Region::domlo;
Array<Real> Region::domhi;

Region::Region (std::string r_name, 
		std::string r_purpose, 
		std::string r_type)
  : name(r_name),
    purpose(r_purpose),
    type(r_type)
{
    if (geometry_eps<0) {
        BoxLib::Abort("Static variable Region::geometry_eps must be set before Region ctr");
    }
    if (domlo.size() == 0) {
        BoxLib::Abort("Static variable Region::domlo must be set before Region ctr");
    }
    if (domhi.size() == 0) {
        BoxLib::Abort("Static variable Region::domhi must be set before Region ctr");
    }
}

void
Region::setVal(FArrayBox&         fab, 
               const Array<Real>& val, 
               const Real*        dx, 
               int                ng,
               int                scomp,
               int                ncomp) const
{
  int nval = val.size();
  BL_ASSERT(ncomp-scomp <= nval);

  Array<Real> x(BL_SPACEDIM);
  const Box& box = fab.box();
  const IntVect& se = box.smallEnd();
  const IntVect& be = box.bigEnd();
  for (IntVect idx=se; idx<=be; box.next(idx)) {
      for (int d=0; d<BL_SPACEDIM; ++d) {
          x[d] = std::min(std::max(domlo[d],domlo[d]+dx[d]*(idx[d]+0.5)),domhi[d]);
      }
      if (inregion(x))
      {
          for (int n=scomp; n<scomp+ncomp;n++) {
              fab(idx,n) = val[n-scomp];
          }
      }
  }
}


void
Region::setVal(FArrayBox&  fab,
               const Real  val, 
               const int   idx,
               const Real* dx, 
               int         ng) const
{
    setVal(fab,Array<Real>(1,val),dx,ng,idx,1);
}

bool 
pointRegion::inregion (const Array<Real>& x) const
{
    bool inflag = true;
    for (int d=0; d<BL_SPACEDIM; ++d) {
        inflag &= (x[d]>=coor[d]);
    }
    return inflag;
}

bool 
boxRegion::inregion(const Array<Real>& x) const
{
    bool inflag = true;
    for (int d=0; d<BL_SPACEDIM; ++d) {
        inflag &= (x[d]>=lo[d] && x[d]<=hi[d]);
    }
    return inflag;
}

colorFunctionRegion::colorFunctionRegion (std::string r_name,
                                          std::string r_purpose,
                                          std::string r_type,
					  std::string file_name,
                                          int         color_val)
    : boxRegion(r_name,r_purpose,r_type,
                Array<Real>(BL_SPACEDIM),Array<Real>(BL_SPACEDIM)), 
    dx(BL_SPACEDIM),
    m_color_val(color_val),
    m_file(file_name)
{
    set_color_map();
}

void
colorFunctionRegion::set_color_map()
{
  Array<char> fileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(m_file, fileCharPtr);
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
    is >> lo[k];
    if (!is.good()) {
      BoxLib::Abort("color function: error reading CORNERLO record");
    }
  }
  for (int k = 0; k < BL_SPACEDIM; ++k) {
    is >> hi[k];
    if (!is.good()) {
      BoxLib::Abort("color function: error reading CORNERHI record");
    }
  }
  for (int k = 0; k < BL_SPACEDIM; ++k) {
    dx[k] = (hi[k] - lo[k]) / length[k];
    if (dx[k]<=0) {
        std::cout << "color function: CORNERLO  (" << lo[k]
                  <<") > CORNERHI (" << hi[k] << ")";
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
colorFunctionRegion::atIndex(const Array<Real> x) const
{
    IntVect idx;
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
        idx[d] = (int) ((x[d] - lo[d])/dx[d]);
    }
    return idx;
}


bool
colorFunctionRegion::inregion (const Array<Real>& x) const
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

std::ostream& operator<< (std::ostream& os, const Region& rhs)
{
    rhs.operator<<(os);
}

std::ostream&
Region::operator<< (std::ostream& os) const
{
    os << "Region:\n";
    os << "  name:    " << name << '\n';
    os << "  purpose: " << purpose << '\n';
    os << "  type:    " << type << '\n';
}

std::ostream&
pointRegion::operator<< (std::ostream& os) const
{
    Region::operator<<(os);
    os << "    coor: ";
    for (int i=0; i<coor.size(); ++i) {
        os << coor[i] << " ";
    }
    os << '\n';
}

std::ostream&
boxRegion::operator<< (std::ostream& os) const
{
    Region::operator<<(os);
    os << "    lo: ";
    for (int i=0; i<lo.size(); ++i) {
        os << lo[i] << " ";
    }
    os << '\n';
    os << "    hi: ";
    for (int i=0; i<hi.size(); ++i) {
        os << hi[i] << " ";
    }
    os << '\n';
}

std::ostream&
colorFunctionRegion::operator<< (std::ostream& os) const
{
    boxRegion::operator<<(os);
    os << "    color_map from file: " << m_file << '\n';
    os << "    FAB data on " << m_color_map->box() << '\n';
    os << "    dx = ";
    for (int i=0; i<dx.size(); ++i) {
        os << dx[i] << " ";
    }
    os << " color val = " << m_color_val << '\n';
}

#if 1
bool 
allRegion::inregion(const Array<Real>& x) const
{
  return true;
}

static
std::string
pick_name(int dir, int lo_or_hi)
{
    std::string name;
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
    return name;
}

static
std::string
pick_purpose(int dir, int lo_or_hi)
{
    std::string purpose;
    if (dir == 0 && lo_or_hi == 0)
        purpose = "xlobc";
    else if (dir == 0 && lo_or_hi == 1)
        purpose = "xhibc";
    else if (dir == 1 && lo_or_hi == 0)
        purpose = "ylobc";
    else if (dir == 1 && lo_or_hi == 1)
        purpose = "yhibc";
#if BL_SPACEDIM == 3
    else if (dir == 2 && lo_or_hi == 0)
        purpose = "zlobc";
    else if (dir == 2 && lo_or_hi == 1)
        purpose = "zhibc";
#endif  
    return purpose;
}


allBCRegion::allBCRegion (int dir, int lo_or_hi,
                          const Array<Real>& lo_,
                          const Array<Real>& hi_)
    : boxRegion(pick_name(dir,lo_or_hi),pick_purpose(dir,lo_or_hi),"box",lo_,hi_)
{
  p_dir  = dir;
  p_lohi = lo_or_hi;

  if (p_lohi == 0)
      hi[p_dir] = lo[p_dir];
  else 
      lo[p_dir] = hi[p_dir];
}

bool 
allBCRegion::inregion(const Array<Real>& x) const
{
    return true;
}
#endif

