/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <winstd.H>
#include "iostream"

#include "Region.H"

Real Region::geometry_eps = -1;
Array<Real> Region::domlo;
Array<Real> Region::domhi;

Region::Region (const std::string& r_name,
		const std::string& r_purpose,
		const std::string& r_type)
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
Region::setVal(BaseFab<Real>&     fab,
               const Array<Real>& val,
               const Real*        dx,
               int                ng,
               int                scomp,
               int                ncomp) const
{
  setVal_DoIt(fab,val,dx,ng,scomp,ncomp);
}

void
Region::setVal(BaseFab<int>&     fab,
               const Array<int>& val,
               const Real*       dx,
               int               ng,
               int               scomp,
               int               ncomp) const
{
  setVal_DoIt(fab,val,dx,ng,scomp,ncomp);
}

template <class T>
void
Region::setVal_DoIt(BaseFab<T>&     fab,
                    const Array<T>& val,
                    const Real*     dx,
                    int             ng,
                    int             scomp,
                    int             ncomp) const
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
      if (inRegion(x))
      {
          for (int n=scomp; n<scomp+ncomp;n++) {
              fab(idx,n) = val[n-scomp];
          }
      }
  }
}

void
Region::setVal(BaseFab<Real>& fab,
               const Real&    val,
               const int      idx,
               const Real*    dx,
               int            ng) const
{
    setVal(fab,Array<Real>(1,val),dx,ng,idx,1);
}

void
Region::setVal(BaseFab<int>& fab,
               const int&    val,
               const int     idx,
               const Real*   dx,
               int           ng) const
{
    setVal(fab,Array<int>(1,val),dx,ng,idx,1);
}

BoxArray
PointRegion::approximate_bounds(const Array<Real>& plo, const Array<Real>& dx) const
{
  IntVect iv;
  for (int d=0; d<BL_SPACEDIM; ++d) {
    iv[d] = (int)((coor[d]-plo[d])/dx[d]);
  }
  return BoxArray(Box(iv,iv));
}

bool
PointRegion::inRegion (const Array<Real>& x) const
{
  return false; // This test is simply not appropriate for pt regions
}

void
PointRegion::setVal(BaseFab<Real>&     fab,
                    const Array<Real>& val,
                    const Real*        dx,
                    int                ng,
                    int                scomp,
                    int                ncomp) const
{
  setVal_DoIt(fab,val,dx,ng,scomp,ncomp);
}

void
PointRegion::setVal(BaseFab<int>&     fab,
                    const Array<int>& val,
                    const Real*       dx,
                    int               ng,
                    int               scomp,
                    int               ncomp) const
{
  setVal_DoIt(fab,val,dx,ng,scomp,ncomp);
}

template<class T>
void
PointRegion::setVal_DoIt(BaseFab<T>&     fab,
                         const Array<T>& val,
                         const Real*     dx,
                         int             ng,
                         int             scomp,
                         int             ncomp) const
{
  int nval = val.size();
  BL_ASSERT(ncomp-scomp <= nval);

  const Box& box = fab.box();
  const IntVect& se = box.smallEnd();
  const IntVect& be = box.bigEnd();
  Real xlo[BL_SPACEDIM] = {D_DECL(domlo[0]  +  se[0]*dx[0],domlo[1]  +  se[1]*dx[1],domlo[2]  +  se[2]*dx[2])};
  Real xhi[BL_SPACEDIM] = {D_DECL(domlo[0]+(be[0]+1)*dx[0],domlo[1]+(be[1]+1)*dx[1],domhi[2]+(be[2]+1)*dx[2])};

  bool not_below = D_TERM(coor[0]>=xlo[0], && coor[1]>=xlo[1], && coor[2]>=xlo[2]);
  bool not_above = D_TERM(coor[0]<=xhi[0], && coor[1]<=xhi[1], && coor[2]<=xhi[2]);

  bool in_domain = not_below && not_above;
  if (in_domain) {
    IntVect idx(D_DECL((coor[0]-domlo[0])/dx[0],(coor[1]-domlo[1])/dx[1],(coor[2]-domlo[2])/dx[2]));
    if (fab.box().contains(idx)) {
      for (int n=scomp; n<scomp+ncomp;n++) {
        fab(idx,n) = val[n-scomp];
      }
    }
  }
}


void
PointRegion::setVal(BaseFab<Real>& fab,
                    const Real&    val,
                    const int      idx,
                    const Real*    dx,
                    int            ng) const
{
  setVal(fab,Array<Real>(1,val),dx,ng,idx,1);
}

void
PointRegion::setVal(BaseFab<int>& fab,
                    const int&    val,
                    const int     idx,
                    const Real*   dx,
                    int           ng) const
{
  setVal(fab,Array<int>(1,val),dx,ng,idx,1);
}

BoxArray
BoxRegion::approximate_bounds(const Array<Real>& plo, const Array<Real>& dx) const
{
  IntVect ivlo, ivhi;
  for (int d=0; d<BL_SPACEDIM; ++d) {
    ivlo[d] = (int)((lo[d]-plo[d]-geometry_eps)/dx[d]);
    ivhi[d] = (int)((hi[d]-plo[d]+geometry_eps)/dx[d]);
  }
  return BoxArray(Box(ivlo,ivhi));
}

bool
BoxRegion::inRegion(const Array<Real>& x) const
{
    bool inflag = true;
    for (int d=0; d<BL_SPACEDIM; ++d) {
        inflag &= (x[d]>=lo[d] && x[d]<=hi[d]);
    }
    return inflag;
}

static bool
pnpoly(const Array<Real>& vertx,
       const Array<Real>& verty,
       const Array<Real>& x)
{
  int i, j;
  bool ret = false;
  int nvert = vertx.size();
  BL_ASSERT(verty.size()==nvert);
  BL_ASSERT(x.size()==2);
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>x[1]) != (verty[j]>x[1])) &&
	 (x[0] < (vertx[j]-vertx[i]) * (x[1]-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ) {
      ret = !ret;
    }
  }
  return ret;
}

#if BL_SPACEDIM==2
PolygonRegion::PolygonRegion (const std::string& r_name,
                              const std::string& r_purpose,
                              const Array<Real>& v1vector,
                              const Array<Real>& v2vector)
  : Region(r_name,r_purpose,"polygon")
{
  BL_ASSERT(v1vector.size()==v2vector.size());
  int N = v1vector.size();
  vector1.resize(N);
  vector2.resize(N);
  for (int i=0; i<N; ++i) {
    vector1[i] = v1vector[i];
    vector2[i] = v2vector[i];
  }
}


bool
PolygonRegion::inRegion (const Array<Real>& x) const
{
  return pnpoly(vector1,vector2,x);
}

BoxArray
PolygonRegion::approximate_bounds(const Array<Real>& plo, const Array<Real>& dx) const
{
  IntVect ivlo, ivhi;
  int N = vector1.size();
  ivlo[0] = (int)((vector1[0]-geometry_eps-plo[0])/dx[0]);
  ivlo[1] = (int)((vector2[0]-geometry_eps-plo[1])/dx[1]);
  ivhi[0] = (int)((vector1[0]+geometry_eps-plo[0])/dx[0]);
  ivhi[1] = (int)((vector2[0]+geometry_eps-plo[1])/dx[1]);

  for (int i=1; i<N; ++i) {
    ivlo[0] = std::min(ivlo[0], (int)((vector1[i]-geometry_eps-plo[0])/dx[0]) );
    ivlo[1] = std::min(ivlo[1], (int)((vector2[i]-geometry_eps-plo[1])/dx[1]) );

    ivhi[0] = std::max(ivhi[0], (int)((vector1[i]+geometry_eps-plo[0])/dx[0]) );
    ivhi[1] = std::max(ivhi[1], (int)((vector2[i]+geometry_eps-plo[1])/dx[1]) );
  }
  return BoxArray(Box(ivlo,ivhi));
}

std::ostream&
PolygonRegion::print(std::ostream& os) const
{
  Region::print(os);
  os << "    polygon nodes ";
  for (int i=0; i<vector1.size(); ++i) {
    os << "(" << vector1[i] << ", " << vector2[i] << ") ";
  }
  os << '\n';
  return os;
}

EllipseRegion::EllipseRegion (const std::string& r_name,
                              const std::string& r_purpose,
                              const Array<Real>& r_center,
                              const Array<Real>& r_radius)
  : Region(r_name,r_purpose,"ellipse")
{
  BL_ASSERT(r_center.size()==BL_SPACEDIM);
  BL_ASSERT(r_radius.size()==BL_SPACEDIM);
  center.resize(BL_SPACEDIM);
  radius.resize(BL_SPACEDIM);
  for (int i=0; i<BL_SPACEDIM; ++i) {
    center[i] = r_center[i];
    radius[i] = r_radius[i];
  }
}


bool
EllipseRegion::inRegion (const Array<Real>& x) const
{
  Real this_rsqrd = 0;
  for (int i=0; i<BL_SPACEDIM; ++i) {
    Real d = (x[i] - center[i])/radius[i];
    this_rsqrd += d*d;
  }
  return this_rsqrd <= 1;
}

BoxArray
EllipseRegion::approximate_bounds(const Array<Real>& plo, const Array<Real>& dx) const
{
  IntVect ivlo, ivhi;
  for (int d=0; d<BL_SPACEDIM; ++d) {
    ivlo[d] = (int)((center[d]-radius[d]-geometry_eps-plo[d])/dx[d]);
    ivhi[d] = (int)((center[d]+radius[d]+geometry_eps-plo[d])/dx[d]);
  }
  return BoxArray(Box(ivlo,ivhi));
}

std::ostream&
EllipseRegion::print(std::ostream& os) const
{
  Region::print(os);
  os << "    ellipse radii: (";
  for (int i=0; i<radius.size(); ++i) {
    os << radius[i] << " ";
  }
  os << ") , center: (";
  for (int i=0; i<center.size(); ++i) {
    os << center[i] << " ";
  }
  os << ")\n";
  return os;
}

#else

SweptPolygonRegion::SweptPolygonRegion (const std::string& r_name,
                                        const std::string& r_purpose,
                                        const Array<Real>& v1vector,
                                        const Array<Real>& v2vector,
                                        const PLANE&       _plane,
                                        const Array<Real>& _extent)
  : Region(r_name,r_purpose,"polygon")
{
  BL_ASSERT(v1vector.size()==v2vector.size());
  int N = v1vector.size();
  vector1.resize(N);
  vector2.resize(N);
  for (int i=0; i<N; ++i) {
    vector1[i] = v1vector[i];
    vector2[i] = v2vector[i];
  }
  plane = _plane;
  extent.resize(2);
  for (int i=0; i<2; ++i) {
    extent[i] = _extent[i];
  }
}


bool
SweptPolygonRegion::inRegion (const Array<Real>& x) const
{
  Array<Real> xproj(2);
  bool ret = true;
  if (plane == XYPLANE) {
    xproj[0] = x[0];
    xproj[1] = x[1];
    ret &= x[2] >= extent[0] & x[2] < extent[1];
  }
  else if (plane == YZPLANE) {
    xproj[0] = x[1];
    xproj[1] = x[2];
    ret &= x[0] >= extent[0] & x[0] < extent[1];
  }
  else {
    xproj[0] = x[0];
    xproj[1] = x[2];
    ret &= x[1] >= extent[0] & x[1] < extent[1];
  }
  if (!ret) {
    return ret;
  }
  return pnpoly(vector1,vector2,xproj);
}

BoxArray
SweptPolygonRegion::approximate_bounds(const Array<Real>& plo,
                                       const Array<Real>& dx) const
{
  return BoxArray(0);
}

std::ostream&
SweptPolygonRegion::print(std::ostream& os) const
{
  Region::print(os);
  os << "    polygon nodes ";
  std::string plStr = (plane==XYPLANE ? "XY" : (plane == YZPLANE ? "YZ" : "XZ"));
  os << "in the " << plStr << " plane: ";
  for (int i=0; i<vector1.size(); ++i) {
    os << "(" << vector1[i] << ", " << vector2[i] << ") ";
  }
  os << '\n';
  os << "    extent: (" << extent[0] << ", " << extent[1] << ")";
  os << '\n';
  return os;
}


RotatedPolygonRegion::RotatedPolygonRegion (const std::string& r_name,
                                            const std::string& r_purpose,
                                            const Array<Real>& v1vector,
                                            const Array<Real>& v2vector,
                                            const PLANE&       _plane,
                                            const Array<Real>& _reference_pt,
                                            const std::string& _axis)
  : Region(r_name,r_purpose,"rotated_polygon")
{
  BL_ASSERT(v1vector.size()==v2vector.size());
  int N = v1vector.size();
  vector1.resize(N);
  vector2.resize(N);
  for (int i=0; i<N; ++i) {
    vector1[i] = v1vector[i];
    vector2[i] = v2vector[i];
  }
  plane = _plane;
  reference_pt.resize(BL_SPACEDIM);
  for (int i=0; i<reference_pt.size(); ++i) {
    reference_pt[i] = _reference_pt[i];
  }
  axis = _axis;
  BL_ASSERT(axis=="X" || axis=="Y" || axis=="Z");
  n = (axis=="X" ?  0 :  (axis=="Y" ? 1 : 2));
  if (plane==XYPLANE) {
    t1 = (n==0 ? 1 : 0);
    t2 = 2;
  }
  else if (plane==YZPLANE) {
    t1 = 0;
    t2 = (n==1 ? 2 : 1);
  }
  else  {
    BL_ASSERT(plane==XZPLANE);
    t1 = 1;
    t2 = (n==2 ? 0 : 2);
  }
}


bool
RotatedPolygonRegion::inRegion (const Array<Real>& x) const
{
  Array<Real> xproj(2);
  Real d1 = x[t1]-reference_pt[t1];
  Real d2 = x[t2]-reference_pt[t2];
  Real d3 = x[n]-reference_pt[n];

  xproj[0] = std::sqrt( d1*d1 + d2*d2 );
  xproj[1] = d3;

  return pnpoly(vector1,vector2,xproj);
}

BoxArray
RotatedPolygonRegion::approximate_bounds(const Array<Real>& plo,
                                         const Array<Real>& dx) const
{
  return BoxArray(0);
}

std::ostream&
RotatedPolygonRegion::print(std::ostream& os) const
{
  Region::print(os);
  os << "    polygon nodes ";
  std::string plStr = (plane==XYPLANE ? "XY" : (plane == YZPLANE ? "YZ" : "XZ"));
  os << "in the " << plStr << " plane: ";
  for (int i=0; i<vector1.size(); ++i) {
    os << "(" << vector1[i] << ", " << vector2[i] << ") ";
  }
  os << '\n';
  return os;
}

#endif



ColorFunctionRegion::ColorFunctionRegion (const std::string& r_name,
                                          const std::string& r_purpose,
					  const std::string& file_name,
                                          int                color_val)
    : BoxRegion(r_name,r_purpose,
                Array<Real>(BL_SPACEDIM),Array<Real>(BL_SPACEDIM)),
    dx(BL_SPACEDIM),
    m_color_val(color_val),
    m_file(file_name)
{
  type = "color_function";
  set_color_map();
}

void
ColorFunctionRegion::set_color_map()
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
ColorFunctionRegion::atIndex(const Array<Real> x) const
{
    IntVect idx;
    for (int d=0; d<BL_SPACEDIM; ++d)
    {
        idx[d] = (int) ((x[d] - lo[d])/dx[d]);
    }
    return idx;
}


bool
ColorFunctionRegion::inRegion (const Array<Real>& x) const
{
    if (! BoxRegion::inRegion(x)) {
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

BoxArray
ColorFunctionRegion::approximate_bounds(const Array<Real>& plo, const Array<Real>& dx) const
{
  return BoxRegion::approximate_bounds(plo,dx);
}

std::ostream& operator<< (std::ostream& os, const Region& rhs)
{
  return rhs.print(os);
}

std::ostream&
Region::print (std::ostream& os) const
{
    os << "Region:\n";
    os << "  name:    " << name << '\n';
    os << "  purpose: " << purpose << '\n';
    os << "  type:    " << type << '\n';
    return os;
}

std::ostream&
PointRegion::print (std::ostream& os) const
{
    Region::print(os);
    os << "    coor: ";
    for (int i=0; i<coor.size(); ++i) {
        os << coor[i] << " ";
    }
    os << '\n';
    return os;
}

std::ostream&
BoxRegion::print (std::ostream& os) const
{
    Region::print(os);
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
    return os;
}

std::ostream&
ColorFunctionRegion::print (std::ostream& os) const
{
    BoxRegion::print(os);
    os << "    color_map from file: " << m_file << '\n';
    os << "    FAB data on " << m_color_map->box() << '\n';
    os << "    dx = ";
    for (int i=0; i<dx.size(); ++i) {
        os << dx[i] << " ";
    }
    os << " color val = " << m_color_val << '\n';
    return os;
}

#if 1
bool
AllRegion::inRegion(const Array<Real>& x) const
{
  return true;
}

BoxArray
AllRegion::approximate_bounds(const Array<Real>& plo,
                              const Array<Real>& dx) const
{
  return BoxArray(0);
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


AllBCRegion::AllBCRegion (int dir, int lo_or_hi)
    : BoxRegion(pick_name(dir,lo_or_hi),pick_purpose(dir,lo_or_hi),domlo,domhi)
{
  p_dir  = dir;
  p_lohi = lo_or_hi;

  if (p_lohi == 0)
      hi[p_dir] = lo[p_dir];
  else
      lo[p_dir] = hi[p_dir];
}

bool
AllBCRegion::inRegion(const Array<Real>& x) const
{
    return true;
}
#endif

ComplementRegion::ComplementRegion (const std::string& r_name,
                                    const std::string& r_purpose,
                                    Region*            _exclude_region)
  : Region(r_name,r_purpose,"complement")
{
  exclude_region = _exclude_region;
}

bool
ComplementRegion::inRegion (const Array<Real>& x) const
{
  return (exclude_region==0 ? true : !exclude_region->inRegion(x));
}

BoxArray
ComplementRegion::approximate_bounds(const Array<Real>& plo,
                                     const Array<Real>& dx) const
{
  return BoxArray(0); // Have not yet optimized this
}



std::ostream&
ComplementRegion::print (std::ostream& os) const
{
  Region::print(os);
  os << "    contains entire domain except the region: " << exclude_region->name << '\n';
  return os;
}


UnionRegion::UnionRegion (const std::string&    r_name,
                          const std::string&    r_purpose,
                          const Array<const Region*>& _regions)
  : Region(r_name,r_purpose,"union")
{
  SetUnionRegions(_regions);
}

void
UnionRegion::SetUnionRegions(const Array<const Region*>& _regions)
{
  regions.resize(_regions.size());
  for (int i=0; i<_regions.size(); ++i) {
    regions[i] = _regions[i];
  }
}

bool
UnionRegion::inRegion (const Array<Real>& x) const
{
  if (regions.size()==0) {
    return false;
  }

  bool ret = false;
  for(int i=0; !ret && i<regions.size(); ++i) {
    ret |= regions[i]->inRegion(x);
  }
  return ret;
}

BoxArray
UnionRegion::approximate_bounds(const Array<Real>& plo,
                                const Array<Real>& dx) const
{
  BoxArray ba = regions[0]->approximate_bounds(plo,dx);
  if (ba.ok()) {
    BoxList bl(ba);
    for(int i=1; i<regions.size(); ++i) {
      ba = regions[i]->approximate_bounds(plo,dx);
      if (ba.ok()) {
        bl.join(BoxList(ba));
      }
    }
    ba = BoxArray(bl);
    ba.removeOverlap();
  }
  return ba;
}

std::ostream&
UnionRegion::print (std::ostream& os) const
{
  os << "The union of the following regions: \n";
  for(int i=0; i<regions.size(); ++i) {
    regions[i]->print(os);
  }
  return os;
}



SubtractionRegion::SubtractionRegion (const std::string&    r_name,
                                      const std::string&    r_purpose,
                                      const Array<const Region*>& _regions)
  : Region(r_name,r_purpose,"union")
{
  SetRegions(_regions);
}

void
SubtractionRegion::SetRegions(const Array<const Region*>& _regions)
{
  regions.resize(_regions.size());
  for (int i=0; i<_regions.size(); ++i) {
    regions[i] = _regions[i];
  }
}

bool
SubtractionRegion::inRegion (const Array<Real>& x) const
{
  if (regions.size() == 0) return true; // Default master region to All, remove no other regions
  if (! (regions[0]->inRegion(x)) ) return false;

  for(int i=1; i<regions.size(); ++i) {
    if (regions[i]->inRegion(x)) {
      return false;
    }
  }
  return true;
}

BoxArray
SubtractionRegion::approximate_bounds(const Array<Real>& plo,
                                      const Array<Real>& dx) const
{
  if (regions.size() == 0) return BoxArray(0); // Default master region to All, remove no other regions
  return regions[0]->approximate_bounds(plo,dx); // HAve not yet optimized to remove subtracted pieces
}

std::ostream&
SubtractionRegion::print (std::ostream& os) const
{
  if (regions.size()>0) {
    os << "Within the region: \n";
    regions[0]->print(os);
    if (regions.size()>1) {
      os << "  but outside the region: \n";
      for(int i=1; i<regions.size(); ++i) {
        regions[i]->print(os);
      }
    }
  }
  return os;
}
