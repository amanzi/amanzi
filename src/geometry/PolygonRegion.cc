/**
 * @file   PolygonRegion.cc
 * @author Rao Garimella
 * @date Fri Jul 29 12:28:10 2011
 * 
 * @brief  Implementation of PolygonRegion class 
 * 
 * 
 */

#include "PolygonRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class PolygonRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// Polygon:: constructors / destructor
// -------------------------------------------------------------
PolygonRegion::PolygonRegion(const std::string name, const unsigned int id,
                             const unsigned int num_points, 
                             const std::vector<Point>& points,
                             const LifeCycleType lifecycle)
  : Region(name,id,points[0].dim(),lifecycle), num_points_(num_points), 
    points_(points),normal_(points[0].dim()),elim_dir_(0)
{
  init();
}

PolygonRegion::PolygonRegion(const char *name, const unsigned int id,
                 const unsigned int num_points,
                 const std::vector<Point>& points,
                 const LifeCycleType lifecycle)
  : Region(name,id,points[0].dim(),lifecycle), num_points_(num_points), 
    points_(points),normal_(points[0].dim()),elim_dir_(0)
{
  init();
}

PolygonRegion::PolygonRegion(const PolygonRegion& old)
  : Region(old), num_points_(old.num_points_), points_(old.points_),
  normal_(old.normal_), elim_dir_(old.elim_dir_)
{
  // empty
}

PolygonRegion::~PolygonRegion(void)
{
  
}

void PolygonRegion::init() {

  if (num_points_ < dimension()) {
    std::stringstream tempstr;
    tempstr << "\nDimension " << dimension() << 
      " regions need to be specified by at least " << dimension() << 
      " points\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }
#ifdef ENABLE_DBC
  if (dimension() == 2 && num_points_ > 2) {
    std::cerr << "\nDimension " << dimension() << 
      " regions specified by more points (" << num_points_ << ") than needed\n";
    std::cerr << "Using only the first two\n";
  }
#endif
  
  if (dimension() == 2) {
    Point vec = points_[1] - points_[0];
    vec /= norm(vec);
    normal_.set(vec[1],-vec[0]);

    elim_dir_ = (vec[0] < vec[1]) ? 0 : 1;
  }
  else if (dimension() == 3) {
    Point vec0 = points_[2]-points_[1];
    Point vec1 = points_[0]-points_[1];
    normal_ = vec0^vec1;
    normal_ /= norm(normal_);

#ifdef ENABLE_DBC    
    for (int i = 3; i < num_points_; i++) {
      vec0 = points_[(i+1)%num_points_]-points_[i];
      vec1 = points_[(i-1+num_points_)%num_points_]-points_[i];
      Point nrml = vec0^vec1;
      nrml /= norm(nrml);

      double dp = nrml*normal_;
      if (fabs(dp-1.0) > 1.0e-06) {
        Errors::Message mesg("Polygon region is not exactly planar");
        Exceptions::amanzi_throw(mesg);
      }
    }
#endif

    /* Determine which direction to eliminate while projecting to one
       of the coordinate planes - to do this we have to find the
       direction in which the polygon is the smallest or in other
       words, the direction in which the normal to the polygon is the
       largest */

    int dmax = -1; 
    double maxlen = -1;
    for (int i = 0; i < 3; i++)
      if (normal_[i] > maxlen) {
        maxlen = normal_[i];
        dmax = i;
      }
       
    elim_dir_ = dmax;
  }

#ifdef ENABLE_DBC
  else {
    std::stringstream tempstr;
    tempstr << "Cannot handle regions of dimension " << dimension() << "\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }
#endif
}

// -------------------------------------------------------------
// PolygonRegion::inside
// -------------------------------------------------------------
bool
PolygonRegion::inside(const Point& p) const
{

#ifdef ENABLE_DBC
  if (p.dim() != points_[0].dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in corner dimension of Polygon \"" << Region::name() << "\" and query point.\n Perhaps the region is improperly defined?\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }
#endif

  /* First check if the point is on the infinite line/plane */

  double d(0.0), res(0.0);

  for (int i = 0; i < p.dim(); ++i) {
    res += normal_[i]*p[i];
    d += normal_[i]*points_[0][i];
  }
  res -= d;

  if (fabs(res) > 1.0e-12)
    return false;


  bool result(false);
  if (dimension() == 2) {
    /* Now check if it lies in the line segment */
    Point vec0 = points_[0]-p;
    Point vec1 = points_[1]-p;
    if (vec0*vec1 < 0.0)
      result = true;
  }
  else {
    /* Now check if the point is in the polygon */

    /* 
       Find the indices of coordinates on the projection plane
       
       if elim_dir_ is 0, then d0 = 1, d1 = 2 (YZ plane)
       if elim_dir_ is 1, then d0 = 2, d1 = 0 (XZ plane)
       if elim_dir_ is 2, then d0 = 0, d1 = 1 (XY plane)
    */
    
    double d0 = (elim_dir_+1)%3;
    double d1 = (elim_dir_+2)%3;
    
    /* Now apply the Jordan curve theorem to do the in/out test */
    /* odd number of crossings - point is inside                */
    
    double u, v;
    u = p[d0]; v = p[d1];
    
    for (int i = 0; i < num_points_; i++) {
      int iplus1 = (i+1)%num_points_;
      double u_i = points_[i][d0];
      double v_i = points_[i][d1];
      double u_iplus1 = points_[iplus1][d0];
      double v_iplus1 = points_[iplus1][d1];
      
      double slope = (u_iplus1-u_i)/(v_iplus1-v_i);
      
      if (((v_i > v && v_iplus1 <= v) || (v_iplus1 > v && v_i <= v)) &&
          (u <= (u_i + (v-v_i)*slope)))
        result = !result;
    }
  }
  return result;
}


} // namespace AmanziGeometry
} // namespace Amanzi
