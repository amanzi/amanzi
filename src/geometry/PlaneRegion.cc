/**
 * @file   PlaneRegion.cc
 * @author Rao Garimella
 * @date 
 * 
 * @brief  Implementation of PlaneRegion class
 * 
 * 
 */

#include "PlaneRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class PlaneRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// PlaneRegion:: constructors / destructor
// -------------------------------------------------------------
PlaneRegion::PlaneRegion(const std::string name, 
			 const unsigned int id,
			 const Point& p, const Point& normal)
  : Region(name,id,p.dim()-1), p_(p), n_(normal)
{

#ifdef ENABLE_DBC
  if (p_.dim() != n_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in point and normal dimensions of PlaneRegion " << Region::name() << "Perhaps the region is improperly defined?\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }
#endif

}

PlaneRegion::PlaneRegion(const char *name, const unsigned int id,
			 const Point& p, const Point& normal)
  : Region(name,id,p.dim()-1), p_(p), n_(normal)
{

#ifdef ENABLE_DBC
  if (p_.dim() != n_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in point and normal dimensions of PlaneRegion " << Region::name() << "Perhaps the region is improperly defined?\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }
#endif

}

PlaneRegion::PlaneRegion(const PlaneRegion& old)
  : Region(old), p_(old.p_), n_(old.n_)
{
  // empty
}

PlaneRegion::~PlaneRegion(void)
{
  
}

// -------------------------------------------------------------
// PlaneRegion::inside -- check if point is on plane
// -------------------------------------------------------------
bool
PlaneRegion::inside(const Point& p) const
{
#ifdef ENABLE_DBC
  if (p.dim() != p_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in point dimension of PlaneRegion \"" << Region::name() << "\" and query point.\n Perhaps the region is improperly defined?\n";
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }
#endif

  bool result(true);

  double d(0.0), res(0.0);

  for (int i = 0; i < p.dim(); ++i) 
    {
      res += n_[i]*p[i];
      d += n_[i]*p_[i];
    }
  res -= d;

  if (fabs(res) <= 1.0e-12)
    result = true;
  else
    result = false;

  return result;
}

} // namespace AmanziGeometry
} // namespace Amanzi
