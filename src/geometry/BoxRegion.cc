/**
 * @file   BoxRegion.cc
 * @author Rao Garimella, William A. Perkins
 * @date Fri Jul 29 12:28:10 2011
 * 
 * @brief  Implementation of BoxRegion class (Adapted from RectangularRegion.cc)
 * 
 * 
 */

#include "BoxRegion.hh"
#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
//  class BoxRegion
// -------------------------------------------------------------

// -------------------------------------------------------------
// BoxRegion:: constructors / destructor
// -------------------------------------------------------------
BoxRegion::BoxRegion(const Set_Name& name, const Set_ID id,
                     const Point& p0, const Point& p1,
                     const LifeCycleType lifecycle,
                     const VerboseObject *verbobj)
  : Region(name,id,p0.dim(),lifecycle,verbobj), p0_(p0), p1_(p1)
{

  if (p0_.dim() != p1_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in dimensions of corner points of BoxRegion \"" << Region::name() << "\"\n Perhaps the region is improperly defined?\n";
    if (verbobj && verbobj->os_OK(Teuchos::VERB_MEDIUM)) {
      Teuchos::OSTab tab = verbobj->getOSTab();
      *(verbobj->os()) << tempstr.str();
    }
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }
  
  // Check if this is a reduced dimensionality box (e.g. even though
  // it is in 3D space it is a 2D box)

  int dim = p0.dim();
  for (int i = 0; i < p0.dim(); i++)
    if (p0[i] == p1[i]) dim--;
  
  if (dim < p0.dim()) Region::set_dimension(dim);
}


BoxRegion::BoxRegion(const BoxRegion& old)
  : Region(old), p0_(old.p0_), p1_(old.p1_)
{
  // empty
}

BoxRegion::~BoxRegion(void)
{
  
}

// -------------------------------------------------------------
// BoxRegion::between_
// -------------------------------------------------------------
bool
BoxRegion::between_(const double& x, const double& x0, const double& x1)
{
  double tol = 1.0e-08;
  return  (x0+tol >= x && x >= x1-tol) || (x1+tol >= x && x >= x0-tol);
}

// -------------------------------------------------------------
// BoxRegion::inside
// -------------------------------------------------------------
bool
BoxRegion::inside(const Point& p) const
{

  if (p.dim() != p0_.dim()) {
    std::stringstream tempstr;
    tempstr << "\nMismatch in corner dimension of BoxRegion \"" << Region::name() << "\" and query point.\n Perhaps the region is improperly defined?\n";
    const VerboseObject *verbobj = Region::verbosity_obj();
    if (verbobj && verbobj->os_OK(Teuchos::VERB_MEDIUM)) {
      Teuchos::OSTab tab = verbobj->getOSTab();
      *(verbobj->os()) << tempstr.str();
    }
    Errors::Message mesg(tempstr.str());
    Exceptions::amanzi_throw(mesg);
  }

  bool result(true);
  for (int i = 0; i < p.dim(); ++i) {
    result = result && between_(p[i], p0_[i], p1_[i]);
  }
  return result;
}

// -------------------------------------------------------------
// BoxRegion::is_degenerate (also indicate in how many dimensions)
// -------------------------------------------------------------
bool
BoxRegion::is_degenerate(int *ndeg) const
{
  *ndeg = 0;
  for (int i = 0; i < p0_.dim(); ++i) {
    if (p0_[i] == p1_[i]) (*ndeg)++;    
  }
  if (*ndeg) return true;

  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
