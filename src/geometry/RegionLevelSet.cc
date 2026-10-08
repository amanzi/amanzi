/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  A region defined by a level set function.

*/

#include <cmath>

#include "dbc.hh"
#include "errors.hh"
#include "Point.hh"

#include "RegionLevelSet.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// RegionLevelSet:: constructors / destructor
// -------------------------------------------------------------
RegionLevelSet::RegionLevelSet(const std::string& name,
                               const int id,
                               const int dim,
                               const std::string& formula,
                               const LifeCycleType lifecycle)
  : Region(name, id, true, RegionType::LEVELSET, dim, dim, lifecycle)
{
  exprtk_ = std::make_shared<Utils::ExprTK>();
  if (!exprtk_->Initialize(dim + 1, formula)) {
    Errors::Message m;
    m << "Initialization of expression has failed.";
    Exceptions::amanzi_throw(m);
  }
}


// -------------------------------------------------------------
// RegionLevelSet::inside -- check if point is inside region
// -------------------------------------------------------------
bool
RegionLevelSet::inside(const Point& p) const
{
  int d = p.dim();
  std::vector<double> args(4, 0.0);

  args[0] = 0.0;
  for (int i = 0; i < d; ++i) args[i + 1] = p[i];
  double dist = (*exprtk_)(args);

  return dist >= -tol_ ? true : false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
