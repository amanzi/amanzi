#include "mesh_function.hh"

#include <algorithm>
#include "errors.hh"

namespace Amanzi {

/* ****************************************************************
* Populates internal array with function values.
**************************************************************** */
void MeshFunction::Compute(double t)
{
  int dim = (*mesh_).space_dimension();
  double *args = new double[1+dim];
  double *xargs = args+1;
  args[0] = t;

  for (SpecList::const_iterator s = spec_list_.begin(); s != spec_list_.end(); ++s) {
    // Here we could specialize on the argument signature of the function: 
    // time-independent functions need only be evaluated at each face on the
    // first call; space-independent functions need only be evaluated once per
    // call and the value used for all faces; etc. Right now we just assume
    // the most general case.
    const Domain &domain = s->first;
    for (Domain::const_iterator d=domain.begin(); d != domain.end(); ++d) {
      AmanziGeometry::Point xc = (*mesh_).face_centroid(*d);
      for (int i=0; i<dim; ++i) xargs[i] = xc[i];
      value_[*d] = (*(s->second))(args);
    }
  }
  delete [] args;
}


/* ****************************************************************
* Similar to Define() but uses a single string.
**************************************************************** */
void MeshFunction::DefineFromString(const std::string region, const Teuchos::RCP<const Function>& f)
{
  std::vector<std::string> regions(1, region);
  Define(regions, f);
}

} // namespace Amanzi
