/*
This is the flow component of the Amanzi code. 
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <set>

#include "errors.hh"
#include "Flow_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* TODO: Verify that a BC has been applied to every boundary face.
* Right now faces without BC are considered no-mass-flux.                                         
****************************************************************** */
void Flow_PK::validate_boundary_conditions(
    BoundaryFunction* bc_pressure, BoundaryFunction* bc_head, BoundaryFunction* bc_flux) const
{
  // Create sets of the face indices belonging to each BC type.
  std::set<int> pressure_faces, head_faces, flux_faces;
  Amanzi::Iterator bc;
  for (bc = bc_pressure->begin(); bc != bc_pressure->end(); ++bc) pressure_faces.insert(bc->first);
  for (bc = bc_head->begin(); bc != bc_head->end(); ++bc) head_faces.insert(bc->first);
  for (bc = bc_flux->begin(); bc != bc_flux->end(); ++bc) flux_faces.insert(bc->first);

  std::set<int> overlap;
  std::set<int>::iterator overlap_end;
  int local_overlap, global_overlap;

  // Check for overlap between pressure and static head BC.
  std::set_intersection(pressure_faces.begin(), pressure_faces.end(),
                        head_faces.begin(), head_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->get_comm()->SumAll(&local_overlap, &global_overlap, 1);  // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "DarcyProblem: static head BC overlap Dirichlet BC on "
        << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }

  // Check for overlap between pressure and flux BC.
  overlap.clear();
  std::set_intersection(pressure_faces.begin(), pressure_faces.end(),
                        flux_faces.begin(), flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->get_comm()->SumAll(&local_overlap, &global_overlap, 1);  // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "DarcyProblem: flux BC overlap Dirichlet BC on "
        << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }

  // Check for overlap between static head and flux BC.
  overlap.clear();
  std::set_intersection(head_faces.begin(), head_faces.end(),
                        flux_faces.begin(), flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->get_comm()->SumAll(&local_overlap, &global_overlap, 1);  // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "DarcyProblem: flux BC overlap static head BC on "
        << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

