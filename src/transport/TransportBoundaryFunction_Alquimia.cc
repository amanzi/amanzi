/*
  This is the transport component of the Amanzi code. 

  License: see $AMANZI_DIR/COPYRIGHT
  Author (v1): Neil Carlson
         (v2): Ethan Coon
*/

#include "TransportBoundaryFunction_Alquimia.hh"

#ifdef ALQUIMIA_ENABLED

namespace Amanzi {
namespace AmanziTransport {

/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Define(const std::vector<std::string> &regions)
{
}


/* ******************************************************************
* Internal subroutine that defines a boundary function.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Define(std::string region) 
{
  RegionList regions(1,region);
}


/* ******************************************************************
* Evaluate values at time.
****************************************************************** */
void TransportBoundaryFunction_Alquimia::Compute(double time) 
{
  // Loop over sides and evaluate values.
  for (int n = 0; n < faces_.size(); n++) {
    int f = faces_[n];
    AmanziGeometry::Point xc = mesh_->face_centroid(f);

    for (int i = 0; i < tcc_index_.size(); i++) {
      values_[n][i] = 0.0;
    }
  }
}

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

