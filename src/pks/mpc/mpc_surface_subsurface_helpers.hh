#ifndef PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_

#include "CompositeVector.hh"

namespace Amanzi {

void
CopySurfaceToSubsurface(const CompositeVector& surf,
                        const Teuchos::Ptr<CompositeVector>& sub);

void
CopySubsurfaceToSurface(const CompositeVector& sub,
                        const Teuchos::Ptr<CompositeVector>& surf);

void
MergeSubsurfaceAndSurfacePressure(const CompositeVector& kr_surf,
				  const Teuchos::Ptr<CompositeVector>& sub_p,
				  const Teuchos::Ptr<CompositeVector>& surf_p);
double
GetDomainFaceValue(const CompositeVector& sub_p, int f);

void
SetDomainFaceValue(CompositeVector& sub_p, int f, double value);  


} // namespace


#endif

