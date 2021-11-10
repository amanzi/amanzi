/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#include "evaluator/EvaluatorPrimary.hh"

namespace Amanzi {

// registry of method
template <>
Utils::RegisteredFactory<Evaluator, EvaluatorPrimary<double>>
  EvaluatorPrimary<double>::fac_("primary variable double");

template <>
Utils::RegisteredFactory<
  Evaluator, EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>
  EvaluatorPrimary<CompositeVector, CompositeVectorSpace>::fac_(
    "primary variable");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/EvaluatorIndependentFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFunction>
  EvaluatorIndependentFunction::fac_("independent variable");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/EvaluatorIndependentTensorFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentTensorFunction>
  EvaluatorIndependentTensorFunction::fac_("tensor independent variable");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/EvaluatorIndependentFromFile.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentFromFile>
  EvaluatorIndependentFromFile::fac_("independent variable from file");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/EvaluatorIndependentConstant.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorIndependentConstant>
  EvaluatorIndependentConstant::fac_("independent variable constant value");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/Evaluator_OperatorApply.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_OperatorApply>
  Evaluator_OperatorApply::fac_("operator application");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/Evaluator_PDE_Diffusion.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_PDE_Diffusion>
  Evaluator_PDE_Diffusion::fac_("diffusion operator");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/Evaluator_PDE_Accumulation.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_PDE_Accumulation>
  Evaluator_PDE_Accumulation::fac_("accumulation operator");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/Evaluator_BCs.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_BCs>
  Evaluator_BCs::fac_("boundary condition aggregator");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/Evaluator_BCsFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, Evaluator_BCsFunction>
  Evaluator_BCsFunction::fac_("boundary condition function");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/EvaluatorSecondaryMonotypeAdditive.hh"
namespace Amanzi {

template <>
Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeAdditive<
                                      CompositeVector, CompositeVectorSpace>>
  EvaluatorSecondaryMonotypeAdditive<CompositeVector,
                                     CompositeVectorSpace>::fac_("additive");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/EvaluatorSecondaryMonotypeMultiplicative.hh"
namespace Amanzi {

template <>
Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeMultiplicative<
                                      CompositeVector, CompositeVectorSpace>>
  EvaluatorSecondaryMonotypeMultiplicative<
    CompositeVector, CompositeVectorSpace>::fac_("multiplicative");

} // namespace Amanzi
/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "evaluator/EvaluatorSecondaryMonotypeFromFunction.hh"
namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMonotypeFromFunction>
  EvaluatorSecondaryMonotypeFromFunction::fac_(
    "secondary variable from function");

} // namespace Amanzi
#include "evaluator/EvaluatorCellToFace.hh"

Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorCellToFace>
Amanzi::EvaluatorCellToFace::reg_("cell-to-face");
#include "evaluator/EvaluatorCellToFaceUpwind.hh"

Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorCellToFaceUpwind>
Amanzi::EvaluatorCellToFaceUpwind::reg_("cell-to-face-upwind");
#include "evaluator/EvaluatorSecondaryMeshedQuantity.hh"

namespace Amanzi {

template<>
Utils::RegisteredFactory<Evaluator, EvaluatorCellVolume> EvaluatorCellVolume::fac_("cell volume");
template<>
Utils::RegisteredFactory<Evaluator, EvaluatorFaceArea> EvaluatorFaceArea::fac_("face area");
template<>
Utils::RegisteredFactory<Evaluator, EvaluatorMeshElevation> EvaluatorMeshElevation::fac_("meshed elevation");
template<>
Utils::RegisteredFactory<Evaluator, EvaluatorMeshSlopeMagnitude> EvaluatorMeshSlopeMagnitude::fac_("meshed slope magnitude");
}
#include "evaluator/EvaluatorPrimaryStaticMesh.hh"

namespace Amanzi {
Utils::RegisteredFactory<Evaluator, EvaluatorPrimaryStaticMesh> EvaluatorPrimaryStaticMesh::fac_("static mesh");
} 
