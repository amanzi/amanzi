/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
   ATS

   Adding an ATS to Thyra converter.

   Non-member functions for conversion between Tree/CompositeVectors and Thyra
   data structures.

*/

#ifndef AMANZI_THYRA_WRAPPERS_HH_
#define AMANZI_THYRA_WRAPPERS_HH_

#include "Teuchos_RCP.hpp"
#include "Thyra_VectorBase.hpp"

#include "CompositeVector.hh"
#include "TreeVector.hh"

namespace Amanzi {
namespace ThyraWrappers {

// Converts a TreeVector into a ThyraVector.
Teuchos::RCP<const Thyra::VectorBase<double>>
CreateThyraVector(const Teuchos::RCP<const TreeVector>& tv);

// Converts a TreeVector into a ThyraVector, nonconst version.
Teuchos::RCP<Thyra::VectorBase<double>>
CreateThyraVector(const Teuchos::RCP<TreeVector>& tv);

// Converts a CompositeVector into a ThyraVector.
Teuchos::RCP<const Thyra::VectorBase<double>>
CreateThyraVector(const Teuchos::RCP<const CompositeVector>& cv);

// Converts a CompositeVector into a ThyraVector, nonconst version.
Teuchos::RCP<Thyra::VectorBase<double>>
CreateThyraVector(const Teuchos::RCP<CompositeVector>& cv);

// Creates a Thyra vector space associated with a TreeVector's structure.
Teuchos::RCP<Thyra::VectorSpaceBase<double>>
CreateThyraVectorSpace(const Teuchos::RCP<const TreeVector>& tv);

// Creates a Thyra vector space associated with a TreeVector's structure, nonconst version.
Teuchos::RCP<Thyra::VectorSpaceBase<double>>
CreateThyraVectorSpace(const Teuchos::RCP<TreeVector>& tv);

// Creates a Thyra vector space associated with a CompositeVector's structure.
Teuchos::RCP<Thyra::VectorSpaceBase<double>>
CreateThyraVectorSpace(const Teuchos::RCP<const CompositeVector>& cv);

// Creates a Thyra vector space associated with a CompositeVector's structure, nonconst version.
Teuchos::RCP<Thyra::VectorSpaceBase<double>>
CreateThyraVectorSpace(const Teuchos::RCP<CompositeVector>& cv);

// Create a TreeVector with structure from the structurevector and data
// from the Thyra vector.
Teuchos::RCP<TreeVector>
CreateTreeVector(std::string name,
                 const Teuchos::RCP<Thyra::VectorBase<double>>& vec,
                 const Teuchos::RCP<const TreeVector>& structure);

void
ViewThyraVectorAsTreeVector(const Teuchos::RCP<Thyra::VectorBase<double>>& vec,
                            const Teuchos::RCP<TreeVector>& tv);

// Create a CompositeVector with structure from the structurevector and data
// from the Thyra vector.
Teuchos::RCP<CompositeVector>
CreateCompositeVector(const Teuchos::RCP<Thyra::VectorBase<double>>& vec,
                      const Teuchos::RCP<const CompositeVector>& structure);

void
ViewThyraVectorAsCompositeVector(const Teuchos::RCP<Thyra::VectorBase<double>>& vec,
                                 const Teuchos::RCP<CompositeVector>& cv);

} // namespace ThyraWrappers
} // namespace Amanzi

#endif
