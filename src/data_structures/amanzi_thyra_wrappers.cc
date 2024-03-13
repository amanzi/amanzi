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

*/

#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "amanzi_thyra_wrappers.hh"

namespace Amanzi {
namespace ThyraWrappers {

// Converts a TreeVector into a ThyraVector.
Teuchos::RCP<Thyra::VectorBase<double>>
CreateThyraVector(const Teuchos::RCP<TreeVector>& tv)
{
  // A tree vector may have EITHER subvecs or CompositeVector data.
  if (tv->data() != Teuchos::null) {
    // has data, simply return the CV's vector.
    return CreateThyraVector(tv->data());
  }

  // If no data, grab subvecs.
  std::vector<Teuchos::RCP<TreeVector>> tv_subvecs = tv->SubVectors();

  // If no data and no subvecs, return NULL
  if (tv_subvecs.size() == 0) { return Teuchos::null; }

  // Loop over subvecs, collecting the Thyra Vectors and VectorSpaces
  // associated with each subvec.
  std::vector<Teuchos::RCP<Thyra::VectorBase<double>>> subvecs;
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces;
  for (std::vector<Teuchos::RCP<TreeVector>>::iterator tv_subvec = tv_subvecs.begin();
       tv_subvec != tv_subvecs.end();
       ++tv_subvec) {
    Teuchos::RCP<Thyra::VectorBase<double>> subvec = CreateThyraVector(*tv_subvec);
    if (subvec != Teuchos::null) {
      subvecs.push_back(subvec);
      subspaces.push_back(subvec->space());
    }
  }

  // Convert the list of subspaces into a ProductVectorSpace
  Teuchos::ArrayView<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces_arrayview(
    subspaces);
  Teuchos::RCP<Thyra::DefaultProductVectorSpace<double>> space =
    Thyra::productVectorSpace<double>(subspaces_arrayview);

  // Convert the list of vectors into a ProductVector
  Teuchos::ArrayView<Teuchos::RCP<Thyra::VectorBase<double>>> subvecs_arrayview(subvecs);
  Teuchos::RCP<Thyra::DefaultProductVector<double>> vec =
    Thyra::defaultProductVector<double>(space, subvecs_arrayview);
  return vec;
};


// Converts a TreeVector into a ThyraVector, const version.
Teuchos::RCP<const Thyra::VectorBase<double>>
CreateThyraVector(const Teuchos::RCP<const TreeVector>& tv)
{
  // A tree vector may have EITHER subvecs or CompositeVector data.
  if (tv->data() != Teuchos::null) {
    // has data, simply return the CV's vector.
    return CreateThyraVector(tv->data());
  }

  // If no data, grab subvecs.
  const std::vector<Teuchos::RCP<TreeVector>> tv_subvecs = tv->SubVectors();

  // If no data and no subvecs, return NULL
  if (tv_subvecs.size() == 0) { return Teuchos::null; }

  // Loop over subvecs, collecting the Thyra Vectors and VectorSpaces
  // associated with each subvec.
  std::vector<Teuchos::RCP<const Thyra::VectorBase<double>>> subvecs;
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces;
  for (std::vector<Teuchos::RCP<TreeVector>>::const_iterator tv_subvec = tv_subvecs.begin();
       tv_subvec != tv_subvecs.end();
       ++tv_subvec) {
    Teuchos::RCP<const Thyra::VectorBase<double>> subvec = CreateThyraVector(*tv_subvec);
    if (subvec != Teuchos::null) {
      subvecs.push_back(subvec);
      subspaces.push_back(subvec->space());
    }
  }

  // Convert the list of subspaces into a ProductVectorSpace
  Teuchos::ArrayView<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces_arrayview(
    subspaces);
  Teuchos::RCP<Thyra::DefaultProductVectorSpace<double>> space =
    Thyra::productVectorSpace<double>(subspaces_arrayview);

  // Convert the list of vectors into a ProductVector
  Teuchos::ArrayView<Teuchos::RCP<const Thyra::VectorBase<double>>> subvecs_arrayview(subvecs);
  Teuchos::RCP<const Thyra::DefaultProductVector<double>> vec =
    Thyra::defaultProductVector<double>(space, subvecs_arrayview);
  return vec;
};


// Converts a CompositeVector into a ThyraVector.
Teuchos::RCP<Thyra::VectorBase<double>>
CreateThyraVector(const Teuchos::RCP<CompositeVector>& cv)
{
  // If the CompositeVector is empty, return null.
  if (cv->num_components() == 0) { return Teuchos::null; }

  // Loop over components, collecting the Thyra Vectors and VectorSpaces
  // associated with each component.
  std::vector<Teuchos::RCP<Thyra::VectorBase<double>>> subvecs;
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces;
  for (CompositeVector::name_iterator name = cv->begin(); name != cv->end(); ++name) {
    // The strategy here is this -- first we wrap the Epetra MV that is the
    // component as a Thyra MultiVector.  This MultiVector can then be
    // converted into a single Thyra Vector on a ProductVectorSpace by taking
    // the outer product of the map's space num_dof times.

    // Create the VectorSpace associated with the component's (non-ghosted) map.
    Teuchos::RCP<const Epetra_Map> map = cv->map(*name, false);
    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> comp_space = Thyra::create_VectorSpace(map);

    // Create a MultiVector using the Epetra_MV data.
    Teuchos::RCP<Thyra::MultiVectorBase<double>> comp_mv =
      Thyra::create_MultiVector(cv->ViewComponent(*name, false), comp_space);

    // Create the component's ProductVectorSpace
    Teuchos::RCP<Thyra::DefaultMultiVectorProductVectorSpace<double>> comp_pv_space =
      Thyra::multiVectorProductVectorSpace<double>(comp_space, cv->num_dofs(*name));

    // Create the ProductVector for the component.
    Teuchos::RCP<Thyra::DefaultMultiVectorProductVector<double>> comp_pv =
      Thyra::multiVectorProductVector<double>(comp_pv_space, comp_mv);

    subvecs.push_back(comp_pv);
    subspaces.push_back(comp_pv_space);
  }

  // Convert the list of subspaces into a ProductVectorSpace
  Teuchos::ArrayView<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces_arrayview(
    subspaces);
  Teuchos::RCP<Thyra::DefaultProductVectorSpace<double>> space =
    Thyra::productVectorSpace<double>(subspaces_arrayview);

  // Convert the list of vectors into a ProductVector
  Teuchos::ArrayView<Teuchos::RCP<Thyra::VectorBase<double>>> subvecs_arrayview(subvecs);
  Teuchos::RCP<Thyra::DefaultProductVector<double>> vec =
    Thyra::defaultProductVector<double>(space, subvecs_arrayview);
  return vec;
};

// Converts a CompositeVector into a ThyraVector, const version.
Teuchos::RCP<const Thyra::VectorBase<double>>
CreateThyraVector(const Teuchos::RCP<const CompositeVector>& cv)
{
  // If the CompositeVector is empty, return null.
  if (cv->num_components() == 0) { return Teuchos::null; }

  // Loop over components, collecting the Thyra Vectors and VectorSpaces
  // associated with each component.
  std::vector<Teuchos::RCP<const Thyra::VectorBase<double>>> subvecs;
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces;
  for (CompositeVector::name_iterator name = cv->begin(); name != cv->end(); ++name) {
    // The strategy here is this -- first we wrap the Epetra MV that is the
    // component as a Thyra MultiVector.  This MultiVector can then be
    // converted into a single Thyra Vector on a ProductVectorSpace by taking
    // the outer product of the map's space num_dof times.

    // Create the VectorSpace associated with the component's (non-ghosted) map.
    Teuchos::RCP<const Epetra_Map> map = cv->map(*name, false);
    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> comp_space = Thyra::create_VectorSpace(map);

    // Create a MultiVector using the Epetra_MV data.
    Teuchos::RCP<const Thyra::MultiVectorBase<double>> comp_mv =
      Thyra::create_MultiVector(cv->ViewComponent(*name, false), comp_space);

    // Create the component's ProductVectorSpace
    Teuchos::RCP<const Thyra::DefaultMultiVectorProductVectorSpace<double>> comp_pv_space =
      Thyra::multiVectorProductVectorSpace<double>(comp_space, cv->num_dofs(*name));

    // Create the ProductVector for the component.
    Teuchos::RCP<const Thyra::DefaultMultiVectorProductVector<double>> comp_pv =
      Thyra::multiVectorProductVector<double>(comp_pv_space, comp_mv);

    subvecs.push_back(comp_pv);
    subspaces.push_back(comp_pv_space);
  }

  // Convert the list of subspaces into a ProductVectorSpace
  Teuchos::ArrayView<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces_arrayview(
    subspaces);
  Teuchos::RCP<const Thyra::DefaultProductVectorSpace<double>> space =
    Thyra::productVectorSpace<double>(subspaces_arrayview);

  // Convert the list of vectors into a ProductVector
  Teuchos::ArrayView<Teuchos::RCP<const Thyra::VectorBase<double>>> subvecs_arrayview(subvecs);
  Teuchos::RCP<const Thyra::DefaultProductVector<double>> vec =
    Thyra::defaultProductVector<double>(space, subvecs_arrayview);
  return vec;
};


// Creates a Thyra vector space associated with a TreeVector's structure.
Teuchos::RCP<Thyra::VectorSpaceBase<double>>
CreateThyraVectorSpace(const Teuchos::RCP<const TreeVector>& tv)
{
  // A tree vector may have EITHER subvecs or CompositeVector data.
  if (tv->data() != Teuchos::null) {
    // has data, simply return the CV's vector.
    return CreateThyraVectorSpace(tv->data());
  }

  // If no data, grab subvecs.
  const std::vector<Teuchos::RCP<TreeVector>> tv_subvecs = tv->SubVectors();

  // If no data and no subvecs, return NULL
  if (tv_subvecs.size() == 0) { return Teuchos::null; }

  // Loop over subvecs, collecting the Thyra Vectors and VectorSpaces
  // associated with each subvec.
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces;
  for (std::vector<Teuchos::RCP<TreeVector>>::const_iterator tv_subvec = tv_subvecs.begin();
       tv_subvec != tv_subvecs.end();
       ++tv_subvec) {
    Teuchos::RCP<Thyra::VectorSpaceBase<double>> subspace = CreateThyraVectorSpace(*tv_subvec);
    if (subspace != Teuchos::null) { subspaces.push_back(subspace); }
  }

  // Convert the list of subspaces into a ProductVectorSpace
  Teuchos::ArrayView<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces_arrayview(
    subspaces);
  Teuchos::RCP<Thyra::DefaultProductVectorSpace<double>> space =
    Thyra::productVectorSpace<double>(subspaces_arrayview);
  return space;
};


// Creates a Thyra vector space associated with a TreeVector's structure.
Teuchos::RCP<Thyra::VectorSpaceBase<double>>
CreateThyraVectorSpace(const Teuchos::RCP<TreeVector>& tv)
{
  Teuchos::RCP<const TreeVector> tv_tmp(tv);
  return CreateThyraVectorSpace(tv_tmp);
};


// Creates a Thyra vector space associated with a CompositeVector's structure.
Teuchos::RCP<Thyra::VectorSpaceBase<double>>
CreateThyraVectorSpace(const Teuchos::RCP<const CompositeVector>& cv)
{
  // If the CompositeVector is empty, return null.
  if (cv->num_components() == 0) { return Teuchos::null; }

  // Loop over components, collecting the Thyra Vectors and VectorSpaces
  // associated with each component.
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces;
  for (CompositeVector::name_iterator name = cv->begin(); name != cv->end(); ++name) {
    // The strategy here is this -- first we wrap the Epetra MV that is the
    // component as a Thyra MultiVector.  This MultiVector can then be
    // converted into a single Thyra Vector on a ProductVectorSpace by taking
    // the outer product of the map's space num_dof times.

    // Create the VectorSpace associated with the component's (non-ghosted) map.
    Teuchos::RCP<const Epetra_Map> map = cv->map(*name, false);
    Teuchos::RCP<const Thyra::VectorSpaceBase<double>> comp_space = Thyra::create_VectorSpace(map);

    // Create the component's ProductVectorSpace
    Teuchos::RCP<Thyra::DefaultMultiVectorProductVectorSpace<double>> comp_pv_space =
      Thyra::multiVectorProductVectorSpace<double>(comp_space, cv->num_dofs(*name));

    subspaces.push_back(comp_pv_space);
  }

  // Convert the list of subspaces into a ProductVectorSpace
  Teuchos::ArrayView<Teuchos::RCP<const Thyra::VectorSpaceBase<double>>> subspaces_arrayview(
    subspaces);
  Teuchos::RCP<Thyra::DefaultProductVectorSpace<double>> space =
    Thyra::productVectorSpace<double>(subspaces_arrayview);
  return space;
};


// Creates a Thyra vector space associated with a CompositeVector's structure.
Teuchos::RCP<Thyra::VectorSpaceBase<double>>
CreateThyraVectorSpace(const Teuchos::RCP<CompositeVector>& cv)
{
  Teuchos::RCP<const CompositeVector> cv_tmp(cv);
  return CreateThyraVectorSpace(cv_tmp);
};


// Create a TreeVector with structure from the structurevector and data
// from the Thyra vector.
Teuchos::RCP<TreeVector>
CreateTreeVector(std::string name,
                 const Teuchos::RCP<Thyra::VectorBase<double>>& vec,
                 const Teuchos::RCP<const TreeVector>& structure)
{
  Teuchos::RCP<TreeVector> tv =
    Teuchos::rcp(new TreeVector(name, *structure, CONSTRUCT_WITHOUT_DATA));
  ViewThyraVectorAsTreeVector(vec, tv);
  return tv;
};

// Create a CompositeVector with structure from the structurevector and data
// from the Thyra vector.
Teuchos::RCP<CompositeVector>
CreateCompositeVector(const Teuchos::RCP<Thyra::VectorBase<double>>& vec,
                      const Teuchos::RCP<const CompositeVector>& structure)
{
  Teuchos::RCP<CompositeVector> cv =
    Teuchos::rcp(new CompositeVector(*structure, CONSTRUCT_WITHOUT_DATA));
  ViewThyraVectorAsCompositeVector(vec, cv);
  return cv;
};


void
ViewThyraVectorAsTreeVector(const Teuchos::RCP<Thyra::VectorBase<double>>& vec,
                            const Teuchos::RCP<TreeVector>& tv)
{
  // A tree vector may have EITHER subvecs or CompositeVector data.
  if (tv->data() != Teuchos::null) {
    // has data, simply return the CV's vector.
    ViewThyraVectorAsCompositeVector(vec, tv->data());
    return;
  }

  // If no data, grab subvecs.
  const std::vector<Teuchos::RCP<TreeVector>> tv_subvecs = tv->SubVectors();

  // If no data and no subvecs, return.
  if (tv_subvecs.size() == 0) { return; }

  // Otherwise the vector is a ProductVector for SubVectors in the TV.
  // Unpack vec as ProductVector
  Teuchos::RCP<Thyra::DefaultProductVector<double>> pvec =
    Teuchos::rcp_dynamic_cast<Thyra::DefaultProductVector<double>>(vec);
  AMANZI_ASSERT(vec != Teuchos::null);

  // Loop over subvectors.
  int i = 0;
  for (std::vector<Teuchos::RCP<TreeVector>>::const_iterator tv_subvec = tv_subvecs.begin();
       tv_subvec != tv_subvecs.end();
       ++tv_subvec) {
    ViewThyraVectorAsTreeVector(pvec->getNonconstVectorBlock(i), *tv_subvec);
    i++;
  }
  return;
};


void
ViewThyraVectorAsCompositeVector(const Teuchos::RCP<Thyra::VectorBase<double>>& vec,
                                 const Teuchos::RCP<CompositeVector>& cv)
{
  // The thyra vec is NOT ghosted -- if the CV is ghosted, we need to make
  // sure it has been created.
  if (cv->ghosted() && !cv->created()) { cv->CreateData(); }

  // Unpack vec as ProductVector
  Teuchos::RCP<Thyra::DefaultProductVector<double>> pvec =
    Teuchos::rcp_dynamic_cast<Thyra::DefaultProductVector<double>>(vec);
  AMANZI_ASSERT(vec != Teuchos::null);

  // Loop over components, assigning the block's vec to the component Epetra_MultiVector
  int i = 0;
  for (CompositeVector::name_iterator name = cv->begin(); name != cv->end(); ++name) {
    Teuchos::RCP<Thyra::VectorBase<double>> block = pvec->getNonconstVectorBlock(i);

    // each block vector is actually a product
    Teuchos::RCP<Thyra::DefaultMultiVectorProductVector<double>> pv_block =
      Teuchos::rcp_dynamic_cast<Thyra::DefaultMultiVectorProductVector<double>>(block);
    AMANZI_ASSERT(pv_block != Teuchos::null);

    // get the MultiVector data
    Teuchos::RCP<Thyra::MultiVectorBase<double>> mv_block = pv_block->getNonconstMultiVector();

    // get an Epetra_MultiVector view of the data and stick it in the CV
    Teuchos::RCP<Epetra_MultiVector> mv =
      Thyra::get_Epetra_MultiVector(*cv->map(*name, false), mv_block);
    cv->SetComponent(*name, mv);
    i++;
  }
  return;
};


} // namespace ThyraWrappers
} // namespace Amanzi
