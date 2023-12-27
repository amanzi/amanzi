/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

//! Visualizes a lifted domain set on the parent mesh.
#include "Key.hh"
#include "StringReducer.hh"
#include "DomainSet.hh"
#include "OutputDomainSet.hh"

namespace Amanzi {


void
OutputDomainSet::write(const Teuchos::ParameterList& attrs, const MultiVector_type& vec) const
{
  std::vector<std::string> names = OutputUtils::names(attrs, vec.getNumVectors());
  AmanziMesh::Entity_kind entity_kind = attrs.get<AmanziMesh::Entity_kind>("location");
  if (entity_kind != AmanziMesh::CELL) return;

  // replace names[0] domain index with a *
  KeyTriple dset_triple;
  Keys::splitDomainSet(names[0], dset_triple);
  Key vis_name = Keys::getKey(std::get<0>(dset_triple), std::get<2>(dset_triple));

  if (!lifted_vectors_.count(vis_name)) {
    // create a lifted vector if we don't currently have one
    auto lifted_vec = Teuchos::rcp(new MultiVector_type(
      ds_->getReferencingParent()->getMap(AmanziMesh::Entity_kind::CELL, false),
      vec.getNumVectors()));

    // also create a lifted set of names
    std::vector<std::string> lifted_names;
    for (const auto& name : names) {
      KeyTriple split_name;
      Keys::splitDomainSet(name, split_name);
      lifted_names.emplace_back(Keys::getKey(std::get<0>(split_name), std::get<2>(split_name)));
    }
    lifted_vectors_[vis_name] = std::make_pair(lifted_vec, lifted_names);
  }

  // copy from the domain-set vector into the lifted vector
  MultiVector_type& lifted_vec = *lifted_vectors_[vis_name].first;
  ds_->doImport(Keys::getDomainInSet(ds_name_, std::get<1>(dset_triple)), vec, lifted_vec);
}


void
OutputDomainSet::write(const Teuchos::ParameterList& attrs, const IntMultiVector_type& vec) const
{
  std::vector<std::string> names = OutputUtils::names(attrs, vec.getNumVectors());
  AmanziMesh::Entity_kind entity_kind = attrs.get<AmanziMesh::Entity_kind>("location");
  if (entity_kind != AmanziMesh::CELL) return;

  // replace names[0] domain index with a *
  KeyTriple dset_triple;
  Keys::splitDomainSet(names[0], dset_triple);
  Key vis_name = Keys::getKey(std::get<0>(dset_triple), std::get<2>(dset_triple));

  if (!lifted_vectors_.count(vis_name)) {
    // create a lifted vector if we don't currently have one
    auto lifted_vec = Teuchos::rcp(new IntMultiVector_type(
      ds_->getReferencingParent()->getMap(AmanziMesh::Entity_kind::CELL, false),
      vec.getNumVectors()));

    // also create a lifted set of names
    std::vector<std::string> lifted_names;
    for (const auto& name : names) {
      KeyTriple split_name;
      Keys::splitDomainSet(name, split_name);
      lifted_names.emplace_back(Keys::getKey(std::get<0>(split_name), std::get<2>(split_name)));
    }
    lifted_int_vectors_[vis_name] = std::make_pair(lifted_vec, lifted_names);
  }

  // copy from the domain-set vector into the lifted vector
  IntMultiVector_type& lifted_vec = *lifted_int_vectors_[vis_name].first;
  ds_->doImport(Keys::getDomainInSet(ds_name_, std::get<1>(dset_triple)), vec, lifted_vec);
}


void
OutputDomainSet::finalizeTimestep()
{
  if (lifted_vector_names_.size() == 0) {
    // have to get a common set of names across all ranks
    std::vector<std::string> my_names;
    for (auto& lv : lifted_vector_names_) { my_names.push_back(lv); }

    Utils::StringReducer<100> reducer(ds_->getComm());
    reducer.checkValidInput(my_names);
    lifted_vector_names_ = reducer.intersectAll(my_names);
  }

  if (lifted_int_vector_names_.size() == 0) {
    // have to get a common set of names across all ranks
    std::vector<std::string> my_names;
    for (auto& lv : lifted_int_vector_names_) { my_names.push_back(lv); }

    Utils::StringReducer<100> reducer(ds_->getComm());
    reducer.checkValidInput(my_names);
    lifted_int_vector_names_ = reducer.intersectAll(my_names);
  }

  // write the lifted vectors
  for (const auto& vecname : lifted_vector_names_) {
    const auto& vecs = lifted_vectors_.at(vecname);
    Teuchos::ParameterList plist(vecname);
    plist.set<Teuchos::Array<std::string>>("subfield names", vecs.second);
    output_->write(plist, *vecs.first);
  }

  // and the int vectors
  for (const auto& vecname : lifted_int_vector_names_) {
    const auto& vecs = lifted_int_vectors_.at(vecname);
    Teuchos::ParameterList plist(vecname);
    plist.set<Teuchos::Array<std::string>>("subfield names", vecs.second);
    output_->write(plist, *vecs.first);
  }

  // clear the lifted vector cache, no need to keep this around as it may be big.
  lifted_vectors_.clear();
  lifted_int_vectors_.clear();

  // keeping these implies that state is the same across usage of this object.
  // If it changes, clear these.
  // lifted_vector_names_.clear();
  // lifted_int_vector_names_.clear();

  // finalize the files
  output_->finalizeTimestep();
}

} // namespace Amanzi
