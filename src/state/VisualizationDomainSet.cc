/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

//! Visualizes a lifted domain set on the parent mesh.
#include "VisualizationDomainSet.hh"
#include "StringReducer.hh"
#include "DomainSet.hh"

namespace Amanzi {


void
VisualizationDomainSet::writeVector_(const Teuchos::ParameterList& attrs,
                                     const MultiVector_type& vec) const
{
  // replace names[0] domain index with a *
  KeyTriple dset_triple;
  Key name = Keys::cleanPListName(attrs);
  Keys::splitDomainSet(name, dset_triple);
  Key vis_name = Keys::getKey(std::get<0>(dset_triple), std::get<2>(dset_triple));

  if (!lifted_vectors_.count(vis_name)) {
    // create a lifted vector if we don't currently have one
    auto lifted_vec = Teuchos::rcp(new MultiVector_type(
      mesh_->getMap(AmanziMesh::Entity_kind::CELL, false), vec.getNumVectors()));

    // also create a lifted set of attrs
    Teuchos::ParameterList lifted_attrs(attrs);
    lifted_attrs.setName(vis_name);
    lifted_vectors_[vis_name] = std::make_pair(lifted_vec, lifted_attrs);
  }

  // copy from the domain-set vector into the lifted vector
  //
  // Note that to get the domain, we use name_ rather than names[0]'s domain
  // name, as this could be an alias.
  MultiVector_type& lifted_vec = *lifted_vectors_[vis_name].first;
  ds_->doImport(Keys::getDomainInSet(name_, std::get<1>(dset_triple)), vec, lifted_vec);
}


void
VisualizationDomainSet::writeVector_(const Teuchos::ParameterList& attrs,
                                     const Vector_type& vec) const
{
  // replace name domain index with a *
  Key name = Keys::cleanPListName(attrs);
  KeyTriple dset_triple;
  Keys::splitDomainSet(name, dset_triple);
  Key vis_name = Keys::getKey(std::get<0>(dset_triple), std::get<2>(dset_triple));

  if (!lifted_vectors_.count(vis_name)) {
    // create a lifted vector if we don't currently have one
    auto lifted_vec =
      Teuchos::rcp(new MultiVector_type(mesh_->getMap(AmanziMesh::Entity_kind::CELL, false), 1));

    // also create a lifted set of attrs
    Teuchos::ParameterList lifted_attrs(attrs);
    lifted_attrs.setName(vis_name);
    lifted_vectors_[vis_name] = std::make_pair(lifted_vec, lifted_attrs);
  }

  // copy from the domain-set vector into the lifted vector
  MultiVector_type& lifted_vec = *lifted_vectors_[vis_name].first;
  ds_->doImport(Keys::getDomainInSet(name_, std::get<1>(dset_triple)), vec, lifted_vec);
}

void
VisualizationDomainSet::finalizeTimestep()
{
  // construct this the first time, then it is fixed
  if (lifted_vector_names_.size() == 0) {
    // have to get a common set of names across all ranks
    std::vector<std::string> my_names;
    for (auto& lv : lifted_vectors_) { my_names.push_back(lv.first); }

    Utils::StringReducer<100> reducer(mesh_->getComm());
    reducer.checkValidInput(my_names);
    lifted_vector_names_ = reducer.intersectAll(my_names);
  }

  // write the lifted vectors
  for (const auto& vecname : lifted_vector_names_) {
    const auto& vecs = lifted_vectors_.at(vecname);
    if (vecs.first->getNumVectors() == 1) {
      Visualization::write(vecs.second, *vecs.first->getVector(0));
    } else {
      Visualization::write(vecs.second, *vecs.first);
    }
  }

  // clear the lifted vector cache, no need to keep this around as it may be big.
  lifted_vectors_.clear();

  // finalize the files
  Visualization::finalizeTimestep();
}

} // namespace Amanzi
