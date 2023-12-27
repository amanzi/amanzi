/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

//! Visualizes a lifted domain set on the parent mesh.
/*!

Domain sets are collections of subdomains that, most commonly, partition a
paraent domain.  Visualizing this collection of objects can be tricky, but is
made much easier by first migrating all the domain sets on to the parent mesh
and visualizing on that mesh.

*/

#pragma once

#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"

#include "Output.hh"

namespace Amanzi {

namespace AmanziMesh {
class DomainSet;
}

class OutputDomainSet : public Output {
 public:
  OutputDomainSet(const std::string& ds_name,
                  const Teuchos::RCP<const AmanziMesh::DomainSet>& ds,
                  const Teuchos::RCP<Output>& output)
    : ds_name_(ds_name), ds_(ds), output_(output)
  {}

  // public interface for data clients
  virtual void
  write(const Teuchos::ParameterList& attrs, const MultiVector_type& vec) const override;
  virtual void
  write(const Teuchos::ParameterList& attrs, const IntMultiVector_type& vec) const override;
  virtual void finalizeTimestep() override;

 protected:
  // note this is lazily constructed, so must be mutable
  // Note that this relies on map being ORDERED!
  mutable std::map<std::string, std::pair<Teuchos::RCP<MultiVector_type>, std::vector<std::string>>>
    lifted_vectors_;
  mutable std::vector<std::string> lifted_vector_names_;
  mutable std::map<std::string,
                   std::pair<Teuchos::RCP<IntMultiVector_type>, std::vector<std::string>>>
    lifted_int_vectors_;
  mutable std::vector<std::string> lifted_int_vector_names_;

  Teuchos::RCP<const AmanziMesh::DomainSet> ds_;
  std::string ds_name_;

  Teuchos::RCP<Output> output_;
};

} // namespace Amanzi
