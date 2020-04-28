/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//!

/*!

This is a simple evaluator for a non-deforming mesh.  It can provide a target
for mesh dependencies of vectors (e.g. cell_volume, elevation, others) that
might potentially change in a deforming mesh case, but (in this case) don't.

*/

#pragma once

#include "EvaluatorPrimary.hh"

namespace Amanzi {

class EvaluatorPrimaryStaticMesh
    : public EvaluatorPrimary<double> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors, assignement operators, etc
  // ---------------------------------------------------------------------------
  EvaluatorPrimaryStaticMesh(Teuchos::ParameterList& plist)
      : EvaluatorPrimary<double>(plist),
        domain_(plist.get<std::string>("mesh name")) {}

  EvaluatorPrimaryStaticMesh(const EvaluatorPrimaryStaticMesh& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EvaluatorPrimaryStaticMesh(*this));
  }

  virtual std::string name() const override { return "static mesh"; }

 protected:
  std::string domain_;
  
 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorPrimaryStaticMesh> fac_;
};

} // namespace Amanzi

