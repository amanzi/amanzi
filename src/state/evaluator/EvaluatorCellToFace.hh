/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/
//! Calculates a face value from cell values.

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "EvaluatorSecondary.hh"
#include "Factory.hh"

namespace Amanzi {

class EvaluatorCellToFace : public EvaluatorSecondary {
 public:
  explicit EvaluatorCellToFace(Teuchos::ParameterList& plist);

  EvaluatorCellToFace(const EvaluatorCellToFace& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EvaluatorCellToFace(*this));
  }

  virtual Evaluator& operator=(const Evaluator& other) override;
  EvaluatorCellToFace& operator=(const EvaluatorCellToFace& other) = default;

  virtual void EnsureCompatibility(State& S) override;
  
  virtual void Update_(State& S) override;

  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) override;

protected:


  std::string algorithm_;

 private:

  static Utils::RegisteredFactory<Evaluator, EvaluatorCellToFace> reg_;
  
};
  
} // namespace Amanzi
