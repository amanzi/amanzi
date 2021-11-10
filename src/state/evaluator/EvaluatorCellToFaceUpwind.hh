/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/
//! Calculates a face value from cell values using upwind algorithm.

#pragma once

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "EvaluatorSecondary.hh"
#include "Factory.hh"
#include "Key.hh"

namespace Amanzi {

class EvaluatorCellToFaceUpwind : public EvaluatorSecondary {
 public:
  explicit EvaluatorCellToFaceUpwind(Teuchos::ParameterList& plist);

  EvaluatorCellToFaceUpwind(const EvaluatorCellToFaceUpwind& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EvaluatorCellToFaceUpwind(*this));
  }

  virtual Evaluator& operator=(const Evaluator& other) override;
  EvaluatorCellToFaceUpwind& operator=(const EvaluatorCellToFaceUpwind& other) = default;

  virtual void EnsureCompatibility(State& S) override;
  
  virtual void Update_(State& S) override;

  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) override;

  static Utils::RegisteredFactory<Evaluator, EvaluatorCellToFaceUpwind> reg_;

 private:
  Key flux_key_;
};
  
} // namespace Amanzi
