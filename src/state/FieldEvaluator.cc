/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface class for a FieldEvaluator.  A FieldEvaluator is a node in the Phalanx-like
dependency tree.

------------------------------------------------------------------------- */
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "FieldEvaluator.hh"


namespace Amanzi {

FieldEvaluator::FieldEvaluator(Teuchos::ParameterList& plist) : plist_(plist) {
  Teuchos::readVerboseObjectSublist(&plist_,this);
  verbosity_ = getVerbLevel();
  out_ = getOStream();
};



FieldEvaluator::FieldEvaluator(const FieldEvaluator& other) :
    plist_(other.plist_),
    verbosity_(other.verbosity_),
    out_(other.out_) {}

} // namespace
