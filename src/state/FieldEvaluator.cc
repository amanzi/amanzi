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

FieldEvaluator::FieldEvaluator(Teuchos::ParameterList& plist) :
    plist_(plist),
    type_(EvaluatorType::UNKNOWN) {
  vo_ = Teuchos::rcp(new VerboseObject(Keys::cleanPListName(plist.name()), plist));
};



FieldEvaluator::FieldEvaluator(const FieldEvaluator& other) :
    plist_(other.plist_),
    type_(other.type_),
    vo_(other.vo_) {}

std::ostream&
operator<<(std::ostream& os, const FieldEvaluator& self) {
  return os << self.WriteToString();
}

} // namespace
