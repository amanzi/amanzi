/*
  The {evalNameString} evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
{docDict}  
  
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "{evalName}_evaluator.hh"
#include "{evalName}_model.hh"

namespace Amanzi {{
namespace {namespace} {{
namespace Relations {{

// Constructor from ParameterList
{evalClassName}Evaluator::{evalClassName}Evaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{{
  Teuchos::ParameterList& sublist = plist_.sublist("{evalName} parameters");
  model_ = Teuchos::rcp(new {evalClassName}Model(sublist));
  InitializeFromPlist_();
}}


// Copy constructor
{evalClassName}Evaluator::{evalClassName}Evaluator(const {evalClassName}Evaluator& other) :
    SecondaryVariableFieldEvaluator(other),
{keyCopyConstructorList}    
    model_(other.model_) {{}}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
{evalClassName}Evaluator::Clone() const
{{
  return Teuchos::rcp(new {evalClassName}Evaluator(*this));
}}


// Initialize by setting up dependencies
void
{evalClassName}Evaluator::InitializeFromPlist_()
{{
  // Set up my dependencies
  // - defaults to prefixed via domain
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);

  std::string my_key_first("{myKeyFirst}");
  if (domain_name == my_key_first) {{
    domain_name = std::string("");
  }} else {{
    domain_name = domain_name+std::string("_");
  }}

  // - pull Keys from plist
{keyInitializeList}
}}


void
{evalClassName}Evaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{{
{keyCompositeVectorList}

{evaluateModel}
}}


void
{evalClassName}Evaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{{
{keyCompositeVectorList}

{evaluateDerivs}
}}


}} //namespace
}} //namespace
}} //namespace
