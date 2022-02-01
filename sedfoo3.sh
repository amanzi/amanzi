function go3 () {
    SED=gsed
    SED_ARGS=-i''

    echo "Running go3 on $1"

    # # RequireField("myfield") --> Require<CV,CVS>("myfield",Tags::NEXT)
    # ${SED} ${SED_ARGS} -E 's|RequireField\(([^,\)]*)\)->|Require<CompositeVector,CompositeVectorSpace>(\1, Tags::NEXT).|g' "$1"
    # ${SED} ${SED_ARGS} -E 's|RequireField\(([^,\)]*)\)|Require<CompositeVector,CompositeVectorSpace>(\1, Tags::NEXT)|g' "$1"

    # # RequireField("myfield", "myowner") --> Require<CV,CVS>("myfield",Tags::NEXT,"myowner")
    # ${SED} ${SED_ARGS} -E 's|RequireField\(([^,\)]*),([^,\)]*)\)->|Require<CompositeVector,CompositeVectorSpace>(\1, Tags::NEXT, \2).|g' "$1"
    # ${SED} ${SED_ARGS} -E 's|RequireField\(([^,\)]*),([^,\)]*)\)|Require<CompositeVector,CompositeVectorSpace>(\1, Tags::NEXT, \2)|g' "$1"

    # # GetFieldData() --> GetPtr()
    # ${SED} ${SED_ARGS} -E 's|GetFieldData\(([^,\)]*)\)->|Get<CompositeVector>(\1).|g' "$1"
    # ${SED} ${SED_ARGS} -E 's|GetFieldData\(([^,\)]*)\)|GetPtr<CompositeVector>(\1)|g' "$1"
    # ${SED} ${SED_ARGS} -E 's|GetFieldData\(([^,\)]*),([^,\)]*)\)->|GetW<CompositeVector>(\1,\2).|g' "$1"
    # ${SED} ${SED_ARGS} -E 's|GetFieldData\(([^,\)]*),([^,\)]*)\)|GetPtrW<CompositeVector>(\1,\2)|g' "$1"    

    # secondary_variable_field_eval
    ${SED} ${SED_ARGS} 's|FieldEvaluator_Factory\.hh|Evaluator_Factory.hh|' "$1"
    ${SED} ${SED_ARGS} 's|secondary_variable_field_evaluator\.hh|EvaluatorSecondaryMonotype.hh|' "$1"
    ${SED} ${SED_ARGS} 's|secondary_variables_field_evaluator\.hh|EvaluatorSecondaryMonotype.hh|' "$1"
    ${SED} ${SED_ARGS} 's|SecondaryVariablesFieldEvaluator|EvaluatorSecondaryMonotypeCV|' "$1"
    ${SED} ${SED_ARGS} 's|SecondaryVariableFieldEvaluator|EvaluatorSecondaryMonotypeCV|' "$1"
    ${SED} ${SED_ARGS} 's|<FieldEvaluator>|<Evaluator>|' "$1"
    ${SED} ${SED_ARGS} 's|my_keys_\.push_back|my_keys_\.emplace_back(|' "$1"

    
}
