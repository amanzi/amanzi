function go_tpetra () {
    echo "Running go_tpetra on $1"

    # comm changes
    gsed -i 's|->MyPID()|->getRank()|g' "$1"
    gsed -i 's|\.MyPID()|.getRank()|g' "$1"
    gsed -i 's|->NumProcs()|->getSize()|g' "$1"
    gsed -i 's|\.NumProcs()|.getSize()|g' "$1"
    
    # changes all reductions to Teuchos reductions
    gsed -i 's|\([^ ]*\)->MaxAll(\([^,]*\),\([^,]*\),\(.*\));|Teuchos::reduceAll<int>(*\1,Teuchos::REDUCE_MAX,\4,\2,\3);|g' "$1"
    gsed -i 's|\([^ ]*\)\.MaxAll(\([^,]*\),\([^,]*\),\(.*\));|Teuchos::reduceAll<int>(\1,Teuchos::REDUCE_MAX,\4,\2,\3);|g' "$1"
    gsed -i 's|\([^ ]*\)->MinAll(\([^,]*\),\([^,]*\),\(.*\));|Teuchos::reduceAll<int>(*\1,Teuchos::REDUCE_MIN,\4,\2,\3);|g' "$1"
    gsed -i 's|\([^ ]*\)\.MinAll(\([^,]*\),\([^,]*\),\(.*\));|Teuchos::reduceAll<int>(\1,Teuchos::REDUCE_MIN,\4,\2,\3);|g' "$1"
    gsed -i 's|\([^ ]*\)->SumAll(\([^,]*\),\([^,]*\),\(.*\));|Teuchos::reduceAll<int>(*\1,Teuchos::REDUCE_SUM,\4,\2,\3);|g' "$1"
    gsed -i 's|\([^ ]*\)\.SumAll(\([^,]*\),\([^,]*\),\(.*\));|Teuchos::reduceAll<int>(\1,Teuchos::REDUCE_SUM,\4,\2,\3);|g' "$1"

    # capitalization changes
    gsed -i 's|->HasComponent|->hasComponent|g' "$1"
    gsed -i 's|->NumVectors|->getNumVectors|g' "$1"
    gsed -i 's|->MyLength|->getLocalLength|g' "$1"
    gsed -i 's|\.HasComponent|.hasComponent|g' "$1"
    gsed -i 's|->GetComponent(|->getComponent(|g' "$1"
    gsed -i 's|\.GetComponent(|.getComponent(|g' "$1"
    gsed -i 's|\.Mesh()|.getMesh()|g' "$1"
    gsed -i 's|->Mesh()|->getMesh()|g' "$1"
    gsed -i 's|\.Map(|.getMap(|g' "$1"
    gsed -i 's|->Map(|->getMap(|g' "$1"

    # tpetra vector changes
    gsed -i 's|\.NumVectors|.getNumVectors|g' "$1"
    gsed -i 's|\.MyLength|.getLocalLength|g' "$1"

    # tpetra map changes
    gsed -i 's|->GID|->getGlobalElement|g' "$1"
    gsed -i 's|\.GID|.getGlobalElement|g' "$1"
    gsed -i 's|->LID|->getLocalElement|g' "$1"
    gsed -i 's|\.LID|.getLocalElement|g' "$1"
}    
    
