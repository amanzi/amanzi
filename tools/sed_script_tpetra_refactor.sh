#set -x
SED=gsed
SED_ARGS="-i "

echo "Running sed on $1"

${SED} ${SED_ARGS} -E 's|ScatterMasterToGhosted\(|scatterMasterToGhosted\(|g' "$1"
${SED} ${SED_ARGS} -E 's|GatherGhostedToOwned\(|gatherGhostedToOwned\(|g' "$1"
${SED} ${SED_ARGS} -E 's|PutScalar\(|putScalar\(|g' "$1"

${SED} ${SED_ARGS} -E 's|ViewComponent\((.*)\)->PutScalar|getComponent\(\1\)->putScalar|g' "$1"
${SED} ${SED_ARGS} -E 's|viewComponent\((.*)\)->putScalar|getComponent\(\1\)->putScalar|g' "$1"

${SED} ${SED_ARGS} -E 's|ViewComponent\(|viewComponent\(|g' "$1"
${SED} ${SED_ARGS} -E 's|HasComponent\(|hasComponent\(|g' "$1"
${SED} ${SED_ARGS} -E 's|\.Mesh\(|.getMesh\(|g' "$1"
${SED} ${SED_ARGS} -E 's|->Mesh\(|->getMesh\(|g' "$1"


