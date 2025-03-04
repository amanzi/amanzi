#set -x
SED=sed
SED_ARGS="-i "

echo "Running sed on $1"

${SED} ${SED_ARGS} -E 's|Entity_ID_List|Entity_ID_View|g' "$1"
${SED} ${SED_ARGS} -E 's|Entity_GID_List|Entity_GID_View|g' "$1"
${SED} ${SED_ARGS} -E 's|Entity_Direction_List|Entity_Direction_View|g' "$1"
${SED} ${SED_ARGS} -E 's|Point_List|Point_View|g' "$1"
${SED} ${SED_ARGS} -E 's|Double_List|Double_View|g' "$1"

${SED} ${SED_ARGS} -E 's|MemSpace_type|MemSpace_kind|g' "$1"
${SED} ${SED_ARGS} -E 's|Parallel_type|Parallel_kind|g' "$1"
${SED} ${SED_ARGS} -E 's|Cell_type|Cell_kind|g' "$1"
${SED} ${SED_ARGS} -E 's|Partitioner_type|Partitioner_kind|g' "$1"
${SED} ${SED_ARGS} -E 's|AccessPattern[_kind]*|AccessPattern_kind|g' "$1"

#${SED} ${SED_ARGS} -E 's|std::vector[[:space:]]*<[[:space:]]*std::vector[[:space:]]*<[[:space:]]*.*Entity_ID[[:space:]]*>[[:space:]]*>|std::vector<AmanziMesh::Entity_ID_List>|g' "$1"
#${SED} ${SED_ARGS} -E 's|std::vector[[:space:]]*<[[:space:]]*std::vector[[:space:]]*<[[:space:]]*.*Point[[:space:]]*>[[:space:]]*>|std::vector<AmanziMesh::Point_List>|g' "$1"

#${SED} ${SED_ARGS} -E 's|std::vector[[:space:]]*<[[:space:]]*.*Entity_ID[[:space:]]*>|AmanziMesh::Entity_ID_List|g' "$1"
#${SED} ${SED_ARGS} -E 's|std::vector[[:space:]]*<[[:space:]]*double[[:space:]]*>|AmanziMesh::Double_List|g' "$1"
#${SED} ${SED_ARGS} -E 's|std::vector[[:space:]]*<[[:space:]]*.*Point[[:space:]]*>|AmanziMesh::Point_List|g' "$1"

git -P diff --shortstat $1