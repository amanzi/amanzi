#!/usr/bin/env bash

# Generic ATS build script that invokes Amanzi's bootstrap.sh 
# with options configured for an ATS build

# see INSTALL_ATS.md for more information

if [ ${ATS_BUILD_TYPE} == Debug ]
then
    dbg_option=--debug
elif [ ${ATS_BUILD_TYPE} == RelWithDebInfo ]
then
    dbg_option=--relwithdebinfo
else
    dbg_option=--opt
fi

echo "Running Amanzi Boostrap with:"
echo "  AMANZI_SRC_DIR: ${AMANZI_SRC_DIR}"
echo "  AMANZI_BUILD_DIR: ${AMANZI_BUILD_DIR}"
echo "  AMANZI_INSTALL_DIR: ${AMANZI_DIR}"
echo "  AMANZI_TPLS_BUILD_DIR: ${AMANZI_TPLS_BUILD_DIR}"
echo "  AMANZI_TPLS_DIR: ${AMANZI_TPLS_DIR}"
echo " with build type: ${ATS_BUILD_TYPE}"
echo ""


${AMANZI_SRC_DIR}/bootstrap.sh \
   ${dbg_option} \
   --with-mpi=${OPENMPI_DIR} \
   --enable-shared \
   --enable-clm \
   --disable-structured  --enable-unstructured \
   --disable-stk_mesh --enable-mstk_mesh \
   --enable-hypre \
   --enable-silo \
   --disable-petsc \
   --disable-amanzi_physics \
   --enable-ats_physics \
   --enable-ats_dev \
   --disable-geochemistry \
   --amanzi-install-prefix=${AMANZI_DIR} \
   --amanzi-build-dir=${AMANZI_BUILD_DIR} \
   --tpl-install-prefix=${AMANZI_TPLS_DIR} \
   --tpl-build-dir=${AMANZI_TPLS_BUILD_DIR} \
   --tpl-download-dir=${ATS_BASE}/amanzi-tpls/Downloads \
   --tools-download-dir=${ATS_BASE}/amanzi-tpls/Downloads \
   --tools-build-dir=${ATS_BASE}/build \
   --tools-install-prefix=${ATS_BASE}/install \
   --with-cmake=`which cmake` \
   --with-ctest=`which ctest` \
   --branch_ats=${ATS_VERSION} \
   --parallel=6


# If TPLs have already been built, and you don't want to go 
# through that long process again, replace
#   --tpl-install-prefix=${AMANZI_TPLS_DIR} \
#   --tpl-build-dir=${AMANZI_TPLS_BUILD_DIR} \
#   --tpl-download-dir=${ATS_BASE}/amanzi-tpls/Downloads \
# with
#   --tpl-config-file=${AMANZI_TPLS_DIR}/share/cmake/amanzi-tpl-config.cmake \


