#!/usr/bin/env bash

rm -rf $ATS_DIR
rm -rf $ATS_BUILD_DIR
mkdir -p $ATS_DIR
mkdir -p $ATS_BUILD_DIR
cd $ATS_BUILD_DIR

cmake \
    -D Amanzi_DIR=${AMANZI_DIR}/lib \
    -D CMAKE_INSTALL_PREFIX=${ATS_DIR} \
    -D CMAKE_BUILD_TYPE=${ATS_BUILD_TYPE} \
  ${ATS_SRC_DIR}

make -j10 install




