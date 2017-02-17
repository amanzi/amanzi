#! /bin/bash

spack_binary=$1
spack_pkg=$2

${spack_binary} install ${spack_pkg}
${spack_binary} find -pe ${spack_pkg} >> ${spack_pkg}_install_path.in
sed -i -e 's#.*${spack_pkg}@.*  ##g' ${spack_pkg}_install_path.in
pwd
 
