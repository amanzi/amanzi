#! /bin/bash

# grab the arguments
spack_binary=$1
spack_pkg=$2
install_dir=$3

# install the package
${spack_binary} install ${spack_pkg}

# create symlinks between where spack installs the pkg
# and where you really want to install it
${spack_binary} view symlink ${install_dir} ${spack_pkg}


#${spack_binary} location -i ${spack_pkg} >> ${spack_pkg}_install_path.in
#cat ${spack_pkg}_install_path.in | $3

 
