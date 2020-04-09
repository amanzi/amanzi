Spack usage.

# Add repository
spack repo add amanzi/spack-repo

# Using environement
spack env create amanzi
spack env activate amanzi

# Install
spack install amanzi@0.89 +mstk

# Install in the current directory, for development build
spack dev-build amanzi@0.89 +mstk

# Options
TODO: add description of options