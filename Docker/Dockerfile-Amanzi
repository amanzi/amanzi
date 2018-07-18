FROM metsi/amanzi-tpls:latest
LABEL Description="This image builds Amanzi from the latest version of the Amanzi TPLs contained in a Docker image"

# Switch to amanzi_user
USER amanzi_user

# Change the Working Directory and update amanzi
WORKDIR /home/amanzi_user/amanzi
RUN git pull
  
RUN ./bootstrap.sh --prefix=${AMANZI_PREFIX} \
   --amanzi-build-dir=/home/amanzi_user/amanzi_builddir/amanzi \
   --tpl-config-file=${AMANZI_TPLS_DIR}/share/cmake/amanzi-tpl-config.cmake \
   --parallel=4 --opt \
   --disable-structured \
   --enable-alquimia --enable-pflotran --enable-crunchtope \
   --with-mpi=/usr \
   --enable-shared \
   && cd /home/amanzi_user/amanzi_builddir/amanzi \
   && make test \
   && rm -r /home/amanzi_user/amanzi_builddir/
