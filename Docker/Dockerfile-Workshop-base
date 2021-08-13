#FROM continuumio/miniconda3
FROM ubuntu:20.04

LABEL Description="Amanzi/ATS/WatershedWorkflow/TINerator Docker container"

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 \
    https_proxy=$https_proxy \
    http_proxy=$http_proxy \
    DIR_BUILD=${HOME}/.build \
    DIR_INSTALL=${HOME}/.install

ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y \

    # Shared
    apt-utils \
    patch \
    g++ \
    gfortran \
    git \
    make \
    cmake \
    emacs \
    vim \
    curl \
    wget \
    libz-dev \
    openssl \
    m4 \
    bzip2 \
    ca-certificates \

    # Amanzi/ATS
    groff \
    mpich \
    libblas-dev \
    liblapack-dev \
    libmpich-dev \
    libssl-dev \
    zlib1g-dev \
    libcurl4-openssl-dev \

    # TINerator/Watershed Workflow
    libgdal-dev \
    unzip \
    bison \
    libgl1-mesa-glx \
    xvfb && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# =================================================== #
# Build Miniconda =================================== #
# =================================================== #
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

ENTRYPOINT [ "/usr/bin/tini", "--" ]

RUN . /opt/conda/etc/profile.d/conda.sh && conda activate base

# =================================================== #
# Build Watershed Workflow ========================== #
# =================================================== #
WORKDIR ${DIR_BUILD}
RUN pip install GDAL==`gdal-config --version` && \
    conda install -c conda-forge boost && \
    pip --no-cache-dir install --upgrade \
       fiona \
       rasterio \
       shapely \
       pyproj \
       && \
    git clone https://github.com/ecoon/meshpy.git --depth 1 && \
    cd meshpy && \
    ./configure.py \
        --python-exe=${CONDA_DIR}/bin/python \
        --boost-inc-dir=${CONDA_DIR}/include \
        --boost-lib-dir=${CONDA_DIR}/lib \
        --boost-python-libname=boost_python39 \
        --disable-static \
        --enable-shared \
        && \
    python setup.py build && \
    python setup.py install && \
    cd .. && \

    # IDEAL-OCTO-WAFFLE
    git clone https://github.com/ecoon/ideal-octo-waffle.git --depth 1

ENV PYTHONPATH=${DIR_BUILD}:${PYTHONPATH}

# =================================================== #
# Build TINerator =================================== #
# =================================================== #
RUN git clone https://github.com/lanl/tinerator.git ${DIR_BUILD}/tinerator --depth 1 && \
    cd ${DIR_BUILD}/tinerator && \
    pip install -r requirements.txt && \
    ./tpls/build-tpls.sh -A -M

ENV PYTHONPATH=${DIR_BUILD}/tpls/seacas/install/lib:${PYTHONPATH}

# =================================================== #
# Build Amanzi/ATS ================================== #
# =================================================== #
ENV AMANZI_PREFIX=${DIR_INSTALL} \
    AMANZI_INSTALL_DIR=${DIR_INSTALL}/amanzi \
    AMANZI_TPLS_DIR=${DIR_INSTALL}/tpls \
    AMANZI_SRC_DIR=${DIR_BUILD}/amanzi \
    ATS_SRC_DIR=${DIR_BUILD}/amanzi/src/physics/ats
    #AMANZI_PETSC_LIBS=$AMANZI_TPLS_DIR/petsc-${petsc_ver}/lib \
    #AMANZI_TRILINOS_LIBS=$AMANZI_TPLS_DIR/trilinos-${trilinos_ver}/lib \
    #AMANZI_SEACAS_LIBS=$AMANZI_TPLS_DIR/SEACAS/lib

RUN git clone https://github.com/amanzi/amanzi.git ${DIR_BUILD}/amanzi --depth 1 && \
    cd ${DIR_BUILD}/amanzi/
RUN ./bootstrap.sh --prefix=${AMANZI_PREFIX} \
   --parallel=4 \
   # --opt \
   --amanzi-build-dir=${DIR_BUILD}/amanzi \
   --tpl-build-dir=${DIR_BUILD}/tpls \
   --tpl-download-dir=${DIR_BUILD}/tpls/Downloads \
   # --disable-build_amanzi \
   # --disable-build_user_guide \
   # --enable-shared \
   # --enable-structured \
   # --enable-unstructured \
   # --enable-silo \
   # --enable-alquimia \
   # --enable-pflotran \
   # --enable-crunchtope \
   # --enable-pf_clm \
   --with-mpi=/usr \
   --with-c-compiler=$(which mpicc) \
   --with-cxx-compiler=$(which mpicxx) \
   --with-fort-compiler=$(which mpif90) \
   && git checkout master \
   && rm -r ${DIR_BUILD}

ENV PATH=${DIR_INSTALL}/install/tools/bin:${PATH}
CMD ["/bin/bash"]
