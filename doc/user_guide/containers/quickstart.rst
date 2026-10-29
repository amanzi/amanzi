.. _Container-Quickstart:

Quickstart: Running Amanzi or ATS in a Container
================================================

The fastest way to run Amanzi or ATS without building from source is to pull
a prebuilt image from Docker Hub and run your input file inside a container.

Prerequisites
-------------

* A working `Docker <https://docs.docker.com/get-docker/>`_ installation
  (Docker Desktop on macOS/Windows, or Docker Engine on Linux).
* An Amanzi (``.xml``) or ATS input file to run.

The Image Stack
---------------

The Amanzi-ATS Docker images are organized as a three-image stack:

#. ``metsi/amanzi-tpls`` -- a base image with compilers, CMake, and the
   compiled third-party libraries (TPLs).
#. ``metsi/amanzi`` -- a compiled version of Amanzi, built on the base image.
#. ``metsi/ats`` -- a compiled version of ATS, built on the base image.

Images 2 and 3 are rebuilt automatically by GitHub Actions CI and pushed to
Docker Hub.  The ``metsi/amanzi`` image is built by the CI in the
``amanzi/amanzi`` repository, while the ``metsi/ats`` image is built by the CI
in the ``amanzi/ats`` repository.  See :ref:`Container-Image-Tags` for the
tagging scheme.

Pulling an Image
----------------

Pull the most recent image built from ``master``:

.. code-block:: console

   docker pull metsi/amanzi:master-latest
   docker pull metsi/ats:master-latest

Running a Simulation
--------------------

A container has its own isolated filesystem, so to run your own input file you
must *mount* a directory from the host machine into the container.  The
recommended pattern is to run from the top level of your working directory and
mount it into the container's work directory.

To run Amanzi on four MPI ranks:

.. code-block:: console

   HOST_MNT=`pwd -P`
   CONT_MNT=/home/amanzi_user/work

   docker run --rm -v $HOST_MNT:$CONT_MNT:delegated -w $CONT_MNT \
       metsi/amanzi:master-latest mpirun -n 4 amanzi input.xml

The flags used above are:

* ``--rm`` -- delete the container automatically on exit.
* ``-v host:container`` -- bind-mount a host directory into the container.
* ``-w`` -- set the working directory inside the container.

To run ATS, substitute the ``metsi/ats`` image and the ``ats`` executable.

Helper Scripts
--------------

Rather than remembering the full command line, the ``Docker`` directory in the
Amanzi source tree provides small wrapper scripts for common tasks:

* ``run-amanzi-docker`` -- run Amanzi on an input file in the current directory.
* ``run-ats-docker`` -- start an ATS container.
* ``run-ats-shell-docker`` -- start an interactive ATS shell (see
  :ref:`Container-Development`).
* ``run-ats-jupyter-docker`` -- start Jupyter Lab in an ATS container (see
  :ref:`Container-Jupyter`).

These are simple starting points; copy and customize them for your own
workflow.

.. note::

   The ``metsi/amanzi`` images on Docker Hub are multi-architecture: the
   ``metsi/amanzi:master-latest`` tag is a manifest covering both ``amd64``
   (``x86_64``) and ``arm64`` builds.  Docker automatically selects the build
   matching your machine, so on Apple computers with M-series (``arm64``) chips
   the native ``arm64`` image is pulled and runs without emulation.  See
   :ref:`Container-Multiarch` for more on the multi-architecture build process.
