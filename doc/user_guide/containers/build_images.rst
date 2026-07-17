.. _Container-Build-Images:

Building the Images Locally
===========================

Most users can rely on the prebuilt images from Docker Hub (see
:ref:`Container-Quickstart`).  You may want to build the images yourself when
you need a native ``arm64`` image, are testing a branch, or are updating the
third-party libraries.

The Dockerfiles and helper scripts live in the ``Docker`` directory of the
Amanzi source tree.  The stack is built in the order TPLs, then Amanzi, then
ATS, because each image builds on the previous one.

.. note::

   Building the TPLs from scratch is time-consuming, and doing so through
   emulation on a different architecture can take several hours.

Cleaning the Build Environment
------------------------------

To remove old build artifacts and start from a clean state:

.. code-block:: console

   cd amanzi/Docker
   docker system prune -a

.. warning::

   ``docker system prune -a`` deletes all stopped containers and all unused
   images on your system.  Make sure you do not need any of them before
   running it.

Building the TPLs Image
-----------------------

Update the paths in ``deploy-tpls-docker.sh`` to match your system, then run
it:

.. code-block:: console

   cd amanzi/Docker
   ./deploy-tpls-docker.sh

The script tags the resulting image
``metsi/amanzi-tpls:<version>-<mpi>-<base>-<ver_tag>`` (for example
``metsi/amanzi-tpls:0.98.9-mpich-ubuntu-jammy``), where ``<version>`` is read
from ``config/SuperBuild/TPLVersions.cmake``.  If you are not pushing the image
to the ``metsi/amanzi-tpls`` Docker Hub repository, this is as far as you need
to go.  To publish it (this requires access to the repository), re-run the
script with ``--push``:

.. code-block:: console

   ./deploy-tpls-docker.sh --push

.. warning::

   Pushing overwrites any image already in the repository that has the same
   tag.

Building the Amanzi Image
-------------------------

``Dockerfile-Amanzi`` builds Amanzi on top of a ``metsi/amanzi-tpls`` image.  It
selects its base image from the ``amanzi_tpls_ver`` and ``mpi_flavor`` build
arguments as ``metsi/amanzi-tpls:${amanzi_tpls_ver}-${mpi_flavor}`` (for example
``metsi/amanzi-tpls:0.98.9-mpich``).  Rather than invoking ``docker build``
directly, use the helper script, which fills in these arguments and tags the
result for you:

.. code-block:: console

   cd amanzi/Docker
   ./deploy-amanzi-docker.sh

The script derives the Amanzi version from the most recent ``amanzi-*`` git tag
and the current commit hash, and tags the resulting image both
``metsi/amanzi:<version>`` and ``metsi/amanzi:latest``.  It uses the ``mpich``
TPLs image by default; pass ``--mpi_flavor`` to change this, and
``--amanzi_tpls_ver`` to build against a different TPLs version.  If a matching
local ``metsi/amanzi-tpls`` image already exists, the build uses it instead of
downloading from Docker Hub, which is convenient for testing an Amanzi build
against TPLs built on a non-``master`` branch.

If your system is behind a proxy, export the ``http_proxy`` and ``https_proxy``
environment variables and pass the script's ``--use_proxy`` option to forward
them to the build.

To publish (requires repository access), add ``--push``:

.. code-block:: console

   ./deploy-amanzi-docker.sh --push

The ATS image is built in a similar way with ``Dockerfile-ATS-build`` (or the
``deploy-ats-docker.sh`` helper).  Note that, unlike ``Dockerfile-Amanzi`` --
which bases on the two-part ``<version>-<mpi>`` TPLs tag --
``Dockerfile-ATS-build`` bases on the full four-part
``<version>-<mpi>-<base>-<ver_tag>`` TPLs tag (for example
``metsi/amanzi-tpls:0.98.9-mpich-ubuntu-jammy``).

Building Amanzi Interactively for Debugging
-------------------------------------------

To debug a build, start an interactive container from the TPLs image and build
Amanzi by hand:

.. code-block:: console

   docker run -it metsi/amanzi-tpls:<tag>

You start in ``/home/amanzi_user/amanzi``.  Build using the same
``bootstrap.sh`` arguments as the final step of ``Dockerfile-Amanzi``:

.. code-block:: console

   ./bootstrap.sh \
       --prefix=/home/amanzi_user/install \
       --amanzi-build-dir=/home/amanzi_user/amanzi_builddir/amanzi \
       --tpl-config-file=/home/amanzi_user/install/tpls/share/cmake/amanzi-tpl-config.cmake \
       --parallel=4 --opt \
       --enable-structured \
       --disable-build_user_guide \
       --enable-alquimia --enable-pflotran --enable-crunchtope \
       --with-mpi=/usr \
       --enable-shared

Type ``exit`` (or press ``Ctrl-C``) to leave the container.
