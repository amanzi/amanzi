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

If you are not pushing the image to the ``metsi/amanzi-tpls`` Docker Hub
repository, this is as far as you need to go.  To publish it, tag it with a
version number and push (this requires access to the repository):

.. code-block:: console

   docker tag metsi/amanzi-tpls:latest metsi/amanzi-tpls:<version>
   docker push metsi/amanzi-tpls

.. warning::

   Pushing overwrites any image already in the repository that has the same
   tag.

Building the Amanzi Image
-------------------------

``Dockerfile-Amanzi`` pulls whatever ``metsi/amanzi-tpls`` image is tagged
``latest`` and builds Amanzi on top of it.  To build against an older TPLs
version, change the tag on the first line of ``Dockerfile-Amanzi`` from
``latest`` to the desired version number (and reset it to ``latest`` before
merging to ``master``).

If a local image tagged ``metsi/amanzi-tpls:latest`` already exists, the build
uses it instead of downloading from Docker Hub.  This is convenient for
testing an Amanzi build against TPLs built on a non-``master`` branch.

Build the image:

.. code-block:: console

   cd amanzi/Docker
   docker build -f Dockerfile-Amanzi -t metsi/amanzi:latest .

The flags are:

* ``-f`` -- the Dockerfile containing the build instructions.
* ``-t`` -- the tag for the resulting image (otherwise it is named by hash).

If your system uses a proxy, pass it as build arguments:

.. code-block:: console

   docker build --build-arg http_proxy=<proxy:port> \
       --build-arg https_proxy=<proxy:port> \
       -f Dockerfile-Amanzi -t metsi/amanzi:latest .

To publish (requires repository access):

.. code-block:: console

   docker tag metsi/amanzi:latest metsi/amanzi:<version>
   docker push metsi/amanzi

The ATS image is built the same way, using ``Dockerfile-ATS-build``.

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
       --prefix=/home/amanzi_user/local \
       --amanzi-build-dir=/home/amanzi_user/amanzi_builddir/amanzi \
       --tpl-config-file=/home/amanzi_user/local/tpls/share/cmake/amanzi-tpl-config.cmake \
       --with-mpi=/usr \
       --enable-shared

Type ``exit`` (or press ``Ctrl-C``) to leave the container.
