.. _Container-Multiarch:

Multi-Architecture Builds
=========================

Docker provides several ways to build images for more than one processor
architecture.  The most common reason to need this is to build native
``arm64`` images for newer Apple computers (the M-series chips) in addition to
the ``amd64``/``x86_64`` images used by Intel/AMD chipsets and most HPC
centers.

There are three broad approaches:

* build multiple architectures on a single system through emulation,
* cross-compile on a single system, or
* build on separate machines and combine the results in a Docker manifest.

Emulation on a Single System
----------------------------

This is the easiest approach, but can be slow.  It uses ``buildx``, a
BuildKit-based tool, to build multi-architecture images through emulation.

In practice, pass the ``--multiarch`` flag to the deploy scripts
(``deploy-tpls-docker.sh``, ``deploy-amanzi-docker.sh``, and
``deploy-ats-docker.sh``).  This flag switches the scripts from a plain
``docker build`` for the local architecture to
``docker buildx build --platform=linux/amd64,linux/arm64``, building both
architectures in one step.  It assumes Docker is already configured to build
multi-architecture images.

.. note::

   Compiling the TPLs through emulation takes several hours.  Plan
   accordingly.

Cross-Compilation on a Single System
------------------------------------

Building with cross-compilers avoids the runtime cost of emulation.  This
approach is under development for the Amanzi and ATS containers.

Combining Per-Architecture Builds in a Manifest
-----------------------------------------------

Another option is to build each architecture on a machine of that
architecture, push both, and then combine them under a single tag using
``docker manifest``.

Build and push each architecture:

.. code-block:: console

   # On an amd64 machine
   docker build -f Dockerfile-TPLs -t metsi/amanzi-tpls:<version>-amd64-mpich-ubuntu-jammy .
   docker push metsi/amanzi-tpls:<version>-amd64-mpich-ubuntu-jammy

   # On an arm64 machine
   docker build -f Dockerfile-TPLs -t metsi/amanzi-tpls:<version>-arm64-mpich-ubuntu-jammy .
   docker push metsi/amanzi-tpls:<version>-arm64-mpich-ubuntu-jammy

With both images on Docker Hub, combine them under a common tag:

.. code-block:: console

   docker manifest create \
       metsi/amanzi-tpls:<version>-mpich-ubuntu-jammy \
       --amend metsi/amanzi-tpls:<version>-amd64-mpich-ubuntu-jammy \
       --amend metsi/amanzi-tpls:<version>-arm64-mpich-ubuntu-jammy \
   && docker manifest push metsi/amanzi-tpls:<version>-mpich-ubuntu-jammy

Clients then pull ``metsi/amanzi-tpls:<version>-mpich-ubuntu-jammy`` and
automatically receive the image matching their architecture.
