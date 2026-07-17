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

In practice, replace ``build`` with ``buildx build --platform=<platforms>`` in
the deploy scripts, for example ``--platform=linux/amd64,linux/arm64``.

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
   docker build -t metsi/amanzi-tpls:mpich-<version>-amd64 .
   docker push metsi/amanzi-tpls:mpich-<version>-amd64

   # On an arm64 machine
   docker build -t metsi/amanzi-tpls:mpich-<version>-arm64 .
   docker push metsi/amanzi-tpls:mpich-<version>-arm64

With both images on Docker Hub, combine them under a common tag:

.. code-block:: console

   docker manifest create \
       metsi/amanzi-tpls:mpich-<version> \
       --amend metsi/amanzi-tpls:mpich-<version>-amd64 \
       --amend metsi/amanzi-tpls:mpich-<version>-arm64 \
   && docker manifest push metsi/amanzi-tpls:mpich-<version>

Clients then pull ``metsi/amanzi-tpls:mpich-<version>`` and automatically
receive the image matching their architecture.
