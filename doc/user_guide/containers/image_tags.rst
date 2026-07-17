.. _Container-Image-Tags:

Images, Tags, and Architectures
===============================

This page describes what images are published, how they are tagged, and which
processor architectures they support.

Published Images
----------------

The Amanzi-ATS image stack consists of three images on Docker Hub:

* ``metsi/amanzi-tpls`` -- base image with compilers, CMake, and the compiled
  third-party libraries.
* ``metsi/amanzi`` -- a compiled version of Amanzi.
* ``metsi/ats`` -- a compiled version of ATS.

The ``metsi/amanzi`` and ``metsi/ats`` images are rebuilt by GitHub Actions CI
on each push and pull request in the ``amanzi/amanzi`` and ``amanzi/ats``
repositories.  The ``metsi/amanzi-tpls`` base image is rebuilt only
occasionally, usually when a TPL version is updated.

Tagging Scheme
--------------

Images are tagged by the branch they were built from, using the ``latest``
image built from the most recent commit on that branch.  For example:

* ``metsi/ats:master-latest`` -- most recent build of the ``master`` branch.
* ``metsi/amanzi:master-latest`` -- most recent build of ``master``.
* ``metsi/amanzi:<branch>-latest`` -- most recent build of a feature or
  release branch.

Only the ``latest`` tag is currently retained for each branch.  The TPLs image
is additionally tagged with a version number when a new set of libraries is
released (see :ref:`Container-Build-Images`).

Architectures
-------------

The images on Docker Hub are currently built only for the ``amd64``/``x86_64``
architecture, which is used by Intel/AMD chipsets and most HPC centers.

On Apple computers with M-series chips (``arm64`` architecture), these images
run through emulation and can be very slow, because the code must be
translated at runtime.  Users with M-series chips who need better performance
may want to build native ``arm64`` images locally.  See
:ref:`Container-Multiarch` for the available approaches to multi-architecture
builds.
