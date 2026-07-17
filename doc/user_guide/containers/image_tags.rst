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

Branch names containing a ``/`` are sanitized when forming the tag: each ``/``
is replaced with ``--``.  For example, a branch named ``user/feat`` produces
the tag ``user--feat-latest``.

Only the ``latest`` tag is currently retained for each branch.  The TPLs image
is additionally tagged with a version number when a new set of libraries is
released (see :ref:`Container-Build-Images`).

Architectures
-------------

The ``metsi/amanzi`` image is published as a multi-architecture manifest that
covers both the ``amd64``/``x86_64`` architecture (used by Intel/AMD chipsets
and most HPC centers) and the ``arm64`` architecture (used by Apple computers
with M-series chips).  CI builds each architecture natively on its own runner
and then combines them into a single manifest.

Pulling the branch tag, for example ``metsi/amanzi:master-latest`` or
``metsi/amanzi:<branch>-latest``, automatically selects the image matching the
host architecture.  On Apple Silicon this means you get a native ``arm64``
image and run without emulation.

Per-architecture tags are also published if you need to pull a specific
architecture explicitly:

* ``metsi/amanzi:<branch>-amd64-latest`` -- the ``amd64`` build.
* ``metsi/amanzi:<branch>-arm64-latest`` -- the ``arm64`` build.

See :ref:`Container-Multiarch` for the available approaches to
multi-architecture builds.
