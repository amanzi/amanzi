.. _Containers:

Containers: Running and Developing Amanzi-ATS with Docker
=========================================================

Through GitHub Actions, the Amanzi-ATS project publishes Docker images of
builds on the ``master`` branch, release branches, and feature branches.
These images can be used to simply run Amanzi or ATS locally in a container,
or they can be used for a range of other purposes, such as developing new
features, and identifying and fixing bugs.

This section documents the common container workflows, the essential Docker
and Git commands they rely on, and best practices for building, running, and
extending the images.

.. _Container-Terminology:

Key Terminology
---------------

As with any technology, there is some jargon that must be kept straight,
otherwise things that are relatively simple can seem confusing.  First, it
is important to distinguish the terms *image* and *container*:

* **Image**: A Docker image includes everything needed to run an application
  in a container: the code, its libraries, environment variables, and so on.
  Images are what the project builds and shares through Docker Hub, such as
  the ``metsi/ats:master-latest`` image.

* **Container**: A Docker container is a running instance of an image.  When
  you launch Jupyter Lab, for example, it runs inside a container that is an
  instance of an image.  Many containers can be started from the same image,
  and you can run multiple shells inside a single running container.

Because images and containers are distinct entities, it is important to
understand which Docker commands act on which entity.

Common Commands
---------------

The workflows in this section draw on a small set of Docker and Git commands.

**Docker**

* ``docker image ls`` -- list local images
* ``docker container ls -a`` -- list containers, including stopped ones
* ``docker run`` -- create and start a container from an image
* ``docker exec`` -- start a new process in a running container
* ``docker rename`` -- rename a container
* ``docker commit`` -- save a container's state as a new image
* ``docker push`` / ``docker pull`` -- share images through a registry

**Git**

* ``git config`` -- set identity and options
* ``git remote`` -- manage remote repositories
* ``git fetch`` / ``git merge`` / ``git pull`` -- update from a remote
* ``git add`` / ``git commit`` -- record changes
* ``git checkout`` / ``git branch`` -- manage branches

Contents
--------

.. toctree::
   :maxdepth: 2

   quickstart.rst
   image_tags.rst
   run_jupyter.rst
   develop_in_container.rst
   build_images.rst
   multiarch.rst
