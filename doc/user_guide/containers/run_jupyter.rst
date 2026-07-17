.. _Container-Jupyter:

Running Jupyter Lab in a Container
==================================

Some Amanzi-ATS workflows -- particularly post-processing and short-course
notebooks -- are driven from Jupyter Lab.  This page describes how to run
Jupyter Lab from within a container and reach it from your host web browser.

The essential difference from the basic run command (see
:ref:`Container-Quickstart`) is that a network port must be mapped from the
container to the host so that the browser can connect to the Jupyter server.

Running a Jupyter-Enabled Image
-------------------------------

Assuming you have an image with Jupyter Lab installed (see
:ref:`Installing Jupyter Lab <Container-Install-Jupyter>` below to build one),
start it with a port mapping:

.. code-block:: console

   HOST_MNT=`pwd -P`
   CONT_MNT=/home/amanzi_user/work

   docker run --rm -it --init \
       --mount type=bind,source=$HOST_MNT,target=$CONT_MNT -w $CONT_MNT \
       -p 8899:8899 \
       <user>/amanzi:ats-jupyter \
       /bin/bash -c 'conda activate base; jupyter lab --port 8899'

The additional flags are:

* ``--init`` -- use the default ``tini`` init process inside the container,
  so that Jupyter and its child processes are reaped cleanly.
* ``-p 8899:8899`` -- map container port 8899 to host port 8899.

Once the server starts, it prints a URL with an access token.  Open that URL
(or ``http://localhost:8899``) in your host browser.

The ``run-ats-jupyter-docker`` helper script in the ``Docker`` directory wraps
this command.

.. _Container-Install-Jupyter:

Installing Jupyter Lab in a Container
-------------------------------------

Installing Jupyter Lab inside a container is very similar to installing it as
part of Anaconda on your laptop.  The steps below install it as the
unprivileged ``amanzi_user``.

Start an interactive ATS shell:

.. code-block:: console

   run-ats-shell-docker

Download the Anaconda shell installer on the host and place it in a mounted
folder so it is visible inside the container, for example:

.. code-block:: console

   mv ~/Downloads/Anaconda3-2023.07-2-Linux-x86_64.sh <mounted-downloads-dir>

Inside the container, create a place for the scripts and installation, then
run the installer in batch mode:

.. code-block:: console

   mkdir anaconda && cd anaconda
   mkdir scripts logs

   # -b batch mode (accepts the license), -p installation prefix
   /bin/bash <downloads>/Anaconda3-2023.07-2-Linux-x86_64.sh \
       -b -p $HOME/anaconda/Anaconda3-2023.07-2

Add the conda initialization to ``~/.bashrc`` so ``conda`` is available in new
shells, then configure Jupyter.  The most important settings allow remote
connections (and, if needed, running as root):

.. code-block:: python

   c.ServerApp.allow_root = True
   c.ServerApp.allow_remote_access = True

The ``container-anaconda.sh`` and ``container-jupyter-setup.sh`` scripts in the
``Docker`` directory automate these installation and configuration steps.

Saving the Result as an Image
-----------------------------

Once Jupyter Lab is installed and configured, commit the container to a new
image so the setup can be reused (see :ref:`Container-Development` for more on
``docker commit``):

.. code-block:: console

   docker rename <container-name> ats_jupyter
   docker commit ats_jupyter <user>/amanzi:ats-jupyter

.. note::

   A full Anaconda installation makes the image large (roughly 10 GB).  Using
   Miniconda and installing only the packages you need can reduce this
   substantially, though the image will still be several GB.
