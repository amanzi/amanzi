.. _Container-Development:

Developing in a Container
=========================

Beyond simply running simulations, containers provide a self-contained,
reproducible environment for developing new features, fixing bugs, and adding
system software.  This page covers three common developer workflows:

* installing additional system software in a running container,
* saving that state as a new image, and
* contributing a bug fix or feature back through a fork and pull request.

Starting an Interactive Shell
-----------------------------

Most development starts from an interactive shell.  The
``run-ats-shell-docker`` helper script sets up a suitable ``docker run``
command:

.. code-block:: console

   run-ats-shell-docker

which is equivalent to:

.. code-block:: console

   HOST_MNT=`pwd -P`
   CONT_MNT=/home/amanzi_user/work

   docker run -it --mount type=bind,source=$HOST_MNT,target=$CONT_MNT \
       -w $CONT_MNT metsi/ats:master-latest /bin/bash

Once the container starts, you are the generic, unprivileged user
``amanzi_user``.

Adding System Software as Root
------------------------------

Installing system packages requires root privileges.  Unlike a normal system,
containers do not provide ``sudo``.  Instead, use ``docker exec`` to start a
new shell as root (user id ``0``) in the *already running* container.  First
give the container a memorable name, then exec into it from another host
terminal:

.. code-block:: console

   docker rename <funky_name> ideas_ats_container
   docker exec -it -u 0 ideas_ats_container /bin/bash

You now have a root shell:

.. code-block:: console

   root@bda2ca690491:/home/amanzi_user/work# whoami
   root

The base image is built from Ubuntu using ``apt`` in non-interactive mode.
When working interactively it helps to install a dialog utility first:

.. code-block:: console

   apt-get update
   apt-get install dialog whiptail

You can then install any system tool you need, for example the X11 sample
applications:

.. code-block:: console

   apt-get install x11-apps

Using X11 Applications (macOS)
------------------------------

To display X11 applications from the container on a macOS host:

#. Install the latest `XQuartz <https://www.xquartz.org/>`_ X11 server and
   run it.
#. In XQuartz settings, enable *Allow connections from network clients*.
#. Quit and restart XQuartz to activate the setting.

In the XQuartz terminal, allow connections from localhost:

.. code-block:: console

   xhost + 127.0.0.1

.. warning::

   Allowing X11 connections is insecure for remote hosts.  Restrict it to
   ``127.0.0.1`` (localhost) as shown above.

Then, inside the container, point the display at the host and run an
application:

.. code-block:: console

   export DISPLAY=host.docker.internal:0
   xeyes

Saving Changes as a New Image
-----------------------------

Changes made inside a running container are lost when the container is
removed.  To preserve them, commit the container to a new image:

.. code-block:: console

   docker commit ideas_ats_container <user>/amanzi:ats-xeyes

Verify the new image by running a shell from it:

.. code-block:: console

   docker run -it --rm <user>/amanzi:ats-xeyes /bin/bash

Contributing a Bug Fix or Feature
---------------------------------

You can develop and submit changes entirely from within a container.  Start an
interactive shell as above, then configure Git so your commits record your
identity:

.. code-block:: console

   git config --global user.name "Your Name"
   git config --global user.email "your.name@example.com"
   git config --list | grep user

Work on a purpose-named branch rather than committing directly to ``master``:

.. code-block:: console

   git checkout -b <yourname>/short-description
   git branch --list

Stage and commit your change:

.. code-block:: console

   git add <files>
   git commit -m "Describe the change."

To share the change, first create a fork on ``github.com`` using the web
interface.  Inside the container, rename the existing ``origin`` remote to
``upstream`` and add your fork as ``origin``:

.. code-block:: console

   git remote rename origin upstream
   git remote add origin https://github.com/<youruser>/<yourfork>
   git remote -v

Push your branch to the fork:

.. code-block:: console

   git push origin <yourname>/short-description

Finally, open a pull request from your fork against the upstream repository
using the GitHub web interface.
