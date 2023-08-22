# Docker Workflows in IDEAS-Watersheds 

## Key Terminology

As with any technology there is some jargon that we need to keep straight otherwise things that are relatively simple will seem very confusing.  First, it is important to define the terms image and container:

  * _Image_: A docker image includes everything that is needed to run an application in a docker container, such as, code, libraries, environment variables, etc.  Note images are what we build and share with you through docker hub, such as the metsi/ats-short-course:latest image you pulled down to  your system. 
  
  * _Container_: A docker container is a runtime instance of a docker image.  So when you run jupyterlab it is running in a container that is an instance of the docker image metsi/ats-short-course:latest.  So you could have many containers started from the same image, you can run multiple shells 

Since images and containers are different entities, it's important to understand which docker commands interact with which entity.  

## Common Commands

### Docker
* docker image ls 
* docker container ls -a

* docker run
* docker exec

* docker rename
* docker commit

* docker push
* docker pull

### Git

* git config
* git remote

* git fetch
* git merge
* git pull

* git add
* git commit

* git checkout
* git branch


## Add (system) software to a running container

Run an ATS container with an interactive shell (i.e., bash).  I recommend creating small scripts for this type of task, rather than trying to remember all the command line options you need and/or want.  So I would execute 

```
run-ats-shell-docker
```

which sets up a simple docker run command

```sh
#
#  --mount provides the mount type,source, and target, e.g.,
#          --mount type=bind,source=mount_point_on_host,target=mount_point_on_container
#  -w      sets current working directory on the container
#  -it     -i and -t together provide an interactive terminal (a tty) for the container
#  -p      map host_port:container_port, e.g., -p 8899:8899 
#  --rm    deletes the container on exit
#
#  On the host, we recommend that the current pwd is the top-level of your desired workspace.
#
HOST_MNT=`pwd -P`

#  On the container, we recommend your mount point is the short-course directory and that
#  is also your PWD on startup.
#
CONT_MNT=/home/amanzi_user/work
CONT_PWD=/home/amanzi_user/work

docker run -it --mount type=bind,source=$HOST_MNT,target=$CONT_MNT -w $CONT_PWD metsi/ats:master-latest /bin/bash
```

Once the container starts you are presented with the bash command line as the generic user named "amanzi_user".  However, this is an unprivileged user and you need to be root (the most privileged user on unix systems) to install additional packages using the system level tools.  Unlike regular systems, containers do not allow sudo (i.e., temporarily elevating privileges of the current shell/user).  Instead they provide an exec command to start new processes in a running container.  With exec you can start a new shell as root,

```
docker rename <funky_name>  ideas_ats_container
docker exec -it -u 0 ideas_ats_container /bin/bash
```

Now we have a root shell and we are ready to go

```console
root@bda2ca690491:/home/amanzi_user/work# whoami
root
```

The Ubuntu image was created using apt in a non-interactive mode.  Now if you are working on the command line its best to install a dialog utility first as this enables reporting and interaction with the terminal.

```
apt-get install dialog
apt-get install whiptail
```

Lets install something fun, some basic X11 applications that will highlight how to use X-windows if necessary.  First on the Mac side you need to enable  

 * Install the latest XQuartz X11 server and run it
 * Activate the option ‘Allow connections from network clients’ in XQuartz settings
 * Quit & restart XQuartz (to activate the setting)

In the xterm that opens when you restart XQuartz, allow connections from the localhost by entering:

```
xhost + 127.0.0.1
```
 
This isn't a very secure thing to do for remote connetions, but for localhost it's ok. Now let's install an Ubuntu package

```
apt-get install x11-apps
```

Set the display variable:

```
export DISPLAY=host.docker.internal:0
```

and run the application

```
xeyes
```

You could install xterm, or other python packages, or whatever system level tool you find you're missing for your debugging.

```
docker commit ideas_ats_container moulton/amanzi:ats-xeyes-v2
```

## Contributing a bug fix or feature

Run an ATS container with an interactive shell (i.e., bash), as we did in the first example,

```
run-ats-shell-docker
```

Once it starts you are presented with the bash command line as the generic user named "amanzi_user".  In order to use git, however, you'll want your commit messages to record your name and email.  So the first step is to update the git configuration

```sh
git config --global user.name "David Moulton"
git config --global user.email "moulton@lanl.gov"
```

Confirm you got it right ...

```sh
git config --list | grep user
```

```console
amanzi_user@62d5ae7ce087:~/amanzi/src$ git config --list | grep user
user.name=David Moulton
user.email=moulton@lanl.gov
```

Often with fork you keep changes on the master branch, but generally I prefer to work in branches with specific a purpose.  So let's make a branch

```sh
git checkout -b david/container-pr-demo
git branch --list
```

Then let's add and commit our change on that branch

```sh
git add README.md
git commit -m "Added enlightening info to the README."
```

Now the next step is to figure out how to share it.  First we fork using the web interface on github.com. In this case I created a fork called ats-containers under my personal account (jd-moulton).

```console
amanzi_user@62d5ae7ce087:~/amanzi/src/physics/ats$ git remote -v
origin	https://github.com/amanzi/ats (fetch)
origin	https://github.com/amanzi/ats (push)
amanzi_user@62d5ae7ce087:~/amanzi/src/physics/ats$ git remote rename origin upstream
```

```console
amanzi_user@62d5ae7ce087:~/amanzi/src/physics/ats$ git remote -v
upstream	https://github.com/amanzi/ats (fetch)
upstream	https://github.com/amanzi/ats (push)
```

```sh
git remote add origin https://github.com/jd-moulton/ats-containers
```

```console
amanzi_user@62d5ae7ce087:~/amanzi/src/physics/ats$ git remote -v
origin	https://github.com/jd-moulton/ats-containers (fetch)
origin	https://github.com/jd-moulton/ats-containers (push)
upstream	https://github.com/amanzi/ats (fetch)
upstream	https://github.com/amanzi/ats (push)
```

Push your update to the forked repository ...

```sh
amanzi_user@62d5ae7ce087:~/amanzi/src/physics/ats$ git push origin david/container-pr-demo
Username for 'https://github.com': jd-moulton
Password for 'https://jd-moulton@github.com':
```

You can go to the github.com fork and see the update is there, and you can easily add a pull request.


## Install Jupyter Lab as an unprivileged user

This is very similar to how you would install Jupyter Lab as part of Anaconda on your laptop natively.  And again, rather than try to remember details and command line options, I use simple scripts to do the installation and setup.
The starting point is as above, run an ATS container with an interactive shell (i.e., bash),

```
run-ats-shell-docker
```

Then download the Anaconda shell installer, and save it locally in a mounted folder. For example,

```sh
mv ~/Downloads/Anaconda3-2023.07-2-Linux-x86_64.sh /amanzi/docker/container/downloads
```

where /amanzi/docker is mounted as /home/amanzi_user/work in the container. Next create a place for the scripts and the installation:

```
mkdir anaconda; cd anaconda
mkdir scripts logs
cp -p ~/work/container/container-anaconda.sh ./scripts
```

Run the installation script:

```
./scripts/container-anaconda.sh
```

which just sets up the batch mode installation command:

```
# Run the sh installer with embedded contents
#
# -b batch mode, assumes you agree to the license
# -p prefix for the installation
#
/bin/bash $DOWNLOADS_DIR/$ACONDA_FILE_NAME -b -p $INSTALL_DIR/$ACONDA_NAME > ${LOGI} 2>&1
```

It also adds conda environment initialization to the ~/.bashrc file.  

Lastly we need to provide the default configuration for jupyter in the ~/.jupyter/jupyter_notebook_config.py. Again, I use a script for this, so we'll copy it over first

```
cd ~/anaconda
cp -p ~/work/container/container-jupyter-setup.sh ./scripts
```

and then run it,

```
./scripts/container-jupyter-setup.sh
```

Note, the most important option is to allow remote connections (and possibly allowing root access):

```
c.NotebookApp.allow_root = True" 
c.NotebookApp.allow_remote_access = True"
```

Finally we commit this container to an image so we can run it in the future. First rename it (probably unecessary, but seems cleaner) and then commit it

```
docker rename practical_raman ideas_ats_jupyter_v3
docker commit ideas_ats_jupyter_v3 moulton/amanzi:ats-jupyter-v3
```

Note I have tagged it for my personal repository on dockerhub. Also, it's huge (~10GB), so not very practical to be moving around.  Using miniconda and installing just what you need could bring this down to 5-7GB, which is still large but better.


To test the new image, run the script

```
run-ats-jupyter-docker
```

which sets up the docker run command similar to before, but adds the port mapping needed for Jupyter Lab to be accessed through the web browser:

```
#  --mount provides the mount type,source, and target, e.g.,
#          --mount type=bind,source=mount_point_on_host,target=mount_point_on_container
#  --init  uses the default tini as the "init" for the container
#  -w      sets current working directory on the container
#  -it     -i and -t together provide an interactive terminal (a tty) for the container
#  -p      map host_port:container_port, e.g., -p 8899:8899 
#  --rm    deletes the container on exit
#

docker run --rm -it --init --mount type=bind,source=$HOST_MNT,target=$CONT_MNT -w $CONT_PWD -p 8899:8899 moulton/amanzi:ats-jupyter-v3 /bin/bash -c 'source /home/amanzi_user/anaconda/Anaconda3-2023.07-2/etc/profile.d/conda.sh; conda activate base; jupyter lab --port 8899'
```

The explicit initialization of the conda environment shouldn't be necessary here as we added this to the ~/.bashrc. But it is, and at least it works.

