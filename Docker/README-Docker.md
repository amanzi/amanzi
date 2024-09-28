
# How to use the Docker files contained in this directory

Currently, the Amanzi-ATS Docker image stack consists of three different images:
1) A "base" image containing compilers, CMake, and compiled versions of the third party libraries (metsi/amanzi-tpls)
2) An image containing a compiled verison of amanzi (metsi/amanzi)
3) An image containing a compiled version of ats (metsi/ats)

Images 2 and 3 are created with each push and pull request in the amanzi/amanzi and amanzi/ats repos through
GitHub actions CI, and are pushed to Docker Hub. Currently, only the "latest" tag is kept from each branch,
built from the most recent commits. Image 1, by contrast, only needs to be built occasionally and usually
when there is an update to a version in the TPLs.

Currently, these images are only built on Docker Hub for amd64/x86_64 architectures and the images on Docker Hub
will be very slow for users with Apple Computers with M-series chips (arm64 architecture) since the code will have
to run through emulation. There may be options eventually to have multiarchitecture builds through CI, but in the
interim, users with M-series chips may want to build their own versions of these images if one for arm64 is not
available on Docker Hub.

This directory contains some bash scripts that can be used to build the container stack locally on your system:

1. How to build and update the TPL Docker image

```
cd amanzi/Docker
docker system prune -a
```
   - This will delete all stopped containers and unused images currently on your system.
   - Use this step to remove old build artifacts and have a cleaner build environment.

- update paths as necessary in the ```./deploy_tpls_docker.sh``` script to match your system
- run ```./deploy_tpls_docker.sh```

If you are not pushing your image to the metsi/amanzi-tpls DockerHub repository, this is as far as you need to go for the TPLs image.

```
docker tag metsi/amanzi-tpls:latest metsi/amanzi-tpls:<insert version number here>
docker push metsi/amanzi-tpls
```
   - For this step, you will need access to the DockerHub repository. Login information can be obtained from David Moulton.
   - This last step pushes the newly-rebuilt-from-scratch docker image (with both tagnames) to the repository located here: [https://hub.docker.com/r/metsi/amanzi-tpls]
   - THIS WILL OVERWRITE ANY IMAGES ALREADY IN THE DOCKERHUB REPOSITORY THAT HAVE THE SAME TAGNAMES
   - If you go to the repository URL listed above and click on “tags” you will see all the different images available, one for each TPL version. Here, you can check that your push went through successfully.

2. How to build and update the Amanzi Docker image

```
cd amanzi/Docker
```
Currently, Dockerfile-Amanzi pulls from the metsi/amanzi-tpls repository (URL listed above) whatever image is labeled as “latest”, and this is what is used in the Github Actions build/tests.

If you want to use an older version of the TPLs instead, all you need to do is change the tag in the first line of Dockerfile-Amanzi from “latest” to the version number you would like to use.
   - make sure to reset to "latest" before pulling changes into master

If you already have a Docker image on your local system tagged metsi/amanzi-tpls:latest, the build will use this image instead of downloading from the URL. Similarly for images with version numbers in the tagname.
   - Use this option for testing builds of Amanzi on a branch other than master.
   - If the TPL image was built on a non-master branch, the Amanzi build will use the same branch and pull in the latest commit.

```
docker build -f Dockerfile-Amanzi -t metsi/amanzi:latest .
```
   - If your system uses a proxy, you will have to add additional arguments to the build command like so:
       ```
       docker build --build-arg http_proxy=<proxy:port> --build-arg https_proxy=<proxy:port> -f Dockerfile-Amanzi -t metsi/amanzi:latest .
       ```
   - The -f tells Docker which file contains the build instructions you want to use
   - The -t tags the final image to something useful, otherwise the image name will default to some hash

*** If you are not pushing your image to the metsi/amanzi DockerHub repository, this is as far as you need to go. ***

```
docker tag metsi/amanzi:latest metsi/amanzi:<insert version number here>
docker push metsi/amanzi
```
   - For this step, you will need access to the DockerHub repository. Login information can be obtained from David Moulton or Alexis Perry.
   - This last step pushes the newly-rebuilt-from-scratch docker image (with both tagnames) to the repository located here: https://hub.docker.com/r/metsi/amanzi .
   - THIS WILL OVERWRITE ANY IMAGES ALREADY IN THE DOCKERHUB REPOSITORY THAT HAVE THE SAME TAGNAMES
   - If you go to the repository URL listed above and click on “tags” you will see all the different images available, one for each TPL version. Here, you can check that your push went through successfully.
   - This repository is private and will require logging in to DockerHub with the same login information as above in order to view it.


3. How to run the TPL container and build Amanzi interactively (for debugging purposes)

```
docker run -it metsi/amanzi-tpls:<tagname>
```

You will automatically be in the /home/amanzi_user/amanzi directory. From there, build according to the last step in Dockerfile-Amanzi:

Important bootstrap.sh arguments to run correctly:
   ```
   --prefix=/home/amanzi_user/local
   --amanzi-build-dir=/home/amanzi_user/amanzi_builddir/amanzi
   --tpl-config-file=/home/amanzi_user/local/tpls/share/cmake/amanzi-tpl-config.cmake
   --with-mpi=/usr
   --enable-shared
   ```

exit (or cntl-C)
   - exits the docker container

4. Building docker containers for multiple processor architectures

Docker now provides multiple ways to build containers for multiple
architectures. The most common use case for needing this is building
images for newer Apple computers (the M-series chips, arm64) in addition
to existing x86_64 that Intel/AMD chipsets (and many hpc centers) use.

Images can either be build through emulation on a single system, one
image cross-compiled on a single system, or images can be built on
two different architectures and combined in the docker manifest later on.

a) Building multiarch on a single system via emulation (easy, but can be slow):

This uses a newer tool in built-kit, buildx, to build multiarchitecture images through emulation.
It is quite straightfoward to use buildx to create multiarch images - simply replace 'build' with 'buildx build --platform=<platforms>'
in the deploy scripts. Probably go grab a coffee or a walk; compiling the TPLs
through emulation will take several hours to complete.

b) Building multiarch on a single system using cross compilers:

(under development for amanzi/ats containers)

c) Building images separately on two different systems and combining in manifest:

Another option that can be used to build multi-architecture images is simply to build
the images on multiple machines that have different architectures, and then merge them
using the docker manifest afterwards. For example,

```
# AMD64
docker build -t metsi/amanzi-tpls:mpich-<version>-amd64 .
docker push metsi/amanzi-tpls:mpich-<version>-amd64

# ARM64
docker build -t metsi/amanzi-tpls:mpich-<version>-arm64 .
docker push metsi/amanzi-tpls:mpich-<version>-arm64
```

Now that both of these images are posted to Docker Hub, they can be combined to use
the same tag using docker manifest:
```
docker manifest create \
   metsi/amanzi-tpls:mpich-<version> \
   --amend metsi/amanzi-tpls:mpich-<version>-amd64 \
   --amend metsi/amanzi-tpls:mpich-<version>-arm64 \
&& docker manifest push metsi/amanzi-tpls:mpich-<version>
```