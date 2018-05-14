========================================================
How to use the Docker files contained in this directory
========================================================
--------------------------------------------------------
1. How to build and update the TPL Docker image
--------------------------------------------------------

cd amanzi/Docker

docker system prune -a
   - This will delete all containers and images currently on your system.
   - Use this step to ensure a totally clean build environment

To build from the latest commit to the master branch, leave Dockerfile-TPLs unchanged.  To build from the latest commit to some other branch, make the following change to Dockerfile-TPLs:
   - Line 33 original: RUN git clone https://github.com/amanzi/amanzi.git
   - Line 33 changed to: RUN git clone https://github.com/amanzi/amanzi.git && git checkout <branch name>
   - reset to original before pulling changes into master

docker build -f Dockerfile-TPLs -t metsi/amanzi-tpls:latest .
   - If your system uses a proxy, you will have to add additional arguments to the build command like so:
       docker build --build-arg http_proxy=<proxy:port> --build-arg https_proxy=<proxy:port> -f Dockerfile-TPLs -t metsi/amanzi-tpls:latest .
   - The -f tells Docker which file contains the build instructions you want to use
   - The -t tags the final image to something useful, otherwise the image name will default to some hash

*** If you are not pushing your image to the metsi/amanzi-tpls DockerHub repository, this is as far as you need to go. ***

docker tag metsi/amanzi-tpls:latest metsi/amanzi-tpls:<insert version number here>
 
docker push metsi/amanzi-tpls
   - For this step, you will need access to the DockerHub repository. Login information can be obtained from David Moulton or Alexis Perry.
   - This last step pushes the newly-rebuilt-from-scratch docker image (with both tagnames) to the repository located here: https://hub.docker.com/r/metsi/amanzi-tpls .
   - THIS WILL OVERWRITE ANY IMAGES ALREADY IN THE DOCKERHUB REPOSITORY THAT HAVE THE SAME TAGNAMES
   - If you go to the repository URL listed above and click on “tags” you will see all the different images available, one for each TPL version. Here, you can check that your push went through successfully.


--------------------------------------------------------
2. How to build and update the Amanzi Docker image
--------------------------------------------------------

cd amanzi/Docker

Currently, Dockerfile-Amanzi pulls from the metsi/amanzi-tpls repository (URL listed above) whatever image is labeled as “latest”, and this is what is used in the Travis CI build/tests.

If you want to use an older version of the TPLs instead, all you need to do is change the tag in the first line of Dockerfile-Amanzi from “latest” to the version number you would like to use.
   - make sure to reset to "latest" before pulling changes into master

If you already have a Docker image on your local system tagged metsi/amanzi-tpls:latest, the build will use this image instead of downloading from the URL. Similarly for images with version numbers in the tagname.
   - Use this option for testing builds of Amanzi on a branch other than master.
   - If the TPL image was built on a non-master branch, the Amanzi build will use the same branch and pull in the latest commit.

docker build -f Dockerfile-Amanzi -t metsi/amanzi:latest .
   - If your system uses a proxy, you will have to add additional arguments to the build command like so:
       docker build --build-arg http_proxy=<proxy:port> --build-arg https_proxy=<proxy:port> -f Dockerfile-Amanzi -t metsi/amanzi:latest .
   - The -f tells Docker which file contains the build instructions you want to use
   - The -t tags the final image to something useful, otherwise the image name will default to some hash

*** If you are not pushing your image to the metsi/amanzi DockerHub repository, this is as far as you need to go. ***

docker tag metsi/amanzi:latest metsi/amanzi:<insert version number here>
 
docker push metsi/amanzi
   - For this step, you will need access to the DockerHub repository. Login information can be obtained from David Moulton or Alexis Perry.
   - This last step pushes the newly-rebuilt-from-scratch docker image (with both tagnames) to the repository located here: https://hub.docker.com/r/metsi/amanzi .
   - THIS WILL OVERWRITE ANY IMAGES ALREADY IN THE DOCKERHUB REPOSITORY THAT HAVE THE SAME TAGNAMES
   - If you go to the repository URL listed above and click on “tags” you will see all the different images available, one for each TPL version. Here, you can check that your push went through successfully.
   - This repository is private and will require logging in to DockerHub with the same login information as above in order to view it.


--------------------------------------------------------
3. How to run the TPL container and build Amanzi interactively (for debugging purposes)
--------------------------------------------------------
docker run -it metsi/amanzi-tpls
   - optionally, docker run -it metsi/amanzi-tpls:<tagname>

You will automatically be in the /home/amanzi_user/amanzi directory. From there, build according to the last step in Dockerfile-Amanzi:

Important bootstrap.sh arguments to run correctly:
   --prefix=/home/amanzi_user/local
   --amanzi-build-dir=/home/amanzi_user/amanzi_builddir/amanzi
   --tpl-config-file=/home/amanzi_user/local/tpls/share/cmake/amanzi-tpl-config.cmake
   --with-mpi=/usr
   --enable-shared

Important bootstrap.sh arguments to ensure the Travis CI build does not take too long:
   --disable-structured

exit
   - exits the docker container
