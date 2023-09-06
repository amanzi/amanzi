#!/usr/bin/env bash

# System info can be handy, and we use it for the log file associated with the installation
OS_NAME=Ubuntu-20.0.4-amd64

# Version and archive file info
ACONDA_VER=2023.07-2
ACONDA_NAME=Anaconda3-${ACONDA_VER}
ACONDA_FILE_NAME=${ACONDA_NAME}-Linux-x86_64.sh

# Set top-level installation directory
#
INSTALL_DIR=/home/amanzi_user/anaconda

# Download the source code from anaconda.com
# https://www.anaconda.com/products/distribution#Downloads

# Set download directory to save space in the docker image
DOWNLOADS_DIR=/home/amanzi_user/work/container/downloads

# Set log file names
# 
LOGDIR=./logs
LOGFILE=${OS_NAME}-${ACONDA_NAME}

LOGI=${LOGDIR}/${LOGFILE}.install

# Run the sh installer with embedded contents
#
/bin/bash $DOWNLOADS_DIR/$ACONDA_FILE_NAME -b -p $INSTALL_DIR/$ACONDA_NAME > ${LOGI} 2>&1

# Add conda environment initialization to the shell startup
echo "" >> ~/.bashrc
echo "# Adding conda initializaton to the default shell." >> ~/.bashrc
echo ". /home/amanzi_user/anaconda/Anaconda3-2023.07-2/etc/profile.d/conda.sh" >> ~/.bashrc
echo "conda activate base" >> ~/.bashrc
echo "" >> ~/.bashrc




