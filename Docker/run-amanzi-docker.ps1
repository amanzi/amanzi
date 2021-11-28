#!/bin/pwsh

#
#  Run Amanzi in docker container
#
#  Customize as you like, this is a simple example and assumes you want to run
#  in your current / present working directory:
#
#  - v mount_point_on_host:mount_point_on_container
#  - w sets current working directory on the container
#  --rm deletes the container on exit
#
#  There's not effort here to parse arguments.  Just assumes you are providing an xml file.
#
Param(
[string]$xml_file,[string]$n=2
)

$HOST_MNT="$PWD"
$CONT_MNT="/home/amanzi_user/work"
$CONT_PWD="/home/amanzi_user/work"

Write-Output $HOST_MNT
Write-Output ${xml_file}, ${n}

docker run --rm -v  ${PWD}:${CONT_MNT} -w ${CONT_PWD} metsi/amanzi:amanzi-1.2-latest mpirun -n ${n} amanzi --xml_file=${xml_file}
#docker run --rm -v "$HOST_MNT:$CONT_MNT:delegated" \-w $CONT_PWD metsi/amanzi:amanzi-1.2-latest mpirun -n 2 amanzi --xml_file=$xml_file
