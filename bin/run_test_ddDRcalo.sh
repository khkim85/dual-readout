#!/bin/bash
#
# Simple script to run DD4hep tests
# - sources this${PackageName}.sh first and then
#   calls the command (given as first argument)
#   with all following arguments
#

#----- initialize environment for this package - including DD4hep 
export DD4hepExamplesINSTALL=/home/khkim4/kh-readout;
source /home/khkim4/kh-readout/bin/thisddDRcalo.sh;
#----- parse command line - first argument is the 
#      test to run
command=$1
theargs=""
shift
for i in "$@" ; do
    if [ $i != ${command} ] ; then 
	theargs="${theargs} $i"
    fi
done

echo " #### LD_LIBRARY_PATH = :  ${LD_LIBRARY_PATH}"

echo " ### running test :  '${command} ${theargs}'"
exec ${command} ${theargs}
