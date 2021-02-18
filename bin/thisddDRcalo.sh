#!/bin/bash
#################################################################################
#
#  Environment script for DD4hep examples - initializes DD4hep (and ROOT)
#  for package: ddDRcalo
# 
#  @author F.Gaede, DESY, 2013
#  @author M.Frank, CERN, 2015
#
#################################################################################
# Default of DD4hep is the primary installation directory
if [ ! ${DD4hep_DIR} ]; then
    export DD4hep_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a_FCC_1/DD4hep/01-12-01/x86_64-centos7-gcc8-opt;
fi;
if [ /cvmfs/sft.cern.ch/lcg/releases/clhep/2.4.1.3-78165/x86_64-centos7-gcc8-opt/lib/CLHEP-2.4.1.3 ]; then
    export CLHEP_DIR=/cvmfs/sft.cern.ch/lcg/releases/clhep/2.4.1.3-78165/x86_64-centos7-gcc8-opt/lib/CLHEP-2.4.1.3;
fi;
source ${DD4hep_DIR}/bin/thisdd4hep.sh;
#
dd4hep_parse_this ${BASH_ARGV[0]} ddDRcalo;
#
#----PATH---------------------------------------------------------------------
dd4hep_add_path    PATH ${THIS}/bin;
#----PYTHONPATH---------------------------------------------------------------
dd4hep_add_path    PYTHONPATH ${THIS}/;
#----ROOT_INCLUDE_PATH--------------------------------------------------------
dd4hep_add_path    ROOT_INCLUDE_PATH ${THIS}/include;
#----LIBRARY_PATH-------------------------------------------------------------
dd4hep_add_library_path ${THIS}/lib;
# -- need to extend dynamic search path for all external libraries:
if [  ]; then
    for lp in ; do
	dd4hep_add_library_path ${lp};
    done;
fi;
