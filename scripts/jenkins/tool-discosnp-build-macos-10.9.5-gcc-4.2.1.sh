#!/bin/bash
#--------------------------------------------------------------#
#         Continuous integration script for Jenkins            #
#--------------------------------------------------------------#
#
# Default mode :
# This script will exit with error (exit code 1) if any of its steps fails.
# To change this behaviour, choose DO_NOT_STOP_AT_ERROR in Jenkins (see below).
#--------------------------------------------------------------#
set +xv

echo "
-----------------------------------------
 Miscellaneous information 
-----------------------------------------
date      : `date`
hostname  : `hostname`
pwd       : `pwd`

-----------------------------------------
 Jenkins build parameters (user defined)
-----------------------------------------
BRANCH_TO_BUILD      : ${BRANCH_TO_BUILD}
INRIA_FORGE_LOGIN    : ${INRIA_FORGE_LOGIN}
DO_NOT_STOP_AT_ERROR : ${DO_NOT_STOP_AT_ERROR}

-----------------------------------------
 Jenkins build parameters (built in)
-----------------------------------------
BUILD_NUMBER         : ${BUILD_NUMBER}
JENKINS_HOME         : ${JENKINS_HOME}
WORKSPACE            : ${WORKSPACE}
"

error_code () { [ "$DO_NOT_STOP_AT_ERROR" = "true" ] && { return 0 ; } }

[ "$DO_NOT_STOP_AT_ERROR" != "true" ] && { set -e ; } || { echo "(!) DEBUG mode, the script will NOT stop..." ; echo; }
set -xv

# quick look at resources
#-----------------------------------------------
sw_vers -productVersion
#-----------------------------------------------
system_profiler SPSoftwareDataType
#-----------------------------------------------
lstopo
#-----------------------------------------------
top -l 1|head -15
#-----------------------------------------------


################################################################
#                       COMPILATION                            #
################################################################

gcc --version
g++ --version

[ `gcc -dumpversion` = 4.2.1 ] && { echo "GCC 4.2.1"; } || { echo "GCC version is not 4.2.1, we exit"; exit 1; }

JENKINS_TASK=tool-${TOOL_NAME}-build-macos-10.9.5-gcc-4.2.1
GIT_DIR=/builds/workspace/$JENKINS_TASK/gatb-${TOOL_NAME}
#N.B. /scratchdir not yet mounted on the osx slave (ciosx).
#     as soon as /scratchdir is created, one has to update TEST procedure, below.
#     refer to linux build target to see how to do that
BUILD_DIR=$GIT_DIR/build

rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR

#-----------------------------------------------
# we need gatb-core submodule to be initialized
cd $GIT_DIR
git submodule init
git submodule update

#-----------------------------------------------
cd $BUILD_DIR

#-----------------------------------------------
cmake -Wno-dev -DJENKINS_TAG=${BRANCH_TO_BUILD} $GIT_DIR

#-----------------------------------------------
make -j 2 || error_code

################################################################
#                       TEST                                   #
################################################################

cd ../test
./simple_test.sh || error_code
cd ../build

################################################################
#                       PACKAGING                              #
################################################################

# Prepare and upload bin and source bundle to the forge
if [ $? -eq 0 ] && [ "$INRIA_FORGE_LOGIN" != none ] && [ "$DO_NOT_STOP_AT_ERROR" != true ]; then
    make package
    make package_source

    # make both tar.gz available as Jenkins build artifacts    
    scp ${ARCHIVE_NAME}-${BRANCH_TO_BUILD}-bin-Darwin.tar.gz ${WORKSPACE}/
    scp ${ARCHIVE_NAME}-${BRANCH_TO_BUILD}-Source.tar.gz ${WORKSPACE}/
fi


