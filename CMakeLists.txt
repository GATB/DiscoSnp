project(DiscoSNP)

cmake_minimum_required(VERSION 2.6)

################################################################################
# Define cmake modules directory
################################################################################
FOREACH (path "cmake" "../cmake"  "thirdparty/gatb-core/cmake"  "../../thirdparty/gatb-core/gatb-core/cmake")
IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/")
SET (CMAKE_MODULE_PATH  "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
ENDIF()
ENDFOREACH(path)

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
#SET (GATB_CORE_EXCLUDE_TOOLS     1)
#SET (GATB_CORE_EXCLUDE_TESTS     1)
#SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
#  TOOLS GENERATION 
################################################################################
ADD_SUBDIRECTORY(tools)
