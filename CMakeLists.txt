project(DiscoSNP++)

cmake_minimum_required(VERSION 2.6)

################################################################################
# The version number.
################################################################################
SET (gatb-tool_VERSION_MAJOR 2)
SET (gatb-tool_VERSION_MINOR 2)
SET (gatb-tool_VERSION_PATCH 1)

IF (DEFINED MAJOR)
    SET (gatb-tool_VERSION_MAJOR ${MAJOR})
ENDIF()
IF (DEFINED MINOR)
    SET (gatb-tool_VERSION_MINOR ${MINOR})
ENDIF()
IF (DEFINED PATCH)
    SET (gatb-tool_VERSION_PATCH ${PATCH})
ENDIF()

set (gatb-tool-version ${gatb-tool_VERSION_MAJOR}.${gatb-tool_VERSION_MINOR}.${gatb-tool_VERSION_PATCH})

# However, continuous integration has priority over local compilation
IF (DEFINED JENKINS_TAG)
    SET (gatb-tool-version ${JENKINS_TAG})
ENDIF()

################################################################################
# Define cmake modules directory
################################################################################
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
#  TOOLS GENERATION
################################################################################

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/bin)
#set(PROJECT_BINARY_DIR bin)
ADD_SUBDIRECTORY(tools)

################################################################################
#  PACKAGING
################################################################################

SET (CPACK_PACKAGE_DESCRIPTION_SUMMARY  "gatb-tool ${PROJECT_NAME}")
SET (CPACK_PACKAGE_VENDOR               "Genscale team (INRIA)")
SET (CPACK_PACKAGE_VERSION_MAJOR        "${gatb-tool_VERSION_MAJOR}")
SET (CPACK_PACKAGE_VERSION_MINOR        "${gatb-tool_VERSION_MINOR}")
SET (CPACK_PACKAGE_VERSION_PATCH        "${gatb-tool_VERSION_PATCH}")
SET (CPACK_PACKAGE_VERSION              "${gatb-tool-version}")

# We chose the kind of archive we want to generate
SET (CPACK_GENERATOR            "TGZ")
SET (CPACK_SOURCE_GENERATOR     "TGZ")

# We copy other stuff
# For creating the BINARY package we include the files we want
#INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools                         DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc                            DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data_sample                    DESTINATION .)
#INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/output_analyses                DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts                        DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/compile_discoSnp++.sh          DESTINATION .)
#INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt                 DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/run_discoSnp++.sh              DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/run_VCF_creator.sh             DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/README                         DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE                        DESTINATION .)

# For the source archive, we exclude some unwanted directories.
# For creating the SOURCE package we should ignore some files
SET (CPACK_SOURCE_IGNORE_FILES  ${CPACK_SOURCE_IGNORE_FILES} ".DS_Store" "/build" "../../cmake/Delivery.cmake" "genere_distrib.sh"  ".gitignore$" "doc/discoSnp_user_guide.odt" "doc/vcf_creator_user_guide.odt" "my_local_data" "donot_release" "simple_test.sh")

# To be done at the end.
INCLUDE (CPack)
