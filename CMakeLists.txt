project(DiscoSNP++)

cmake_minimum_required(VERSION 2.6)



#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
################################################################################
# Define cmake modules directory
################################################################################
FOREACH (path "cmake" "../../cmake")
    IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
        SET (CMAKE_MODULE_PATH  "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/${path}")
    ENDIF()
ENDFOREACH(path)


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
set(PROJECT_BINARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)
#set(PROJECT_BINARY_DIR bin)
ADD_SUBDIRECTORY(tools)



################################################################################
#  PACKAGING 
################################################################################

# We chose the kind of archive we want to generate
SET (CPACK_GENERATOR            "TGZ")
SET (CPACK_SOURCE_GENERATOR     "TGZ")



# We copy other stuff 
# For creating the BINARY package we include the files we want
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tools                          DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc                            DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data_sample                    DESTINATION .)
INSTALL (DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/output_analyses                DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/compile_discoSnp++.sh          DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt                 DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/run_discoSnp++.sh              DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/run_VCF_creator.sh             DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/README                         DESTINATION .)
INSTALL (FILES     ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE                        DESTINATION .)

# For the source archive, we exclude some unwanted directories.
# For creating the SOURCE package we should ignore some files
SET (CPACK_SOURCE_IGNORE_FILES  ${CPACK_SOURCE_IGNORE_FILES} ".DS_Store" "/build" "../../cmake/Delivery.cmake" "genere_distrib.sh"  ".gitignore$" "doc/discoSnp_user_guide.odt" "doc/vcf_creator_user_guide.odt" "my_local_data" "donot_release")

# To be done at the end. 
INCLUDE (Delivery)
