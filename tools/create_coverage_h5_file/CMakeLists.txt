project(create_coverage_h5_file)

################################################################################
# TOOL
################################################################################

# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# we give the headers directories from : 
#       - from project source 
#       - from GATB-CORE source
include_directories (src  ${gatb-core-includes})

# we define the files to be compiled
file (GLOB_RECURSE  ProjectFiles  src/*)

# we define the artifact to be built: the project binary
add_executable        (${PROJECT_NAME}  ${ProjectFiles})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIRZ})


# we define which libraries to be linked with project binary
target_link_libraries (${PROJECT_NAME}  ${gatb-core-libraries})


# We copy the project binary to the 'bin' directory
INSTALL (TARGETS ${PROJECT_NAME} DESTINATION bin)