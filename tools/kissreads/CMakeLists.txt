project(kissreads)

# POSSIBILITIES:
# * compiled for generic usage. By default: kmer_spanning (a position is considered as covered if the kmer starting at this position is fully covered by a read)
# * cmake -DLIVIUS_TESTS=1 ..        : output only motifs where au-vb is specific to one datasets and av'-u'b is specific to the other
# * cmake -DCLASSICAL_SPANNING=1 ..  : compiled without kmer (a position is considered as covered as soon as a read maps this position)
# * cmake -DINPUT_FROM_KISSPLICE=1 .. : compiled for dealing with kissplice output and to count separately junctions and central portions

################################################################################

FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    add_definitions ("-Wall -O3 -DOMP")
else ()
    add_definitions ("-Wall -O3")
endif()

# we define the files to be compiled
file (GLOB_RECURSE  ProjectFiles  src/*)

# we set the include directories
include_directories (src  "include")



# we may have extra flags
if (LIVIUS_TESTS)
    add_definitions ("-DLIVIUS_TESTS")
endif()

if (CLASSICAL_SPANNING)
    add_definitions ("-DCLASSICAL_SPANNING")
endif()

if (INPUT_FROM_KISSPLICE)
    add_definitions ("-DINPUT_FROM_KISSPLICE")
endif()

# we define the artifact to be built: the project binary
add_executable  (${PROJECT_NAME}  ${ProjectFiles})

# we define which libraries to be linked with project binary
#target_link_libraries (${PROJECT_NAME}  "-lm -lz -fopenmp")
target_link_libraries (${PROJECT_NAME}  "-lm -lz")


