cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(testNTST)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
##include(UseGeant4_internal)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build field01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(fieldNTST_SCRIPTS
    run2a.mac run2b.mac run2c.mac run2xa.mac run2xb.mac run2xc.mac
  )
# do-runs.csh  do-runs.sh

foreach(_script ${fieldNTST_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

GEANT4_EXECUTABLE(testNTST testNTST.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(NTST-build BUILD testNTST)
GEANT4_ADD_TEST(NTST COMMAND NTST NTST.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )
