#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(be2)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
#option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(be2 be2.cc ${sources} ${headers})
target_link_libraries(be2 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(BE2_SCRIPTS
  atlasbar_20um.mac  
  atlasbar_emv_20um.mac 
  atlasbar_emv.mac
  atlasbar_emx_20um.mac  
  atlasbar_emx.mac  
  atlasbar.mac
  )

foreach(_script ${BE2_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS be2 DESTINATION $ENV{G4BIN}/$ENV{G4SYSTEM} )

#----------------------------------------------------------------------------
# Define tests if benchmarks included in Geant4 project
#
if(CMAKE_PROJECT_NAME STREQUAL Geant4)
  GEANT4_ADD_TEST(bench-electromagnetic-be2-build BUILD be2 LABEL Nightly)
  GEANT4_ADD_TEST(bench-electromagnetic-be2_standard
                  COMMAND be2 atlasbar.mac  
                  LABEL Nightly
                  DEPENDS bench-electromagnetic-be2-build 
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(bench-electromagnetic-be2_emv
                  COMMAND be2  atlasbar_emv.mac
                  LABEL Nightly
                  DEPENDS bench-electromagnetic-be2-build 
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(bench-electromagnetic-be2_emx
                  COMMAND be2 atlasbar_emx.mac
                  LABEL Nightly
                  DEPENDS bench-electromagnetic-be2-build 
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(bench-electromagnetic-be2_standard20
                  COMMAND be2  atlasbar_20um.mac
                  LABEL Nightly
                  DEPENDS bench-electromagnetic-be2-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(bench-electromagnetic-be2_emv20
                  COMMAND be2 atlasbar_emv_20um.mac
                  LABEL Nightly
                  DEPENDS bench-electromagnetic-be2-build 
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(bench-electromagnetic-be2_emx20
                  COMMAND be2 atlasbar_emx_20um.mac
                  LABEL Nightly
                  DEPENDS bench-electromagnetic-be2-build 
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
