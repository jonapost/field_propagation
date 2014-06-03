#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(be1)

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
add_executable(be1 be1.cc ${sources} ${headers})
target_link_libraries(be1 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(BE1_SCRIPTS
  cms10gev_emv.mac
  cms10gev_emx.mac  
  cms10gev.mac  
  igprof.mac
  )

foreach(_script ${BE1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS be1 DESTINATION $ENV{G4BIN}/$ENV{G4SYSTEM} )


#----------------------------------------------------------------------------
# Define tests if benchmarks included in Geant4 project
#
if(CMAKE_PROJECT_NAME STREQUAL Geant4)
  GEANT4_ADD_TEST(bench-electromagnetic-be1-build BUILD be1)
  GEANT4_ADD_TEST(bench-electromagnetic-be1_standard
                  COMMAND be1  cms10gev.mac
                  DEPENDS bench-electromagnetic-be1-build 
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(bench-electromagnetic-be1_emv
                  COMMAND be1  cms10gev_emv.mac
                  DEPENDS bench-electromagnetic-be1-build 
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(bench-electromagnetic-be1_emx
                  COMMAND be1  cms10gev_emx.mac
                  DEPENDS bench-electromagnetic-be1-build 
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()


