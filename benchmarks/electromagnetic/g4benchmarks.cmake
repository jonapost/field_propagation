cmake_minimum_required(VERSION 2.6)

#---Utility Macros----------------------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4macros.cmake)

#---Make sure that VERBOSE is OFF to avoid screwing up the build performance
unset(ENV{VERBOSE})

#---General Configuration---------------------------------------------------
GET_PWD(pwd)
GET_HOST(host)
GET_NCPUS(ncpu)

#---Set the source and build directory--------------------------------------
set(CTEST_SOURCE_DIRECTORY "/users/jjacquem/geant4/geant4-09-04-ref-09")
set(CTEST_BINARY_DIRECTORY "/users/jjacquem/geant4/geant4-09-04-ref-09/build")
#---------------------------------------------------------------------------

set(CTEST_SITE "${host}")
if(WIN32)
  #set(CTEST_CMAKE_GENERATOR "NMake Makefiles")
  set(CTEST_CMAKE_GENERATOR "Visual Studio 9 2008")
else()
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  set(CTEST_BUILD_COMMAND "make -s -i -j${ncpu}") 
endif()
set(CTEST_BUILD_CONFIGURATION "Release")
GET_CONFIGURATION_TAG(tag)
set(CTEST_BUILD_NAME ${tag})

#---CDash settings----------------------------------------------------------
set(CTEST_PROJECT_NAME "Geant4")
set(CTEST_NIGHTLY_START_TIME "00:00:00 CET")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "aidasoft.desy.de")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=Geant4")
set(CTEST_DROP_SITE_CDASH TRUE)

#---Configuratiuon and Build Settings---------------------------------------
set(CTEST_CONFIGURE_COMMAND "cmake -DCMAKE_BUILD_TYPE=${CTEST_BUILD_CONFIGURATION}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DGEANT4_ENABLE_TESTING=ON")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DGEANT4_BUILD_EXAMPLES=OFF")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DGEANT4_INSTALL_DATA=ON")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DGEANT4_USE_GDML=OFF") 
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} -DXERCESC_ROOT_DIR=$ENV{XERCESC_ROOT_DIR}")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND} \"-G${CTEST_CMAKE_GENERATOR}\"")
set(CTEST_CONFIGURE_COMMAND "${CTEST_CONFIGURE_COMMAND}  ${CTEST_SOURCE_DIRECTORY}")
set(CTEST_CONFIGURATION_TYPE "Release")

#---Custom CTest settings---------------------------------------------------
set(CTEST_CUSTOM_TESTS_IGNORE test19 test29 test39 test49 
                              example-ext-biasing-b02 example-adv-eRosita)
if(WIN32)
  set(CTEST_CUSTOM_TESTS_IGNORE ${CTEST_CUSTOM_TESTS_IGNORE} example-ext-geometry-olap)
endif()

set(CTEST_TEST_TIMEOUT 1200)
set(CTEST_NOTES_FILES  ${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME})


#---Runtime environment-----------------------------------------------------
if(WIN32) 
  if(NOT CTEST_CMAKE_GENERATOR MATCHES Makefiles)
    set(_cfg /${CTEST_BUILD_CONFIGURATION})
  endif()
  set(ENV{PATH} "${CTEST_BINARY_DIRECTORY}/outputs/runtime${_cfg};$ENV{PATH}")
endif()

#---CTest commands----------------------------------------------------------
ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
ctest_start("Experimental")
#ctest_update()
ctest_configure()
ctest_build()
ctest_test(PARALLEL_LEVEL ${ncpu} INCLUDE bench-electromagnetic)
ctest_submit()
