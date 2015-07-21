cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(RKTest)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(include)
# add_executable(RKTest RKTest.cc src/*.cc)
# if(GEANT4_BUILD_GRANULAR_LIBS)
#     include(Geant4MacroLibraryTargets)
#     GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
# endif()
GEANT4_EXECUTABLE(RKTestRun RKTestRun.cc src/*.cc include/*.hh)

# GEANT4_EXECUTABLE(RKTest test.cc src/*.cc include/*.hh)

# GEANT4_EXECUTABLE(RKTest RKTest.cc)

# GEANT4_ADD_TEST(testH COMMAND testH BUILD testH)
# GEANT4_ADD_TEST(tPMF COMMAND tPMF BUILD tPMF)
# GEANT4_ADD_TEST(testJ COMMAND testJ BUILD testJ)




