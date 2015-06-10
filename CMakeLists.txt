# (1)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(field_propogation)

# (2)
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# (3)
include(${Geant4_USE_FILE})

add_subdirectory(magneticfield)
add_subdirectory(test)

