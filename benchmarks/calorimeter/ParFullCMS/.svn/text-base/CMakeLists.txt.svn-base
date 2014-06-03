# - CMakeLists.txt for Geant4 MT Example ParN02
#
#-----------------------------------------------------------------------
# Setup the project
set(name ParFullCMS)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(${name})

#-----------------------------------------------------------------------
# Find Geant4MT package, using the 'mt' component to enforce this.
# As with other examples, we activate all available UI and Vis drivers by 
# default.
# You can set WITH_GEANT4_UIVIS to OFF via the command line or 
# ccmake/cmake-gui to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED gdml ui_all vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

#-----------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#-----------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so that they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#-----------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${name} ParmainApplication.cc ${sources} ${headers})
target_link_libraries(${name} ${Geant4_LIBRARIES})

#-----------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build ${name}. This is so that we can run the executable directly 
# because it relies on these scripts being in the current working 
# directory.
#
set(PARFULLCMS_SCRIPTS
 start.rndm
 vis.mac
 cms.gdml
 bench0.g4
 bench1.g4
 README
  )

foreach(_script ${PARFULLCMS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#-----------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${name} DESTINATION bin)

#
# add to testing only for multihreaded 
#
if(GEANT4_BUILD_MULTITHREADED)
	GEANT4_ADD_TEST(bench-calo-ParFullCMS-build
               	 BUILD ${name})

	foreach(pl QGSP_BERT FTFP_BERT)

	  GEANT4_ADD_TEST(bench-calo-ParFullCMS-${pl}
               	 COMMAND ${name}  ${CMAKE_CURRENT_SOURCE_DIR}/bench0.g4 ${pl}
               	 DEPENDS bench-calo-ParFullCMS-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
	endforeach()
endif()
