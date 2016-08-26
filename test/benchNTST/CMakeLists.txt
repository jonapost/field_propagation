cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(name benchNTST)
project(${name})

find_package(Boost REQUIRED)
include_directories( ${Boost_INCLUDE_DIR} )


find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
##include(UseGeant4_internal)

# Copy all scripts to the build directory, i.e. the directory in which we
# build field01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(fieldNTST_SCRIPTS
	bench1.g4
	bench2.g4
	bench3.g4 
	TableST5.dat
	babarevt.out
	SVT.dat
  )
# do-runs.csh  do-runs.sh

foreach(_script ${fieldNTST_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executable, and link it to the Geant4 libraries

add_executable(${name} testNTST.cc ${sources} ${headers})
target_link_libraries(${name} ${Geant4_LIBRARIES})

