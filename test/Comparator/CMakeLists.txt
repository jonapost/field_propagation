# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(pname Comparator)
project(${pname})

find_package(Geant4 REQUIRED)

# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#copy plot.py to bin dir 
configure_file(
	${PROJECT_SOURCE_DIR}/plot.py
    ${PROJECT_BINARY_DIR}/plot.py
COPYONLY)

# Add the executable, and link it to the Geant4 libraries
add_executable(${pname} main.cc ${sources} ${headers})
target_link_libraries(${pname} ${Geant4_LIBRARIES})

target_compile_features(${pname} PRIVATE cxx_range_for)
