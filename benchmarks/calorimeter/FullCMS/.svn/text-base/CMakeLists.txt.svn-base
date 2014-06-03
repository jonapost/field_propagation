cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(name full_cms)
project(${name})
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(${name} EXCLUDE_FROM_ALL ${name}.cc ${sources})
target_link_libraries(${name} ${Geant4_LIBRARIES})

GEANT4_ADD_TEST(bench-calo-FullCMS-build
                PRECMD  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/cms.gdml cms.gdml
                POSTCMD ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/start.rndm start.rndm
                BUILD full_cms)
GEANT4_ADD_TEST(bench-calo-FullCMS-FTFP_BERT
                COMMAND full_cms  ${CMAKE_CURRENT_SOURCE_DIR}/bench1.g4 FTFP_BERT
                DEPENDS bench-calo-FullCMS-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(bench-calo-FullCMS-QGSP_BERT
                COMMAND full_cms  ${CMAKE_CURRENT_SOURCE_DIR}/bench1.g4 QGSP_BERT
                DEPENDS bench-calo-FullCMS-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


