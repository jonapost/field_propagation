cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(name mainHadCalo)
project(${name})
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(${name} EXCLUDE_FROM_ALL ${name}.cc ${sources})
target_link_libraries(${name} ${Geant4_LIBRARIES})

GEANT4_ADD_TEST(bench-calo-HadCalCMS-build
                BUILD ${name})

foreach(pl QGS_BIC QGSP_INCLXX QGSP_FTFP_BERT QGSP_BIC QGSP_BIC_HP QGSP_BIC_EMY
           QGSP_BERT_HP QBBC FTFP_BERT_TRV FTFP_BERT_HP FTF_BIC Shielding)

  GEANT4_ADD_TEST(bench-calo-HadCalCMS-${pl}
                COMMAND ${name}  ${CMAKE_CURRENT_SOURCE_DIR}/run1.g4 ${pl}
                DEPENDS bench-calo-HadCalCMS-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

foreach(pl QGSP_BERT FTFP_BERT)

  GEANT4_ADD_TEST(bench-calo-HadCalCMS-${pl}
                COMMAND ${name}  ${CMAKE_CURRENT_SOURCE_DIR}/run1_100.g4 ${pl}
                DEPENDS bench-calo-HadCalCMS-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

