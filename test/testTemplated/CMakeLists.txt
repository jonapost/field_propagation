project(test_Propagate_Mag_Field)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../magneticfield/include)

file(GLOB_RECURSE HEADER_FILES_TMP "*.hh" "*.inl")
	set(headers ${HEADER_FILES_TMP} PARENT_SCOPE)

file(GLOB_RECURSE SOURCE_FILES_TMP "*.cc")
	set(sources ${SOURCE_FILES_TMP} PARENT_SCOPE)

add_executable( testPropagateMagField testPropagateMagField.cc ${headers} ${sources})

add_dependencies(testPropagateMagField magfield)

target_link_libraries(testPropagateMagField magfield)

target_link_libraries(testPropagateMagField ${Geant4_LIBRARIES})
