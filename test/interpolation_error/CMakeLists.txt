project(interpolation_error)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../magneticfield/include)

file(GLOB_RECURSE HEADER_FILES_TMP "*.hh" "*.inl")
	set(headers ${HEADER_FILES_TMP} )

file(GLOB_RECURSE SOURCE_FILES_TMP "*.cc")
	set(sources ${SOURCE_FILES_TMP})

add_executable( interpolation_error interpolation_error.cc ${headers} ${sources})

add_dependencies(interpolation_error magfield)

target_link_libraries(interpolation_error magfield)

target_link_libraries(interpolation_error ${Geant4_LIBRARIES})
