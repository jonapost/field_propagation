project(interpolate)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../magneticfield/include)

file(GLOB_RECURSE HEADER_FILES_TMP "*.hh" "*.inl")
	set(headers ${HEADER_FILES_TMP} )

file(GLOB_RECURSE SOURCE_FILES_TMP "*.cc")
	set(sources ${SOURCE_FILES_TMP})

add_executable( interpolate interpolate.cc ${headers} ${sources})

add_dependencies(interpolate magfield)

target_link_libraries(interpolate magfield)

target_link_libraries(interpolate ${Geant4_LIBRARIES})
