project(modifiedTestH)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../magneticfield/include)

file(GLOB_RECURSE HEADER_FILES_TMP "*.hh" "*.inl")
	set(headers ${HEADER_FILES_TMP} PARENT_SCOPE)

file(GLOB_RECURSE SOURCE_FILES_TMP "*.cc")
	set(sources ${SOURCE_FILES_TMP} PARENT_SCOPE)

add_executable( modifiedtestH modifiedtestH.cc ${headers} ${sources})

add_dependencies(modifiedtestH magfield)

target_link_libraries(modifiedtestH magfield)

target_link_libraries(modifiedtestH ${Geant4_LIBRARIES})
