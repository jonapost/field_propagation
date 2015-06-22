project(compare_pair_by_time)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../magneticfield/include)

file(GLOB_RECURSE HEADER_FILES_TMP "*.hh" "*.inl")
	set(headers ${HEADER_FILES_TMP} )

file(GLOB_RECURSE SOURCE_FILES_TMP "*.cc")
	set(sources ${SOURCE_FILES_TMP})

add_executable( comparePairByTime comparePairByTime.cc ${headers} ${sources})

add_dependencies(comparePairByTime magfield)

target_link_libraries(comparePairByTime magfield)

target_link_libraries(comparePairByTime ${Geant4_LIBRARIES})
