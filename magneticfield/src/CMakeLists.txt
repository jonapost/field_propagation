
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB HEADER_FILES_TMP ../include/*.hh ../include/*.inl)
	set(headers ${HEADER_FILES_TMP})

file(GLOB SOURCE_FILES_TMP "*.cc")
	set(sources ${SOURCE_FILES_TMP})


add_library( magfield STATIC ${headers} ${sources})


