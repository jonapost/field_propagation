
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB HEADER_FILES_TMP ../include/*.hh ../include/*.inl)
	set(headers ${HEADER_FILES_TMP})

file(GLOB Tfiles ../include/T*.hh ../include/T*.inl)

# message( "Tfiles from Headers: ${Tfiles}" )

list(REMOVE_ITEM headers ${Tfiles})

message( "Headers: ${headers}" )

file(GLOB SOURCE_FILES_TMP *.cc)
	set(sources ${SOURCE_FILES_TMP})

file(GLOB Tfiles T*.cc)

# message( "Tfiles from Sources: ${Tfiles}" )

list(REMOVE_ITEM sources ${Tfiles})

# message( "Sources: ${sources}" )

add_library( magfield STATIC ${headers} ${sources})


