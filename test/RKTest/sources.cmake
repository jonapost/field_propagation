
# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

#
# Define the RKTest Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME RKTest
    HEADERS
        RKTestDriver.hh
        RKTestField.hh
        RKTestTrack.hh
    SOURCES
        RKTestDriver.cc 
        RKTestField.cc 
        RKTestTrack.cc        
    GRANULAR_DEPENDENCIES
        G4globman
    GLOBAL_DEPENDENCIES
        G4global
    LINK_LIBRARIES
)

# List any source specific properties here

