# Build notes for Field Propagation (Nystrom-dev-error branch)
#### Preliminary Note on G4FieldTrack:
I discarded the version that was in the field_propagation project directory and chose to work with the version in the Geant4 install, although I forget the precise reasons why I did this. Everything was compiled against Geant4.10.01 without the most resent patch.

## CMakeLists.txt files:
There was not a lot of helpful documentation on CMake that I found available. The following sketches of CMakeLists.txt files are what I put together to build the project. They are probably not the best way of doing this. 

### CMakeLists.txt (in top level directory)
```
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(field_propogation)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
add_subdirectory(magneticfield)
add_subdirectory(test)
```


### CMakeLists.txt (in magneticfield/ directory)
```
project(magfield)
add_subdirectory(src)
```
### CMakeLists.txt (in test/ directory)
```
add_subdirectory(testTemplated)
add_subdirectory(compareByTime)
add_subdirectory(interpolation_error)
```
### CMakeLists.txt (in magneticfield/src directory)
- The Template classes T* where causing problems in the beginning so I instructed CMake to exclude them, since my project did not depend on them.
```
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
file(GLOB HEADER_FILES_TMP ../include/*.hh ../include/*.inl)
set(headers ${HEADER_FILES_TMP})
file(GLOB Tfiles ../include/T*.hh ../include/T*.inl)
file(GLOB SOURCE_FILES_TMP *.cc)
set(sources ${SOURCE_FILES_TMP})
file(GLOB Tfiles T*.cc)
add_library( magfield STATIC ${headers} ${sources})
```

### CMakeLists.txt (in test/testTemplated directory)
```
project(test_Propagate_Mag_Field)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../magneticfield/include)
file(GLOB_RECURSE HEADER_FILES_TMP "*.hh" "*.inl")
set(headers ${HEADER_FILES_TMP} )
file(GLOB_RECURSE SOURCE_FILES_TMP "*.cc")
set(sources ${SOURCE_FILES_TMP})
add_executable( testPropagateMagField testPropagateMagField.cc ${headers} ${sources})
add_dependencies(testPropagateMagField magfield)
target_link_libraries(testPropagateMagField magfield)
target_link_libraries(testPropagateMagField ${Geant4_LIBRARIES})
```
For the other program directories replace each occurance of `testPropagateMagField` with either `compareByTime` or with `intepolation_error` as the case may be.

## Program Arguments:

The 3 programs are are
   - test/compareByTime/compareByTime
        - Only uses the stepper, not a full stack
   - test/testTemplated/testPropagateMagField
        - Full stack    
   - test/interpolation_error/interpolation_error
        - small program to calculate interpolation error
 
Calling Parameters are:

### For compareByTime 

- stepper type, possible choices:
    - (-1) Murua5459
    - (0) MuruaRKN6459 (currently doesn't work)
    - (1) FineRKNG45
    - (2) FineRKNG34
    - (3) ChawlaSharmaRKNstepper
    - (4) G4CashKarpRKF45
    - (5) G4ClassicalRK4
    - (6) G4SimpleHeum
    - (7) DormandPrince745
    - (8) BogackiShampine45 (which might or might not work with this program at the moment)
- step length  
- Number of steps to perform
- file name for binary output
- filename for meta output (not really usefull for this program)
    
    
### For testPropagateMagField

- stepper type
    - (-1) ChawlaSharmaRKNstepper
    - (0) MuruaRKN5459 (currently doesn't work)
    - (1) FineRKNG45
    - (2) FineRKNG34
    - (3) MuruaRKN6459 (as said, doesn't work)
    - (4) G4CashKarpRKF45
    - (5) G4ClassicalRK4
    - (6) G4SimpleHeum
    - (7) DormandPrince745
    - (8) BogackiShampine45 (which might or might not work with this program at the moment)
- step length (to feed to each call of G4Propagator's `ComputeStep()` method)
- number of `ComputeStep()`'s to take
- largest possible step size (do not pass anything and there will be no effect)

Output comes in the form of the following files:

- tPMF_output1 (binary file of regular step segments)
- tPMF_meta_output1, a 2 line text file containing 2 numbers:
    - 1st number is the number of regular steps taken
    - 2nd number is the number of overshoot steps taken
- tPMF_no_function_calls1 (text file of accumulated function calls per regular step segment)
- tPMF_no_function_calls_overshoot1 (same but for overshoot segments)
- tPMF_intersection_indices1 (not all that important at the moment)
- tPMF_overshoot_segments1 (same as tPMF_output1, but for the overshoot segments)




