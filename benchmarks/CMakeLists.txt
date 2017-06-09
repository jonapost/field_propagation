#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : 
# Package: Geant4.src.G4geometry.
#
# Intermediate level CMakeLists.txt - just process subdirectories
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66356 2012-12-18 09:02:32Z gcosmo $
#
#------------------------------------------------------------------------------

## Changed ( by D. Sorokin ? )

#---Add all subdirectories that contains a CMakeLists.txt file -------------------------------------
file(GLOB_RECURSE files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CMakeLists.txt)
foreach( file ${files} )
  get_filename_component(path ${file} PATH)
  if(path)
    add_subdirectory(${path})
  endif()
endforeach()
