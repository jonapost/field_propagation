to add boost to Geant4 do the following steps:
1. download boost from http://www.boost.org/
2. unpack it in /path_to_boost/
3. set the envaroment variable BOOST_ROOT=/path_to_boost/boost_1_60_0
4. add the following lines to CMakeLists.txt at the top of Geant4 directory
find_package(Boost REQUIRED)
include_directories( ${Boost_INCLUDEDIR})

