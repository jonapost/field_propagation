######################################
  Author: Jean Jacquemier
  Contact: jean.jacquemier@lapp.in2p3.fr
######################################
last modified: March 18 2011

I/ Run CPU Benchmark in standalone

  1/ Configure your system for correct Geant4 release

  2/ recompile be1 and be2 application
  2.1/ prompt%> cd be1/build
  2.2/ prompt%> cmake ..
  2.3/ prompt%> make 
  2.4/ prompt%> make install
  
  3/ Launch production with this exact syntax for the option (release)
  prompt%> source /users/jjacquem/geant4/$G4RELEASE/install-static/bin/geant4.csh
  prompt%> source production.csh geant4-09-05-ref-XX

  4/ create histogram and text file
  prompt%> ./benchmark.py geant4-09-05-ref-XX

  5/ copy results to /afs/cern.ch/sw/geant4/user/vnivanch/verification/electromagnetic/test_cpu
  prompt%> cp *.gif /afs/cern.ch/sw/geant4/user/vnivanch/verification/electromagnetic/test_cpu
  prompt%> cp cpu_result.txt /afs/cern.ch/sw/geant4/user/vnivanch/verification/electromagnetic/test_cpu

II/ Run Benchmark with CTest/CMake/CDash

 1/ checkout benchmark in SAME  directory than Geant4 (svn+ssh://jjacquem@svn.cern.ch/reps/g4tests/trunk/benchmarks/electromagnetic)
  example:
     prompt%>ls /users/jjacquem/geant4/geant4-09-05-ref-09
    benchmarks  build  cmake  CMakeLists.txt  config  Configure  environments examples  install-static  LICENSE  ReleaseNotes  source  tests


  2/ go to benchmark/electromagnetic directory and run ctest
    prompt%> cd /users/jjacquem/geant4/geant4-09-05-ref-XX/benchmarks/electromagnetic 
  2.1/ edit g4benchmarks.cmake to change G4 version to correct one

  2.2/ run Ctest
    prompt%> ctest -S g4benchmarks.cmake

  3/ Result will be send to CDash (http://aidasoft.desy.de/CDash/index.php?project=Geant4)



