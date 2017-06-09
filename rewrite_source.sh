#!bin/sh
#this script will rewrite the source and recompile it
#by Qiuchen Xie

#copy revised magneticfield directory to rewrite the one in /source

cp -r -b -v magneticfield $G4BASE/geometry/ > .copy.out 2> .copy.err

#then compile 
G4BUILD=/home/josh/Documents/geant4.10.0-build/
#G4BUILD=$G4BASE/../../geant4.10.0-build

cd ${G4BUILD}
make -j4 
make install 
