#/bin/csh -v
setenv g4version $1

setenv G4RELEASE ${g4version}
# before release 09-04-ref-6 
#source /users/jjacquem/geant4/${G4RELEASE}/install-static/share/geant4-9.4.0/config/geant4-9.4.0.csh
#source /users/jjacquem/geant4/${G4RELEASE}/install-static/share/Geant4-9.5.0/geant4-environment-setup.csh
setenv G4BIN ${G4WORKDIR}/bin/${G4RELEASE}



setenv G4LEVELGAMMADATA "/users/adotti/g4datafiles/PhotonEvaporation2.1"
setenv G4RADIOACTIVEDATA "/users/adotti/g4datafiles/RadioactiveDecay3.3"
setenv G4LEDATA "/users/adotti/g4datafiles/G4EMLOW6.23"
setenv G4NEUTRONHPDATA "/users/adotti/g4datafiles/G4NDL3.14"
setenv G4ABLADATA "/users/adotti/g4datafiles/G4ABLA3.0"
setenv G4REALSURFACEDATA "/users/adotti/g4datafiles/RealSurface1.0"
setenv G4NEUTRONXSDATA "/users/adotti/g4datafiles/G4NEUTRONXS1.0"
setenv G4PIIDATA "/users/adotti/g4datafiles/G4PII1.2"



mkdir -p $g4version
cd $g4version
rm -f *
setenv g4results `pwd`
source $BENCHMARKS/electromagnetic/be1/run.csh
source $BENCHMARKS/electromagnetic/be2/run.csh
rm -f current*.*
cd ..
#

