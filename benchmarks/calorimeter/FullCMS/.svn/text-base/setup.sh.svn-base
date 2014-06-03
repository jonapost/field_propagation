#--------------------------------------------------------------------
# Last update: 14-Nov-2011
#--------------------------------------------------------------------
#
export G4LIB_BUILD_GDML=1
export G4LIB_USE_GDML=1
###export XERCESCROOT=/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p1/x86_64-slc5-gcc43-opt
export XERCESCROOT=/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p1/x86_64-slc5-gcc45-opt
###export XERCESCROOT=/afs/cern.ch/sw/lcg/external/XercesC/3.1.1p1/x86_64-slc5-gcc46-opt
export LD_LIBRARY_PATH=$XERCESCROOT/lib:$LD_LIBRARY_PATH
#
#
export G4SYSTEM=Linux-g++
export G4INSTALL=/home/ribon/Work/Geant4
export G4LIB=$G4INSTALL/lib
#
###export G4DEBUG=1
###export G4FPE_DEBUG=1   
###export G4PROFILE=1
###export G4OPTIMISE=1
###export G4_NO_VERBOSE=1
###export G4LIB_BUILD_STATIC=1
###export G4LIB_BUILD_SHARED=1
####export CPPVERBOSE=1
#
###export G4DEV=/afs/cern.ch/sw/geant4/dev
#
#--- data ---
DATAG4=$G4INSTALL/Data
export G4LEVELGAMMADATA=$DATAG4/PhotonEvaporation
export G4RADIOACTIVEDATA=$DATAG4/RadioactiveDecay
export G4LEDATA=$DATAG4/G4EMLOW
export NeutronHPCrossSections=$DATAG4/G4NDL
export G4NEUTRONHPDATA=$DATAG4/G4NDL
#
#--- CLHEP ---
#export CLHEP_BASE_DIR=/home/ribon/Work/CLHEP/2.1.0.1/MyInstall
#export CLHEP_INCLUDE_DIR=$CLHEP_BASE_DIR/include
#export CLHEP_LIB_DIR=$CLHEP_BASE_DIR/lib
#export CLHEP_LIB=CLHEP
#
#--- analysis/UI/visualization ---
###export G4UI_BUILD_XM_SESSION=1
###export G4UI_USE_XM=1
###export G4UI_BUILD_XAW_SESSION=1
###export G4UI_USE_XAW=1
export G4UI_USE_TCSH=1
#
export G4VIS_NONE=1
###export G4VIS_BUILD_DAWN_DRIVER=1
###export G4VIS_BUILD_OPENGLX_DRIVER=1
###export G4VIS_BUILD_OPENGLXM_DRIVER=1
###export G4VIS_BUILD_VRML_DRIVER=1
###export G4VIS_USE_DAWN=1
###export G4VIS_USE_OPENGLX=1
###export G4VIS_USE_OPENGLXM=1
###export G4VIS_USE_VRML=1
#
# --- PATH ---
export PATH=$PATH:$G4BIN/$G4SYSTEM
###export PATH=$PI_DIR/bin:$PATH
#
# --- LD_LIBRARY_PATH ---
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$G4LIB/$G4SYSTEM
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLHEP_LIB_DIR
#
# ===================
# === APPLICATION ===
# ===================
#
# --- Setup your Geant4 application environment ---
#
export G4WORKDIR=$PWD
export G4TMP=$G4WORKDIR/tmp
export G4BIN=$G4WORKDIR/bin
#
export PATH=$PATH:$G4WORKDIR/bin/$G4SYSTEM
#
###export G4ANALYSIS_USE=1
#
