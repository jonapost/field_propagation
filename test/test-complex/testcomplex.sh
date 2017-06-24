#!/bin/bash

#################################################################
# GDML Geomtery File 
GEOMFILE="/where/is/your/cms2015.gdml"
#################################################################

#################################################################
# ROOT file that contains the pre-generated primary events
EVENTFILE="/wher/is/your/pp14TeVminbias.root"
#################################################################

#################################################################
# Geant4 macro file with some Geant4 commands
GEANTMACRO="/where/is/your/g4macro_TPHYS.mac" 
#################################################################

#################################################################
# Low energy cut value given in GeV units
LOWENERGYCUT="0.001"
#################################################################

#################################################################
# (OPTIONAL) Physics list name: Tabulated physics is the default
PHYSLISTNAME="TABPHYS"
#################################################################

#################################################################
# Set the required score type (0 is default i.e. no scoring)
SCORETYPEFLAG="0"
#################################################################
 
./testcomplex \
--geomFile $GEOMFILE \
--eventFile $EVENTFILE \
--geantMacro $GEANTMACRO \
--lowEnergyCut $LOWENERGYCUT \
--physListName $PHYSLISTNAME \
--scoreType $SCORETYPEFLAG

