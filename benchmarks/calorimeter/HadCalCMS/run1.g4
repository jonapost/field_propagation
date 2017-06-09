#-------------------------------------------------------
#  Simplified ATLAS HEC calorimeter (Cu-LAr), 
#  30 GeV pi-, B=2T, 5 event.
#-------------------------------------------------------
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0 
#
/gun/particle pi-
#
/gun/energy 30 GeV
#
/mydet/absorberMaterial Copper
/mydet/activeMaterial LiquidArgon
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 9.96
/mydet/calorimeterRadius 4.98
/mydet/activeLayerNumber 60
/mydet/readoutLayerNumber 20
/mydet/activeLayerSize 8.5
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/setField 2.0 tesla
#
/mydet/update
#
/run/beamOn 5 
#
