#---------------------------------------------------------
# Simplified ATLAS Tile calorimeter (Fe-Sci) at 90 degrees
# 30 GeV pi-, B=3T, 500 event.
#---------------------------------------------------------
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0 
#
/gun/particle pi-
#
/gun/energy 30 GeV
#
/mydet/absorberMaterial Iron
/mydet/activeMaterial Scintillator
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 0
/mydet/absorberTotalLength 1680.0
/mydet/calorimeterRadius 840.0
/mydet/activeLayerNumber 120
/mydet/readoutLayerNumber 20
/mydet/activeLayerSize 3.0
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/setField 3.0 tesla
#
/mydet/update
#
/run/beamOn 500 
#
