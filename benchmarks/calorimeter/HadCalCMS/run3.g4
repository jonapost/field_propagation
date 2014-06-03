#---------------------------------------------------
# Simplified hadronic calorimeter: homogeneous PbWO4
# 30 GeV pi-, B=4T, 5 event.
#---------------------------------------------------
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0 
#
/gun/particle pi-
#
/gun/energy 30 GeV
#
/mydet/absorberMaterial PbWO4
/mydet/activeMaterial PbWO4
/mydet/isCalHomogeneous 1
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 10.0
/mydet/calorimeterRadius 5.0
/mydet/activeLayerNumber 100
/mydet/readoutLayerNumber 20
/mydet/activeLayerSize 4.0
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/setField 4.0 tesla
#
/mydet/update
#
/run/beamOn 5 
#
