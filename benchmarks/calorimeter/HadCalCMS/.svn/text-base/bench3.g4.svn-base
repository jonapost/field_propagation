#----------------------------------------------------------------
# Last update: 23-Oct-2007
#
# Adapted from  cmsHCAL.g4  (see it for more information) .
#
# We use here 50 GeV pi- , 4000 events, with 4 T magnetic field,
# and without saving the seed. 
#----------------------------------------------------------------
#
/random/resetEngineFrom start.rndm
#
/run/verbose 1 
/event/verbose 0 
/tracking/verbose 0
#
/gun/particle pi-
#
/gun/energy 50 GeV
#
/mydet/absorberMaterial Copper
/mydet/activeMaterial Scintillator
/mydet/isCalHomogeneous 0
#
/mydet/isUnitInLambda 1
/mydet/absorberTotalLength 9.9625
/mydet/calorimeterRadius 4.98125
/mydet/activeLayerNumber 25
/mydet/activeLayerSize 4.0
/mydet/isRadiusUnitInLambda 1
/mydet/radiusBinSize 0.1
/mydet/radiusBinNumber 10
#
/mydet/setField 4.0 tesla
#
/mydet/update
#
/run/beamOn 4000
#
