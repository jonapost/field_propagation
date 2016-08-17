
#
# Macro file for "testNTST.cc"
# 
# (can be run in batch, without graphic)
#
/control/verbose 2
#/control/saveHistory
/NTST/setDebug 0
#
/run/verbose 1
/event/verbose 0
/tracking/verbose 0
#
# construct the default geometry
# 

/gen/choose evt
#
/run/initialize
/run/beamOn  1000

/NTST/getFieldStats
