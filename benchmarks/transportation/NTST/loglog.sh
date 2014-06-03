
target=testNTST

for i in $(seq 1 15)
do
k=${k}0
for n in 3 4 8 13
do
echo "/control/verbose 2 
/NTST/setDebug 0 
/run/verbose 1
/event/verbose 0 
/tracking/verbose 0 
/run/onlyTransport true
/gen/choose evt
/field/setStepperType $n
/run/initialize
/field/setMinEpsilon 0.00${k}1
/field/setMaxEpsilon 0.${k}1
#/field/setMinStep 0.${k}1
/field/update
/run/beamOn  1000
/NTST/getFieldStats
/field/getChordFinderStats" >loglog.mac

$G4BIN/$G4SYSTEM/$target loglog.mac > ./data/$target.stepper${n}.${i}.out 
echo "stepper${n}, with 0.1^${i}" >> newdata_${n}.out 
cat ./data/$target.stepper${n}.${i}.out | grep "User=" >> newdata_${n}.out 
done

done
exit

