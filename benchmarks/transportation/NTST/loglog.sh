target=testNTST
echo "output from stepper 4, 15, 8, 14" > data.txt
for n in 15
    #14 15
    #4 8
    #4 15 8 14 
do
echo " " >> data.txt
echo "case:$n" >> data.txt
echo " " >> data.txt
for i in $(seq 2 9)
do
k=${k}0
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
/field/update
/run/beamOn 1000 
/NTST/getFieldStats
/field/getChordFinderStats" >loglog.mac

#valgrind --tool=callgrind ./$target loglog.mac > ./data/$target.stepper${n}.${i}.out 
./$target loglog.mac > ./data/$target.stepper${n}.${i}.out 
#echo "stepper${n}, with 0.1^${i}" >> newdata_${n}.out 

cat ./data/$target.stepper${n}.${i}.out | grep -o "User=.*s R" | sed -s 's/User=//' | sed -s 's/s R//' >> data.txt 
done
k=0
done
exit

