
target=testNTST
#1 to 11
k="00"
for i in $(seq 1 6)
do
k=${k}0
for n in 4 15 
#3 4 8 13
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
/field/setMinStep 0.${k}1
/field/update
/run/beamOn  1000
/NTST/getFieldStats
/field/getChordFinderStats" > loglog_i${i}.mac

#valgrind --tool=callgrind --dump-instr=yes --trace-jump=yes --callgrind-out-file=callgrind_stepper${n}_i${i} ./$target loglog_i${i}.mac > ./data_sandro/$target.stepper${n}.${i}.out 
taskset -c 0 ./$target loglog_i${i}.mac > ./data_sandro/$target.stepper${n}.${i}.out 2> ./data_sandro/$target.stepper${n}.${i}.err

echo "stepper${n}, with 0.1^${i}" >> newdata_${n}.out 
cat ./data_sandro/$target.stepper${n}.${i}.out | grep "User=" >> ./data_sandro/newdata_${n}.out 
done

done


