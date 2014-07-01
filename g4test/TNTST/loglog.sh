#compile program
make clean -f Makefile
make -f Makefile
target=testNTST

echo "\n NOTE: TNTST CAN ONLY TEST TEMPLATE MODE CLASSES \n"
echo "\n Applied vectorization... \n" 
echo "Template Version Benchmark is running... \n"
echo "Template Stepper Options: 0 2 3 4 8; loading... \n"

#main script
#testing error stepper
echo "template version output from stepper 4=RK4 3=Heum 2=Runge 0=ExEuler " > data.txt
for n in 4 3 2 0 
do
echo "testing stepper option $n is running... \n "
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

#callgraph mode
#valgrind --tool=callgrind ./$target loglog.mac > ./data/$target.stepper${n}.${i}.out 

./$target loglog.mac > ./data/$target.stepper${n}.${i}.out 
#echo "stepper${n}, with 0.1^${i}" >> newdata_${n}.out 

cat ./data/$target.stepper${n}.${i}.out | grep -o "User=.*s R" | sed -s 's/User=//' | sed -s 's/s R//' >> data.txt 
done
k=0
done
exit

