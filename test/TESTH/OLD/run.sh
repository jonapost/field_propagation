#run.sh
step_len=1.0
no_of_steps=1000
echo "Step_len="$step_len "no_of_steps="$no_of_steps

./testH 0 $step_len $no_of_steps > exact_100.dat
./testH 4 $step_len $no_of_steps > bs45_100.dat
./testH 3 $step_len $no_of_steps > dopri45_100.dat
./testH 2 $step_len $no_of_steps > bs23_100.dat
./testH 1 $step_len $no_of_steps > rkf_100.dat
./testH 5 $step_len $no_of_steps > rk4class_100.dat
./testH 6 $step_len $no_of_steps > sh_100.dat