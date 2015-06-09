#/usr/bin/python
import os

#generating the .dat files :
for i in {0.1,1.0,10.0,100.0,200.0}:
	filename_suffix="_" + str(i) + "_.dat"
	no_of_steps=3000/i
	print "Step_len=",i, " | no_of_steps=", no_of_steps
	os.system("./testH 0 " + str(i) + " " + str(no_of_steps) + " > " + "Dats/exact"+filename_suffix )
	os.system("./testH 4 " + str(i) + " " +str(no_of_steps) + " > " + "Dats/bs45"+filename_suffix )
	os.system("./testH 3 " + str(i) + " " +str(no_of_steps) +" > " + "Dats/dopri45"+filename_suffix )
	os.system("./testH 2 " + str(i) + " " +str(no_of_steps) +" > " + "Dats/bs23"+filename_suffix )
	os.system("./testH 1 " + str(i) + " " +str(no_of_steps) +" > " + "Dats/rkf"+filename_suffix )
	os.system("./testH 5 " + str(i) + " " +str(no_of_steps) +" > " + "Dats/rk4class"+filename_suffix )
	os.system("./testH 6 " + str(i) + " " +str(no_of_steps) +" > " + "Dats/simpleHeum"+filename_suffix )

#calling the Gnuplot with input file (please edit externally) plot.p
# os.system("gnuplot 'ErrPlot.p' ")



# step_len=1.0
# no_of_steps=1000
# echo "Step_len="$step_len "no_of_steps="$no_of_steps
# ./testH 0 $step_len $no_of_steps > bs45_100.dat
# ./testH 4 $step_len $no_of_steps > bs45_100.dat
# ./testH 3 $step_len $no_of_steps > dopri45_100.dat
# ./testH 2 $step_len $no_of_steps > bs23_100.dat
# ./testH 1 $step_len $no_of_steps > rkf_100.dat
# ./testH 5 $step_len $no_of_steps > rk4class_100.dat
# ./testH 6 $step_len $no_of_steps > sh_100.dat

# set term svg
# set output "errcmp.svg"
# plot "dopri45_5.dat" u 2:5 w l ,\
# 	 "dopri45_1.dat" u 2:5 w l ,\
# 	 "dopri45_10.dat" u 2:5 w l,\
# 	 "dopri45_25.dat" u 2:5 w l,\
# 	 "dopri45_50.dat" u 2:5 w l, \
# 	 "dopri45_100.dat" u 2:5 w l
