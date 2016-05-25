#/usr/bin/python
import os
os.system("mkdir Dats")
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
