*********************
	INTRODUCTION
*********************
This folder Contains various tests written so far to collect and analyze data from the steppers of Geant4. Most of these steppers are new (not available in Geant4 out of the box) have been made as a part of the ongoing project (in GSOC 2015)

*********************
		ABOUT
********************* 
Each of the stepper contains a single .cc file named testH.cc which is a single source file for carrying out the test. Along with that the folder may contain a .p file which is a gnuplot script file for plotting the .dat file generated out of the test.

The general flow for running and analyse the result would be :

Compile testH.cc (You may need to adjust the dependencies) 
Run the executable testH using one of the built in script
Open the .svg file in a browser to view the graph