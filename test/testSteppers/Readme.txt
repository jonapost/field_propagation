testH.cc can be compiled after bringing in the necessary dependencies/libraries
Alternatively, it can be embedded in the tests directory by forming a new directory
(testH <dir> , for instance ) with a CMakeLists.txt edited with the line

GEANT4_ADD_EXECUTABLE(testH testH.cc)

It will set up a make target which can easily be built with XCode/Eclipse by adding target testH in that directory

Alternatively the makefile in that directory can also be used

plot.p is a small GNUplot script file used to plot in the svg format the obtained output


************************************************************
In the shell run :
************************************************************

$> ./testH $stepper_no $step_size > $file_name

But plot.p has been configured to use the following ref
stepper_no '1' : rkf.dat (RKF45Stepper)
stepper_no '2' : exact.dat (ExactHelix)

to plot graph for example using stepper_no 2 and step_size


$> ./testH 1 0.1 > rkf.dat
$> ./testH 2 0.1 > exact.dat
$> gnuplot
gnuplot> load 'plot.p'
