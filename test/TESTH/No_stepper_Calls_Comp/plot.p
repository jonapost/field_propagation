set term svg
set output "copa-20.svg"
set xlabel "-log(epsilon) - Accuracy"
set ylabel "No. of Calls to the stepper"
plot "out.dat" u 1:2 w l title "Exact",\
"out.dat" u 1:3 w l title "G4CashKarpRKF45",\
"out.dat" u 1:5 w l title "DormandPrince745",\
"out.dat" u 1:7 w l title "G4ClassicalRK4",\
"out.dat" u 1:8 w l title "G4SimpleHeum",\
"out.dat" u 1:9 w l title "G4ExplicitEuler",\
"out.dat" u 1:4 w l title "BogackiShampine23",\
"out.dat" u 1:6 w l title "BogackiShampine45"
 #set logscale -- ignore