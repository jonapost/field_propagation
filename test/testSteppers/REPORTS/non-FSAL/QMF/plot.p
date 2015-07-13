set term svg
set output "0.01tesla/cm.dat"
set xlabel "-log(epsilon) - Accuracy"
set ylabel "No. of field evaluations"
set logscale y
set key left top
plot "0.01.dat" u 1:3 w l title "G4CashKarpRKF45",\
"0.01.dat" u 1:5 w l title "DormandPrince745",\
"0.01.dat" u 1:2 w l title "G4ClassicalRK4",\
"0.01.dat" u 1:8 w l title "VernerRK67",\
"0.01.dat" u 1:7 w l title "VernerRK56",\
"0.01.dat" u 1:9 w l title "VernerRK78",\
"0.01.dat" u 1:4 w l title "BogackiShampine23",\
"0.01.dat" u 1:6 w l title "BogackiShampine45"
 #set logscale -- ignore
 #set logscale y
 #"0.01.dat" u 1:2 w l title "Exact",\
 #"0.01.dat" u 1:8 w l title "G4SimpleHeum",\
 #"0.01.dat" u 1:9 w l title "G4ExplicitEuler",\

# This is just the template of a GNUplot file for plotting the
# .dat files. Please modify it according to needs.