set term svg
set output "pi_by_12.svg"
set xrange [0: 850]
set yrange [-400:400]
set title "Plot at UFM z:1000 Gauss theta step = pi/12"
plot "pi_by_12.dat"  u 2:3 w lp title "y vs x"
set term svg
set output "pi_by_10.svg"
set xrange [0: 850]
set yrange [-400:400]
set title "Plot at UFM z:1000 Gauss theta step = pi/10"
plot "pi_by_10.dat"  u 2:3 w lp title "y vs x"
set term svg
set output "pi_by_8.svg"
set xrange [0: 850]
set yrange [-400:400]
set title "Plot at UFM z:1000 Gauss theta step = pi/8"
plot "pi_by_8.dat"  u 2:3 w lp title "y vs x"
set term svg
set output "pi_by_6.svg"
set xrange [0: 850]
set yrange [-400:400]
set title "Plot at UFM z:1000 Gauss theta step = pi/6"
plot "pi_by_6.dat"  u 2:3 w lp title "y vs x"
set term svg
set output "pi_by_4.svg"
set xrange [0: 850]
set yrange [-400:400]
set title "Plot at UFM z:1000 Gauss theta step = pi/4"
plot "pi_by_4.dat"  u 2:3 w lp title "y vs x"
set term svg
set output "pi_by_3.svg"
set xrange [0: 850]
set yrange [-400:400]
set title "Plot at UFM z:1000 Gauss theta step = pi/3"
plot "pi_by_3.dat"  u 2:3 w lp title "y vs x"
set term svg
set output "pi_by_2.svg"
set xrange [0: 850]
set yrange [-400:400]
set title "Plot at UFM z:1000 Gauss theta step = pi/2"
plot "pi_by_2.dat"  u 2:3 w lp title "y vs x"
#edit this file as per your requirements