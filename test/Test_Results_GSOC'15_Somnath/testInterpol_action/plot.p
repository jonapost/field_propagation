set term svg
set output "errcmp.svg"
set xrange [0: 850]
set yrange [-400:400]
set title "Plot at UFM z:1000 Gauss step len = 300 mm"
plot "out.dat"  u 2:3 w lp title "y vs x"
#edit this file as per your requirements