unset logscale
set terminal svg
set output "comparingX.svg"
plot "exact.dat" u 0:2 with lines title "exact-x", \
	"rkf.dat" u 0:2 with lines title "rkf45-x", \
	"bs23.dat" u 0:2 with lines title "bs23-x", \
	"bs45.dat" u 0:2 with lines title "bs45-x", \
	"dopri45.dat" u 0:2 with lines title "dopri45-x" 