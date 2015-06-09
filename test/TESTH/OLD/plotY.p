unset logscale
set term svg
set output "comparing_yPlot.svg"
plot "exact.dat" u 0:4 w l title "eyact-y", \
	"rkf.dat" u 0:4  w l title "rkf45-y", \
	"bs23.dat" u 0:4 w l title "bs23-y", \
	"bs45.dat" u 0:4 w l title "bs45-y", \
	"sh.dat" u 0:4 w l title "simpleHeum-y",\
	"rk4class.dat" u 0:4 w l title "ClassicalRk4-y",\
	"dopri45.dat" u 0:5 w l title "dopri45-y" 