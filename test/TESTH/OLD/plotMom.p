unset logscale
set term svg
set output "comparing_45.svg"
plot "exact.dat" u 0:5 w l title "eyact-y", \
	"bs45.dat" u 0:5 w l title "bs45-y", \
	"dopri45.dat" u 0:5 w l title "dopri45-y",\
	"sh.dat" u 0:5 w l title "simpleHeum-y",\
	"rk4class.dat" u 0:5 w l title "ClassicalRk4-y"