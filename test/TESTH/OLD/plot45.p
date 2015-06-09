unset logscale
set term svg
set output "comparing_45.svg"
plot  "exact.dat" u 2:5 w l title "eyact-xy",\
	 "dopri45.dat" u 2:5 w l title "dopri45-xy",\
	 "bs45.dat" u 2:5 w l title "bs45-xy"