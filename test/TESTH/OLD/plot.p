unset logscale
set term svg
set output "comparing_45.svg"
set title "Step Size = 400mm"
plot  "exact.dat" u 1:5 w l title "eyact-y",\
	 "sh.dat" u 1:5 w l title "sh-y",\
	 