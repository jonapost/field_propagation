set term svg
set output "errcmp.svg"
plot "Dats/dopri45_0.1_.dat"  u 2:5 w l ,\
	 "Dats/dopri45_1.0_.dat" u 2:5 w l ,\
	 "Dats/dopri45_10.0_.dat" u 2:5 w l,\
	 "Dats/dopri45_100.0_.dat" u 2:5 w l
