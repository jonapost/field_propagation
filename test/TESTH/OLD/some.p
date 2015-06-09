some.p
"bs45.dat" u 0:2 w l title "bs45-y", \
	"dopri45.dat" u 0:2 w l title "dopri45-y",\
	"sh.dat" u 0:2 w l title "simpleHeum-y",\
	"rk4class.dat" u 0:2 w l title "ClassicalRk4-y"

	"exact.dat" u 0:2 w l title "eyact-x",\
	 "exact.dat" u 0:3 w l title "eyact-y",\

	 	 "bs45.dat" u 0:2 w l title "bs45-x",\
	 "bs45.dat" u 0:3 w l title "bs45-y",\

	 	 "dopri45.dat" u 0:2 w l title "dopri45-x",\
	 "dopri45.dat" u 0:3 w l title "dopri45-y",\

	 "exact.dat" u 2:5 w l title "eyact-xy",\
	 "dopri45.dat" u 2:5 w l title "dopri45-xy",\
	 "rkf.dat" u 2:5 w l title "rkf-xy",\
	 "bs23.dat" u 2:5 w l title "bs23-xy",\
	 "rk4class.dat" u 2:5 w l title "rk4class-xy",\
	 "sh.dat" u 2:5 w l title "sh-xy",\
	 "bs45.dat" u 2:5 w l title "bs45-xy",\

	 "bs23.dat" u 1:5 w l title "bs23-y",\

	 "bs23.dat" u 1:5 w l title "bs23-y",\
	 "rk4class.dat" u 1:5 w l title "rk4class-y",\
	 "sh.dat" u 1:5 w l title "sh-y",\
	 "bs45.dat" u 1:5 w l title "bs45-y",\
	 "exact.dat" u 1:2 w l title "eyact-x",\
	 "dopri45.dat" u 1:2 w l title "dopri45-x",\
	 "rkf.dat" u 1:2 w l title "rkf-x",\
	 "bs23.dat" u 1:2 w l title "bs23-x",\
	 "rk4class.dat" u 1:2 w l title "rk4class-x",\
	 "sh.dat" u 1:2 w l title "sh-x",\
	 "bs45.dat" u 1:2 w l title "bs45-x"

	 "dopri45.dat" u 1:5 w l title "dopri45-y",\
	 "rkf.dat" u 1:5 w l title "rkf-y",\