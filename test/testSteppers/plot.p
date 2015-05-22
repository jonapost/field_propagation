set term svg
set output "cmp.svg"
plot "rkf.dat" using 1:2 title 'xy-rkf' with lines linecolor rgb "red", \
	"exact.dat" using 1:2 title 'xy-exact' with lines linecolor rgb "blue"