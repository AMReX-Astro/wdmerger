reset
set terminal postscript eps enhanced color
set output "plot_freefall.eps"


r0 = 5.7318715E9
pi = 3.14159265359
G  = 6.67428E-8
M = 1.40*1.9891E33

tff = (pi/2.0)*sqrt( r0**3 / (2*G*M) )

set samples 400
set xrange [1:0]
set yrange [0:0.99]

set size 0.65,0.65

set key bottom right

f(x) = (2 / pi) * ( acos( sqrt( x ) ) + sqrt( ( x ) * ( 1 - x ) ) )

set xlabel " r / r_0 "

set ylabel " t / t_{ff} "

plot f(x) title "Analytical Solution" w l lc rgb "blue" lw 3, \
"wdmerger_diag.out" using ( ($19-$18) / r0 ):($1 / tff) title "CASTRO Solution" with points pt 6 lw 1.5

!epstopdf plot_freefall.eps
