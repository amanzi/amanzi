reset
set terminal pdf enhanced color lw 5 
set output "Simplified-F-Area-1D-pflotran.pdf"

system "perl -w -i -p -e 'if ($. == 1 || $. == 2 || $. == 3) {s/^(.*)/#$1/}' pflotran-00*.tec"
system "date"
set time

set title "Pflotran: 5 component UO2 1D problem"

set xlabel "Distance [m]"
set xrange [0:100]
set xtics 10

set multiplot layout 1,2

set yrange [2:8]
set ylabel "pH"
set key bottom right
plot "pflotran-001.tec" u 1:4 w l t "t=0.10 y", \
     "pflotran-002.tec" u 1:4 w l t "t=0.25 y", \
     "pflotran-003.tec" u 1:4 w l t "t=0.50 y", \
     "pflotran-004.tec" u 1:4 w l t "t=1.00 y", \
     "pflotran-005.tec" u 1:4 w l t "t=1.50 y"

set yrange [0:4.0e-5]
set ylabel "Total UO_2^{+2} [mol/L]"
set key bottom right
plot "pflotran-001.tec" u 1:9 w l t "t=0.10 y", \
     "pflotran-002.tec" u 1:9 w l t "t=0.25 y", \
     "pflotran-003.tec" u 1:9 w l t "t=0.50 y", \
     "pflotran-004.tec" u 1:9 w l t "t=1.00 y", \
     "pflotran-005.tec" u 1:9 w l t "t=1.50 y"

unset multiplot
