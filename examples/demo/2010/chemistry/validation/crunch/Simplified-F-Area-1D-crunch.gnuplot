reset
set terminal pdf enhanced color lw 5 
set output "Simplified-F-Area-1D-crunch.pdf"

system "date"
set time

set title "Crunch: 5 component UO2 1D problem"

set xlabel "Distance [m]"
set xrange [0:100]
set xtics 10

set multiplot layout 1,2

set yrange [2:8]
set ylabel "pH"
set key bottom right
plot "pH1.out" u 1:2 w l t "t=0.10 y", \
     "pH2.out" u 1:2 w l t "t=0.25 y", \
     "pH3.out" u 1:2 w l t "t=0.50 y", \
     "pH4.out" u 1:2 w l t "t=1.00 y", \
     "pH5.out" u 1:2 w l t "t=1.50 y"

set yrange [0:4.0e-5]
set ylabel "Total UO_2^{+2} [mol/L]"
set key bottom right
plot "totcon1.out" u 1:6 w l t "t=0.10 y", \
     "totcon2.out" u 1:6 w l t "t=0.25 y", \
     "totcon3.out" u 1:6 w l t "t=0.50 y", \
     "totcon4.out" u 1:6 w l t "t=1.00 y", \
     "totcon5.out" u 1:6 w l t "t=1.50 y"

unset multiplot
