reset
set terminal pdf enhanced color lw 5 
set output "Simplified-F-Area-1D-amanzi.pdf"

system "date"
set time

set title "Amanzi: 5 component UO2 1D problem"

set xlabel "Distance [m]"
set xrange [0:100]
set xtics 10

set multiplot layout 1,2

set yrange [2:8]
set ylabel "pH"
set key bottom right
plot "pH_00020.dat" w l t "t=0.11 y", \
     "pH_00048.dat" w l t "t=0.26 y", \
     "pH_00105.dat" w l t "t=0.56 y", \
     "pH_00191.dat" w l t "t=1.01 y", \
     "pH_00286.dat" w l t "t=1.51 y"

set yrange [0:4.0e-5]
set ylabel "Total UO_2^{+2} [mol/L]"
set key bottom right
plot "conc_4_UO2++_00020.dat" w l t "t=0.11 y", \
     "conc_4_UO2++_00048.dat" w l t "t=0.26 y", \
     "conc_4_UO2++_00105.dat" w l t "t=0.56 y", \
     "conc_4_UO2++_00191.dat" w l t "t=1.01 y", \
     "conc_4_UO2++_00286.dat" w l t "t=1.51 y"

unset multiplot
