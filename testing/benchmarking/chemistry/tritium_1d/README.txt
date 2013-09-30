
./amanzi --xml_file=amanzi-u-1d-tritium.xml 2>/dev/null
vis-amanzi-1d.py -i tritium -o amanzi-tritium.pdf -t "0 10 20 30 40 50" -g major



pflotran -pflotranin 1d-tritium.in
vis-pflotran-1d.py -i pflotran.h5 -n Tritium -s total -t 0,2,13,24,35,46 -o pflotran-tritium.pdf -g major
