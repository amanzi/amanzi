
./amanzi --xml_file=amanzi-u-1d-calcite.xml 2>/dev/null
vis-amanzi-1d.py -i calcite -n Ca++ -s total -o amanzi-ca.pdf -t "0 10 20 30 40 50" -g major
vis-amanzi-1d.py -i calcite -n HCO3- -s total -o amanzi-hco3.pdf -t "0 10 20 30 40 50" -g major
vis-amanzi-1d.py -i calcite -n H+ -s total -o amanzi-h.pdf -t "0 10 20 30 40 50" -g major
vis-amanzi-1d.py -i calcite -n pH -s total -o amanzi-ph.pdf -t "0 10 20 30 40 50" -g major
vis-amanzi-1d.py -i calcite -n Calcite -s mineral_vf -o amanzi-calcite.pdf -t "0 10 20 30 40 50" -g major


pflotran -pflotranin 1d-calcite.in
vis-pflotran-1d.py -i pflotran.h5 -n Ca++ -s total -t 0,2,13,24,35,46 -o pflotran-ca.pdf -g major
vis-pflotran-1d.py -i pflotran.h5 -n HCO3- -s total -t 0,2,13,24,35,46 -o pflotran-hco3.pdf -g major
vis-pflotran-1d.py -i pflotran.h5 -n H+ -s total -t 0,2,13,24,35,46 -o pflotran-h.pdf -g major
vis-pflotran-1d.py -i pflotran.h5 -n pH -s total -t 0,2,13,24,35,46 -o pflotran-ph.pdf -g major
vis-pflotran-1d.py -i pflotran.h5 -n Calcite -s mineral_vf -t 0,2,13,24,35,46 -o pflotran-calcite.pdf -g major
