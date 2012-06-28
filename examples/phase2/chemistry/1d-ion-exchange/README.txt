./amanzi --xml_file=amanzi-u-1d-ion-exchange.xml 2>/dev/null
vis-amanzi-1d.py -i ion-exchange -n Na+ -s total -t "0 10 20 30 40 50" --grid major -o amanzi-na.pdf
vis-amanzi-1d.py -i ion-exchange -n Ca++ -s total -t "0 10 20 30 40 50" --grid major -o amanzi-ca.pdf
vis-amanzi-1d.py -i ion-exchange -n Mg++ -s total -t "0 10 20 30 40 50" --grid major -o amanzi-mg.pdf
vis-amanzi-1d.py -i ion-exchange -n Cl- -s total -t "0 10 20 30 40 50" --grid major -o amanzi-cl.pdf

vis-amanzi-1d.py -i ion-exchange -n Na+ -s sorbed -t "0 10 20 30 40 50" --grid major -o amanzi-na-sorbed.pdf
vis-amanzi-1d.py -i ion-exchange -n Ca++ -s sorbed -t "0 10 20 30 40 50" --grid major -o amanzi-ca-sorbed.pdf
vis-amanzi-1d.py -i ion-exchange -n Mg++ -s sorbed -t "0 10 20 30 40 50" --grid major -o amanzi-mg-sorbed.pdf
vis-amanzi-1d.py -i ion-exchange -n Cl- -s sorbed -t "0 10 20 30 40 50" --grid major -o amanzi-cl-sorbed.pdf


pflotran -pflotranin 1d-ion-exchange.in
vis-pflotran-1d.py -i pflotran.h5 -n Na+ -s total -t "0 2 13 24 35 46" --grid major -o pflotran-na.pdf
vis-pflotran-1d.py -i pflotran.h5 -n Ca++ -s total -t "0 2 13 24 35 46" --grid major -o pflotran-ca.pdf
vis-pflotran-1d.py -i pflotran.h5 -n Mg++ -s total -t "0 2 13 24 35 46" --grid major -o pflotran-mg.pdf
vis-pflotran-1d.py -i pflotran.h5 -n Cl- -s total -t "0 2 13 24 35 46" --grid major -o pflotran-cl.pdf

vis-pflotran-1d.py -i pflotran.h5 -n Na+ -s sorbed -t "0 2 13 24 35 46" --grid major -o pflotran-na-sorbed.pdf
vis-pflotran-1d.py -i pflotran.h5 -n Ca++ -s sorbed -t "0 2 13 24 35 46" --grid major -o pflotran-ca-sorbed.pdf
vis-pflotran-1d.py -i pflotran.h5 -n Mg++ -s sorbed -t "0 2 13 24 35 46" --grid major -o pflotran-mg-sorbed.pdf
vis-pflotran-1d.py -i pflotran.h5 -n Cl- -s sorbed -t "0 2 13 24 35 46" --grid major -o pflotran-cl-sorbed.pdf
