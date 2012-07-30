rm -rf *.xmf surface_complexation_*.h5

./amanzi --xml_file=amanzi-u-1d-surface-complexation.xml 2>/dev/null

vis-amanzi-1d.py -i surface-complexation -n pH -s total -t "0 10 20 30 40 50" --grid major -o amanzi-pH.pdf
vis-amanzi-1d.py -i surface-complexation -n H+ -s total -t "0 10 20 30 40 50" --grid major -o amanzi-h.pdf
vis-amanzi-1d.py -i surface-complexation -n Na+ -s total -t "0 10 20 30 40 50" --grid major -o amanzi-na.pdf
vis-amanzi-1d.py -i surface-complexation -n NO3- -s total -t "0 10 20 30 40 50" --grid major -o amanzi-no3.pdf
vis-amanzi-1d.py -i surface-complexation -n Zn++ -s total -t "0 10 20 30 40 50" --grid major -o amanzi-zn.pdf
vis-amanzi-1d.py -i surface-complexation -n H+ -s sorbed -t "0 10 20 30 40 50" --grid major -o amanzi-h-sorbed.pdf
vis-amanzi-1d.py -i surface-complexation -n Na+ -s sorbed -t "0 10 20 30 40 50" --grid major -o amanzi-na-sorbed.pdf
vis-amanzi-1d.py -i surface-complexation -n NO3- -s sorbed -t "0 10 20 30 40 50" --grid major -o amanzi-no3-sorbed.pdf
vis-amanzi-1d.py -i surface-complexation -n Zn++ -s sorbed -t "0 10 20 30 40 50" --grid major -o amanzi-zn-sorbed.pdf




pflotran -pflotranin 1d-surface-complexation.in

vis-pflotran-1d.py -i pflotran.h5 -s total -n pH -t "0 2 13 24 35 46" --grid major -o pflotran-pH.pdf
vis-pflotran-1d.py -i pflotran.h5 -s total -n H+ -t "0 2 13 24 35 46" --grid major -o pflotran-h.pdf
vis-pflotran-1d.py -i pflotran.h5 -s total -n Na+ -t "0 2 13 24 35 46" --grid major -o pflotran-na.pdf
vis-pflotran-1d.py -i pflotran.h5 -s total -n NO3- -t "0 2 13 24 35 46" --grid major -o pflotran-no3.pdf
vis-pflotran-1d.py -i pflotran.h5 -s total -n Zn++ -t "0 2 13 24 35 46" --grid major -o pflotran-zn.pdf
vis-pflotran-1d.py -i pflotran.h5 -s sorbed -n H+ -t "0 2 13 24 35 46" --grid major -o pflotran-h-sorbed.pdf
vis-pflotran-1d.py -i pflotran.h5 -s sorbed -n Na+ -t "0 2 13 24 35 46" --grid major -o pflotran-na-sorbed.pdf
vis-pflotran-1d.py -i pflotran.h5 -s sorbed -n NO3- -t "0 2 13 24 35 46" --grid major -o pflotran-no3-sorbed.pdf
vis-pflotran-1d.py -i pflotran.h5 -s sorbed -n Zn++ -t "0 2 13 24 35 46" --grid major -o pflotran-zn-sorbed.pdf
