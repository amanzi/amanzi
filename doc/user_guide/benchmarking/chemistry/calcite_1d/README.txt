Run Amanzi
$AMANZI_EXE --xml_file=amanzi-u-1d-calcite.xml 2>/dev/null

Run Amanzi w/Alquimia interfacing for chemistry
$AMANZI_EXE --xml_file=amanzi-u-1d-calcite-alq.xml 2>/dev/null

PFloTran
pflotran -pflotranin pflotran/1d-calcite.in
