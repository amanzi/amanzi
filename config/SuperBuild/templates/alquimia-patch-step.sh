# This is run in Alquimia's source directory.
cp CMakeLists.txt CMakeLists.txt.old
sed s/\$ENV/$/g CMakeLists.txt.old > CMakeLists.txt

cp alquimia/CMakeLists.txt alquimia/CMakeLists.txt.old
sed 's/ALQUIMIA_HAVE_CRUNCHFLOW)/ALQUIMIA_HAVE_CRUNCHFLOW)\n  list\(APPEND ALQUIMIA_TPLS lapack\)/g' \
    alquimia/CMakeLists.txt.old > alquimia/CMakeLists.txt
