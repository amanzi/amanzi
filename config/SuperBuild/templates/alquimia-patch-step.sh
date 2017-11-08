# This is run in Alquimia's source directory.
cp CMakeLists.txt CMakeLists.txt.old
sed s/\$ENV/$/g CMakeLists.txt.old > CMakeLists.txt
