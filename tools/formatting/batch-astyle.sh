#!/usr/bin/env bash

if [ -z "$1" ]; then
    echo  usage: $0 file-to-process file2 file3 ...
    exit 1
fi

astyle_file=`dirname $0`/astylerc.alpha
cpplint_options_file=`dirname $0`/cpplint_options.alpha

cpplint=`dirname $0`/../py_lib/cpplint/cpplint.py
cpplint_command="${cpplint} --filter=`cat ${cpplint_options_file}`"

echo $cpplint_command

for file in $@; do
    echo Processing $file
    astyle --options=${astyle_file} $file
    `${cpplint_command} $file`
done



    