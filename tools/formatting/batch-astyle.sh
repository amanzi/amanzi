if [ -z "$1" ]; then
    echo  usage: $0 file-to-process file2 file3 ...
    exit 1
fi

astyle_file=`dirname $0`/astylerc

for file in $@; do
    echo Processing $file
    astyle --options=${astyle_file} $file
done



    