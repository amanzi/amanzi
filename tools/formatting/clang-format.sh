CLANG_FORMAT_MAJOR_VERSION=`clang-format --version | sed -E 's/.*version ([0-9]*).*/\1/'`

EXPECTED_CLANG_FORMAT_VERSION="19"

if [ "$CLANG_FORMAT_MAJOR_VERSION" = "$EXPECTED_CLANG_FORMAT_VERSION" ]; then
    find . -name \*.hh -print -o -name \*.cc -print | xargs clang-format -i -style=file -verbose
    for i in `find . -name \*.hh -print -o -name \*.cc -print`; do python ../tools/formatting/format_short_blocks.py $i; done
else
    echo "clang-format version must be $EXPECTED_CLANG_FORMAT_VERSION"
fi

