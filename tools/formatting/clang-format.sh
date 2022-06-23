find . -name \*.hh -print -o -name \*.cc -print | xargs clang-format -i -style=file -verbose
