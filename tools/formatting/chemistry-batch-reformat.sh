#!/bin/bash

# Note, the style stuff should be handled by amanzi-fix-buffer, but I
# can't get it to pick up the correct astylerc file, so duplicating it
# here for now....

if [ -z "$1" ]; then
    echo  usage: $0 file-to-indent file2 file3 ...
    exit 1
fi

for file in $@; do
    echo Loading $file
    emacs --batch --load ~/.emacs.d/init.el --file $file \
        --execute '(amanzi-fix-buffer)'
    astyle --indent=spaces=2 \
           --convert-tabs \
           --brackets=attach \
           --pad-oper \
           --unpad-paren \
           --pad-header \
	   --add-brackets \
           --align-pointer=type \
           --lineend=linux \
	   $file
    perl -w -i -p -e 's@} //@}  //@g' $file
    perl -w -i -p -e 's@; //@;  //@g' $file
    perl -w -i -p -e 's@, //@,  //@g' $file
    perl -w -i -p -e 's@//([\S]{1})@// $1@g' $file
done


