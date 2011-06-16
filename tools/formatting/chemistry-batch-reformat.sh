#!/bin/bash

# Notes
#
# - the style stuff should be handled by amanzi-fix-buffer, but I
# can't get it to pick up the correct astylerc file, so duplicating it
# here for now....  
#
# - astyle and emacs treat public, private, protected indentation
# differently. emacs is "correct", so astyle must be called before
# emacs.
#
# - amanzi-fix-buffer does not save it's changes, so (save-buffer)
# needs to be added as the last step of that function, but probably
# don't want it in the interactive mode...

if [ -z "$1" ]; then
    echo  usage: $0 file-to-indent file2 file3 ...
    exit 1
fi

for file in $@; do
    echo Loading $file
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
    emacs --batch --load ~/.emacs.d/init.el --file $file \
        --execute '(amanzi-fix-buffer)'
    perl -w -i -p -e 's@} //@}  //@g' $file
    perl -w -i -p -e 's@; //@;  //@g' $file
    perl -w -i -p -e 's@, //@,  //@g' $file
    perl -w -i -p -e 's@//([\S]{1})@// $1@g' $file
    perl -w -i -p -e 's@(#endif)[ ]{1}(//)@$1  $2@g' $file
done


