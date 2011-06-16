CRUNCHFLOW = ~/projects/crunchflow2007/bin/CrunchFlow2007

TESTS = \
	carbonate-unit \
	carbonate-debye-huckel \
	ca-carbonate-unit \
	ca-carbonate-debye-huckel \
	calcite

#
# tests
#

all : test

tests : test

test : $(TESTS)

carbonate-unit :
	echo $@.cfin > PestControl.ant
	$(CRUNCHFLOW) > $@.cout
	@echo 

carbonate-debye-huckel :
	echo $@.cfin > PestControl.ant
	$(CRUNCHFLOW) > $@.cout
	@echo 

ca-carbonate-unit :
	echo $@.cfin > PestControl.ant
	$(CRUNCHFLOW) > $@.cout
	@echo 

ca-carbonate-debye-huckel :
	echo $@.cfin > PestControl.ant
	$(CRUNCHFLOW) > $@.cout
	@echo 

calcite :
	echo $@.cfin > PestControl.ant
	$(CRUNCHFLOW) > $@.cout
	@echo 

calcite-large-time-step :
	echo $@.cfin > PestControl.ant
	$(CRUNCHFLOW) > $@.cout
	@echo 



#
# clean - Clean up the directory.
#
clean : 
	rm -f *~ *.CKP *.ln *.BAK *.bak .*.bak \
		*.cout \
		core errs \
		,* .emacs_* \
		tags TAGS \
		make.log MakeOut \
		*.out *.cfout *.rst *.ant

#
# clobber - Really clean up the directory.
#
clobber : clean 
	rm -f .Makedepend *.o *.il $(EXECUTABLES) $(LIBRARIES)

#
# FORCE - Null rule to force things to happen.
#
FORCE :
