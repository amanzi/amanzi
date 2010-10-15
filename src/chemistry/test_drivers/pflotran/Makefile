PFLOTRAN = ~/projects/pflotran/src/pflotran/pflotran

PFLOTRAN_FLAGS = -pflotranin

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
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

carbonate-debye-huckel :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

ca-carbonate-unit :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

ca-carbonate-debye-huckel :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

calcite :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
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
