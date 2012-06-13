PFLOTRAN = ~/projects/pflotran/src/pflotran/pflotran

PFLOTRAN_FLAGS = -pflotranin

TESTS = \
	carbonate-unit-activity \
	carbonate-debye-huckel-activity \
	ca-carbonate-unit-activity \
	ca-carbonate-debye-huckel-activity \
	calcite-kinetics \
	calcite-kinetics-volume-fractions \
	general-reaction \
	surface-complexation-1 \
	surface-complexation-2

#
# tests
#

all : test

tests : test

test : $(TESTS)

carbonate-unit-activity :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

carbonate-debye-huckel-activity :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

ca-carbonate-unit-activity :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

ca-carbonate-debye-huckel-activity :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

calcite-kinetics :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

calcite-kinetics-volume-fractions :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

general-reaction :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

ion-exchange-valocchi :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

ionx-valocchi :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

surface-complexation-1 :
	$(PFLOTRAN) $(PFLOTRAN_FLAGS) $@.in > $@.cout
	mv pflotran.out $@.pfout 
	@echo 

surface-complexation-2 :
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
		*.out *.pfout *.rst *.ant *.tec

#
# clobber - Really clean up the directory.
#
clobber : clean 
	rm -f .Makedepend *.o *.il $(EXECUTABLES) $(LIBRARIES)

#
# FORCE - Null rule to force things to happen.
#
FORCE :
