SRCS_SPECIATE = \
	./ThermoDatabase.cpp \
	./geochemistry_test.cpp

OBJS_SPECIATE = \
	$(SRCS_SPECIATE:%.cpp=%.o)


EXECUTABLES = geochemistryTest

LIBS = \
	-L../reactions -lgeochemrxns \
	-L../base_chemistry -lgeochembase \
	-L../solvers -lgeochemsolvers \
	-L../utility -lgeochemutil

DEPEND_LIBS = \
	../reactions/libgeochemrxns.a \
	../base_chemistry/libgeochembase.a \
	../solvers/libgeochemsolvers.a \
	../utility/libgeochemutil.a

INCLUDES = -I. \
	-I../includes 

CXXFLAGS += -g -Wall -ansi -Wextra $(INCLUDES)

DIFF = diff
DFLAGS = --brief

GEOCHEM = ./geochemistryTest
GEOVERB = -v 2

TESTS_CRUNCH = \
	carbonate-unit-activity-coeff-crunch \
	carbonate-debye-huckel-crunch \
	ca-carbonate-unit-activity-coeff-crunch \
	ca-carbonate-debye-huckel-crunch \
	calcite-kinetics-crunch

TESTS_PFLOTRAN = \
	carbonate-unit-activity-coeff-pflotran \
	carbonate-debye-huckel-pflotran \
	ca-carbonate-unit-activity-coeff-pflotran \
	ca-carbonate-debye-huckel-pflotran \
	calcite-kinetics-pflotran

#
# Generic rules to generate executables.
#
% : %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

#
# generic build rule for all compilers
#

geochemistryTest : $(OBJS_SPECIATE) $(DEPEND_LIBS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_SPECIATE) $(LIBS) $(LDFLAGS)

#
# Make everything.
#
all : $(EXECUTABLES)


#
# tests
#

tests : test

test : $(TESTS_CRUNCH) $(TESTS_PFLOTRAN)

test-crunch : $(TESTS_CRUNCH)

test-pflotran : $(TESTS_PFLOTRAN)



#
# carbonate
#
carbonate-unit-activity-coeff-crunch  :
	@echo
	@echo ----------------------------------------
	$(GEOCHEM) -m crunch -t 1 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/crunch/$@.test
	@echo

carbonate-unit-activity-coeff-pflotran  :
	@echo
	@echo ----------------------------------------
	$(GEOCHEM) -m pflotran -t 1 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/pflotran/$@.test
	@echo

carbonate-debye-huckel-crunch  : 
	@echo
	@echo ----------------------------------------
	$(GEOCHEM) -m crunch -t 2 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/crunch/$@.test
	@echo

carbonate-debye-huckel-pflotran  : 
	@echo
	@echo ----------------------------------------
	$(GEOCHEM) -m pflotran -t 2 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/pflotran/$@.test
	@echo

#
# ca-carbonate
#
ca-carbonate-unit-activity-coeff-crunch : 
	@echo
	@echo ----------------------------------------
	$(GEOCHEM) -m crunch -t 3 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/crunch/$@.test
	@echo

ca-carbonate-debye-huckel-crunch : 
	@echo
	@echo ----------------------------------------
	$(GEOCHEM) -m crunch -t 4 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/crunch/$@.test
	@echo

ca-carbonate-unit-activity-coeff-pflotran : 
	@echo
	@echo ----------------------------------------
	$(GEOCHEM) -m pflotran -t 3 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/pflotran/$@.test
	@echo

ca-carbonate-debye-huckel-pflotran : 
	@echo
	@echo ----------------------------------------
	$(GEOCHEM) -m pflotran -t 4 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/pflotran/$@.test
	@echo

#
# calcite kinetics
#
calcite-kinetics-crunch : 
	@echo
	@echo ----- WARNING: These results have not been verified! -----
	$(GEOCHEM) -m crunch -t 5 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/crunch/$@.test
	@echo

calcite-kinetics-pflotran : 
	@echo
	@echo ----- WARNING: These results have not been verified! -----
	$(GEOCHEM) -m pflotran -t 5 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/pflotran/$@.test
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
		make.log MakeOut

#
# clobber - Really clean up the directory.
#
clobber : clean 
	rm -f .Makedepend *.o *.il $(EXECUTABLES) $(LIBRARIES)

#
# FORCE - Null rule to force things to happen.
#
FORCE :
