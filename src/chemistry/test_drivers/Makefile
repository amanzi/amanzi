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

TESTS = \
	carbonate-unit-activity-coeff \
	carbonate-debye-huckel \
	ca-carbonate-unit-activity-coeff \
	ca-carbonate-debye-huckel \
	calcite-kinetics

GEOCHEM = ./geochemistryTest
GEOVERB = -v 2

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

test : $(TESTS)

carbonate-unit-activity-coeff  :
	@echo ----------------------------------------
	$(GEOCHEM) -t 1 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/$@.test
	@echo ----------------------------------------
	@echo

carbonate-debye-huckel  : 
	@echo ----------------------------------------
	$(GEOCHEM) -t 2 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/$@.test
	@echo ----------------------------------------
	@echo

ca-carbonate-unit-activity-coeff : 
	@echo ----------------------------------------
	$(GEOCHEM) -t 3 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/$@.test
	@echo ----------------------------------------
	@echo

ca-carbonate-debye-huckel : 
	@echo ----------------------------------------
	$(GEOCHEM) -t 4 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/$@.test
	@echo ----------------------------------------
	@echo

calcite-kinetics : 
	@echo ----- WARNING: These results have not been verified! -----
	$(GEOCHEM) -t 5 $(GEOVERB) > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/$@.test
	@echo ----------------------------------------
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
