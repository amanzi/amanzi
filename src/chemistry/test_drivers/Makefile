SRCS_SPECIATE = \
	./SimpleCarbonate.cpp \
	./LargeCarbonate.cpp \
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
	simpleCarbonate \
	largeCarbonate

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

simpleCarbonate : 
	./geochemistryTest -t 1 -v 1 > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/$@.test

largeCarbonate : 
	./geochemistryTest -t 2 -v 1 > $@.cout
	-$(DIFF) $(DFLAGS) $@.cout results/$@.test


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
