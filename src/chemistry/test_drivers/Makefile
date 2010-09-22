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

CXXFLAGS += $(INCLUDES)


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
# clean - Clean up the directory.
#
clean : 
	rm -f *~ *.CKP *.ln *.BAK *.bak .*.bak \
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
