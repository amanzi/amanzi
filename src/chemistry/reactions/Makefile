LIBRARIES = \
	libgeochemrxns.a

SRCS = \
	./AqueousEquilibriumComplex.cpp \
	./Beaker.cpp \
	./Geochemistry.cpp

OBJS = \
	$(SRCS:%.cpp=%.o)



EXECUTABLES = simpleMain

EXE_SRCS = \
	./SimpleMain.cpp

EXE_OBJS = \
	./$(EXE_SRCS:%.cpp=%.o)

LIBS = \
	-L./ -lgeochemrxns \
	-L../base_chemistry -lgeochembase \
	-L../solvers -lgeochemsolvers \
	-L../utility -lgeochemutil

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
libgeochemrxns.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)


simpleMain : $(EXE_OBJS) $(LIBRARIES)
	$(CXX) $(CXXFLAGS) -o $@ $(EXE_OBJS) $(LIBS) $(LDFLAGS)

#
# Make everything.
#
all : $(LIBRARIES) $(EXECUTABLES)


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
