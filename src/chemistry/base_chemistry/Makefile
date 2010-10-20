LIBRARIES = \
	libgeochembase.a

SRCS = \
	./ActivityModelFactory.cpp \
	./ActivityModel.cpp \
	./ActivityModelDebyeHuckel.cpp \
	./ActivityModelUnit.cpp \
	./IonExchangeSite.cpp \
	./Species.cpp 

OBJS = \
	$(SRCS:%.cpp=%.o)

INCLUDES = -I. \
	-I../includes 

CXXFLAGS += -g -Wall -ansi -Wextra $(INCLUDES)


#
# Generic rules to generate executables.
#
% : %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

#
# generic build rule for all compilers
#
libgeochembase.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)

#
# Make everything.
#
all : $(LIBRARIES)


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
	rm -f .Makedepend *.o *.il $(LIBRARIES)

#
# FORCE - Null rule to force things to happen.
#
FORCE :
