objdir = ./obj

CXX              = g++
CFLAGS		 = -I ./includes
CPPFLAGS         = 

obj = ${objdir}/Species.o \
	${objdir}/AqueousEquilibriumComplex.o \
	${objdir}/Geochemistry.o \
	${objdir}/LU.o \
	${objdir}/Block.o \
	${objdir}/Glenn.o \
	${objdir}/SimpleMain.o \
	${objdir}/FileIO.o

${objdir}/AqueousEquilibriumComplex.o : \
	./reactions/AqueousEquilibriumComplex.cpp \
	./base_chemistry/Species.cpp \
	./utility/Block.cpp 
	${CXX} ${CPPFLAGS} -c ./reactions/AqueousEquilibriumComplex.cpp \
	${CFLAGS}
	mv AqueousEquilibriumComplex.o ${objdir}

${objdir}/Block.o : ./utility/Block.cpp
	${CXX} ${CPPFLAGS} -c ./utility/Block.cpp \
	${CFLAGS}
	mv Block.o ${objdir}

${objdir}/FileIO.o : ./utility/FileIO.cpp 
	${CXX} ${CPPFLAGS} -c ./utility/FileIO.cpp \
	${CFLAGS}
	mv FileIO.o ${objdir}

${objdir}/Glenn.o : ./utility/Glenn.cpp \
	./reactions/Geochemistry.cpp \
	./reactions/AqueousEquilibriumComplex.cpp \
	./base_chemistry/Species.cpp 
	${CXX} ${CPPFLAGS} -c ./utility/Glenn.cpp \
	${CFLAGS}
	mv Glenn.o ${objdir}

${objdir}/Geochemistry.o : ./reactions/Geochemistry.cpp \
	./reactions/AqueousEquilibriumComplex.cpp \
	./base_chemistry/Species.cpp \
	./solvers/LU.cpp \
	./utility/Block.cpp 
	${CXX} ${CPPFLAGS} -c ./reactions/Geochemistry.cpp \
	${CFLAGS}
	mv Geochemistry.o ${objdir}

${objdir}/LU.o : ./solvers/LU.cpp 
	${CXX} ${CPPFLAGS} -c ./solvers/LU.cpp \
	${CFLAGS}
	mv LU.o ${objdir}

${objdir}/Newton.o : ./solvers/Newton.cpp \
	./solvers/LU.cpp \
	./utility/Block.cpp 
	${CXX} ${CPPFLAGS} -c ./solvers/Newton.cpp \
	${CFLAGS}
	mv Newton.o ${objdir}

${objdir}/SimpleMain.o : ./reactions/SimpleMain.cpp \
	./reactions/Geochemistry.cpp \
	./utility/Glenn.cpp 
	${CXX} ${CPPFLAGS} -c ./reactions/SimpleMain.cpp \
	${CFLAGS}
	mv SimpleMain.o ${objdir}

${objdir}/Species.o : ./base_chemistry/Species.cpp
	${CXX} ${CPPFLAGS} -c ./base_chemistry/Species.cpp \
	${CFLAGS}
	mv Species.o ${objdir}


glenn : $(obj) 
	${CXX}   -o glenn $(obj) 


