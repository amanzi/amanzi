ifeq (${ATS_BUILD_TYPE},Debug)
	CXX_FLAGS = -g -O2 -DENABLE_DBC
else
	CXX_FLAGS = -g -O3
endif

CXX_FLAGS +=  -std=c++11

TPLS_LIB = ${AMANZI_TPLS_DIR}/lib
TPLS_INCLUDE = ${AMANZI_TPLS_DIR}/include
TPLS_LIBS = -lexodus -lnetcdf -lhdf5 -lhdf5_hl -lz


extrude: main.o
	make -C src
	mpicxx -std=c++11 ${CXX_FLAGS} -L${TPLS_LIB} ${TPLS_LIBS} src/extrude.a main.o -o extrude

clean:
	rm -f ./*.o
	rm -f extrude_*
	rm -f .depend
	rm -f ./*.d

%.o: %.cc
	mpicxx -std=c++11 $(CXX_FLAGS) -I src -I $(TPLS_INCLUDE) -c $< -o $@




