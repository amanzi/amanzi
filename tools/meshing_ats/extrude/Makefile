ifeq (${ATS_BUILD_TYPE},Debug)
	CXX_FLAGS = -g -O2 -DENABLE_DBC
else
	CXX_FLAGS = -g -O3
endif

CXX_FLAGS +=  -std=c++11

TPLS_LIB = ${AMANZI_TPLS_DIR}/lib
TPLS_INCLUDE = ${AMANZI_TPLS_DIR}/include
TPLS_LIBS = -lexodus -lnetcdf -lhdf5 -lhdf5_hl -lz

all: extrude_one_layer extrude_uniform extrude_homogeneous_uniform extrude_with_pinchouts

extrude_one_layer: extrude_one_layer.o
	make -C src
	mpicxx -std=c++11 ${CXX_FLAGS} -L${TPLS_LIB} ${TPLS_LIBS} src/extrude.a extrude_one_layer.o -o extrude_one_layer

extrude_uniform: extrude_uniform.o
	make -C src
	mpicxx -std=c++11 ${CXX_FLAGS} -L${TPLS_LIB} ${TPLS_LIBS} src/extrude.a extrude_uniform.o -o extrude_uniform

extrude_homogeneous_uniform: extrude_homogeneous_uniform.o
	make -C src
	mpicxx -std=c++11 ${CXX_FLAGS} -L${TPLS_LIB} ${TPLS_LIBS} src/extrude.a extrude_homogeneous_uniform.o -o extrude_homogeneous_uniform

extrude_with_pinchouts: extrude_with_pinchouts.o
	make -C src
	mpicxx -std=c++11 ${CXX_FLAGS} -L${TPLS_LIB} ${TPLS_LIBS} src/extrude.a extrude_with_pinchouts.o -o extrude_with_pinchouts

clean:
	rm -f ./*.o
	rm -f extrude_homogeneous_uniform extrude_with_pinchouts extrude_one_layer extrude_uniform
	rm -f .depend
	rm -f ./*.d
	make -C src clean
%.o: %.cc
	mpicxx -std=c++11 $(CXX_FLAGS) -I src -I $(TPLS_INCLUDE) -c $< -o $@




