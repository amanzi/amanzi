CXX_FLAGS = -g -O2 -DENABLE_DBC
#CXX_FLAGS = -O3

CXX_FLAGS +=  -std=c++11

TPLS_LIB = ${AMANZI_TPLS_DIR}/lib
TPLS_INCLUDE = ${AMANZI_TPLS_DIR}/include
TPLS_LIBS = -lexodus -lnetcdf -lhdf5 -lhdf5_hl -lz


SRCS = dbc.cc exceptions.cc Mesh.cc readMesh2D.cc writeMesh3D.cc main.cc

OBJS=$(SRCS:%.cc=%.o)
DEPS=$(OBJS:%.o=%.d)

extrude: $(DEPS) $(OBJS)
	mpicxx -std=c++11 ${CXX_FLAGS} -L${TPLS_LIB} ${TPLS_LIBS} -o extrude ${OBJS}

clean:
	rm -f ./*.o
	rm -f extrude
	rm -f .depend

%.d: %.cc
	mpicxx -std=c++11 $(CXX_FLAGS) -I $(TPLS_INCLUDE) -MM -MF $(patsubst %.o,%.d,$@) -o $@ $<

%.o: %.cc
	mpicxx -std=c++11 $(CXX_FLAGS) -I $(TPLS_INCLUDE) -c $< -o $@

-include $(DEPS)



