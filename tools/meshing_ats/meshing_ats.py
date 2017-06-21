"""
Extrudes a 2D mesh to generate an ExodusII 3D mesh.

Works with and assumes all polyhedra cells (and polygon faces).

To see usage, run:
------------------------------------------------------------

python meshing_ats.py -h



Example distributed with this source, to run:
------------------------------------------------------------

$> cd four-polygon-test
$> python ../meshing_ats.py -n 10 -d 1 ./four_polygon.vtk
$> mkdir run0
$> cd run0
$> ats --xml_file=../test1-fv-four-polygon.xml 



Requires building the latest version of Exodus
------------------------------------------------------------
git clone https://github.com/gsjaardema/seacas.git 

Note this requires, in turn, HDF5 and netCDF.  It is highly
recommended that you get these through anaconda or a comparable python
distribution:

conda install h5py
conda install netCDF4

```
SEACAS_BUILD_DIR=/Users/uec/codes/seacas/build-dev
SEACAS_INSTALL_DIR=/Users/uec/codes/seacas/install-dev
SEACAS_SRC_DIR=/Users/uec/codes/seacas/repos/dev

ANACONDA_DIR=/Users/uec/codes/anaconda

CC=`which mpicc`
CXX=`which mpicxx`
FC=`which mpif90`

mkdir -p $SEACAS_BUILD_DIR
mkdir -p $SEACAS_INSTALL_DIR
cd $SEACAS_BUILD_DIR

cmake  \
    -D SEACASProj_ENABLE_SEACASExodus:BOOL=ON \
    -D CMAKE_INSTALL_PREFIX:PATH=${SEACAS_INSTALL_DIR} \
    -D CMAKE_BUILD_TYPE=Debug \
    -D BUILD_SHARED_LIBS:BOOL=ON \
    \
    -D CMAKE_CXX_COMPILER:FILEPATH=${CXX} \
    -D CMAKE_C_COMPILER:FILEPATH=${CC} \
    -D CMAKE_Fortran_COMPILER:FILEPATH=${FC} \
    -D SEACASProj_SKIP_FORTRANCINTERFACE_VERIFY_TEST:BOOL=ON \
    -D TPL_ENABLE_Netcdf:BOOL=ON \
    -D TPL_ENABLE_Matio:BOOL=OFF \
    -D TPL_ENABLE_MPI=ON \
    -D TPL_ENABLE_CGNS:BOOL=OFF \
    \
    -D Netcdf_LIBRARY_DIRS:PATH=${ANACONDA_DIR}/lib \
    -D Netcdf_INCLUDE_DIRS:PATH=${ANACONDA_DIR}/include \
    -D HDF5_ROOT:PATH=${ANACONDA_DIR} \
    -D HDF5_NO_SYSTEM_PATHS=ON \
${SEACAS_SRC_DIR}


```

Excecute the configure script from your SEACAS_BUILD_DIR, then
```make``` and ```make install```
"""

import sys,os
import numpy as np
import collections
sys.path.append(os.path.join(os.environ["SEACAS_DIR"],"lib"))
import exodus
import argparse

class SideSet(object):
    def __init__(self, name, setid, elem_list, side_list):
        assert(type(setid) == int)
        assert(type(elem_list) == list or type(elem_list) == np.ndarray)
        assert(type(side_list) == list or type(side_list) == np.ndarray)

        self.name = name
        self.setid = setid
        self.elem_list = elem_list
        self.side_list = side_list

class LabeledSet(object):
    def __init__(self, name, setid, entity, ent_ids):
        assert entity in ['CELL', 'FACE', 'NODE']
        assert(type(setid) == int)
        assert(type(ent_ids) == list or type(ent_ids) == np.ndarray)

        self.name = name
        self.setid = setid
        self.entity = entity
        self.ent_ids = np.array(ent_ids)


class Mesh2D(object):
    def __init__(self, coords, connectivity, labeled_sets=None):
        """
        Creates a 2D mesh from coordinates and a list cell-to-node connectivity lists.

        coords          : numpy array of shape (NCOORDS, NDIMS)
        connectivity    : list of lists of integer indices into coords specifying a
                          (clockwise OR counterclockwise) ordering of the nodes around
                          the 2D cell
        labeled_sets    : list of LabeledSet objects
        """
        assert type(coords) == np.ndarray
        assert len(coords.shape) == 2

        self.dim = coords.shape[1]

        self.coords = coords
        self.conn = connectivity
        if labeled_sets is not None:
            self.labeled_sets = labeled_sets
        else:
            self.labeled_sets = []

        self.validate()


    def validate(self):
        assert self.coords.shape[1] == 2 or self.coords.shape[1] == 3
        assert type(self.conn) is list
        for f in self.conn:
            assert type(f) is list
            assert len(set(f)) == len(f)
            for i in f:
                assert i < self.coords.shape[0]



        for ls in self.labeled_sets:
            if ls.entity == "NODE":
                size = len(self.coords)
            elif ls.entity == "CELL":
                size = len(self.conn)

            for i in ls.ent_ids:
                assert i < size
        return True

    def num_cells(self):
        return len(self.conn)

    def num_nodes(self):
        return self.coords.shape[0]

    def num_edges(self):
        return len(self.edges())

    @staticmethod
    def edge_hash(i,j):
        return tuple(sorted((i,j)))
    
    def edges(self):
        return self.edge_counts().keys()

    def edge_counts(self):
        try:
            return self._edges
        except AttributeError:
            self._edges = collections.Counter(self.edge_hash(f[i], f[(i+1)%len(f)]) for f in self.conn for i in range(len(f)))
        return self._edges

    def plot(self, color=None, ax=None):
        if color is None:
            import colors
            cm = colors.cm_mapper(0,self.num_cells()-1)
            colors = [cm(i) for i in range(self.num_cells())]
        else:
            colors = color

        verts = [[self.coords[i,0:2] for i in f] for f in self.conn]
        from matplotlib import collections
        gons = collections.PolyCollection(verts, facecolors=colors)
        from matplotlib import pyplot as plt
        if ax is None:
            fig,ax = plt.subplots(1,1)
        ax.add_collection(gons)
        ax.autoscale_view()


    @classmethod
    def read_VTK(cls, filename):
        with open(filename,'r') as fid:
            line = fid.readline()
            while not line.startswith("POINTS"):
                line = fid.readline()
            ncoords = int(line.strip().split()[1])
            coords = np.zeros([ncoords,3],'d')
            i = 0
            while i < ncoords:
                coord_dat = np.array([float(p) for p in fid.readline().strip().split()])
                assert len(coord_dat) % 3 == 0
                ncoords_this_line = len(coord_dat) / 3
                for j in range(ncoords_this_line):
                    coords[i,:] = coord_dat[3*j:3*(j+1)]
                    i += 1

            line = fid.readline()
            while not line.startswith("POLYGONS"):
                line = fid.readline()
            ngons = int(line.strip().split()[1])
            gons = []
            for i in range(ngons):
                line = [int(n) for n in fid.readline().strip().split()[1:]]
                gon = [int(n) for n in line]

                # check handedness -- need normals to point up!
                cross = []
                for i in range(len(gon)):
                    if i == len(gon)-1:
                        ip = 0
                        ipp = 1
                    elif i == len(gon)-2:
                        ip = i+1
                        ipp = 0
                    else:
                        ip = i+1
                        ipp = i+2
                    d2 = coords[gon[ipp]] - coords[gon[ip]]
                    d1 = coords[gon[i]] - coords[gon[ip]]
                    cross.append(np.cross(d2, d1))
                if (np.array([c[2] for c in cross]).mean() < 0):
                    gon.reverse()
                    
                gons.append(gon)

        return cls(coords, gons)
            
    @classmethod
    def from_Transect(cls, x, z):
        """Creates a 2D surface strip mesh from transect data"""
        # coordinates
        y = np.array([0,1])
        Xc, Yc = np.meshgrid(x, y)
        Xc = Xc.flatten()
        Yc = Yc.flatten()

        Zc = np.concatenate([z,z])

        # connectivity
        nsurf_cells = len(x)-1
        conn = []
        for i in range(nsurf_cells):
            conn.append([i, i+1, nsurf_cells + i + 2, nsurf_cells + i + 1])

        coords = np.array([Xc, Yc, Zc])
        return cls(coords.transpose(), conn)
    

class Mesh3D(object):
    def __init__(self, coords, face_to_node_conn, elem_to_face_conn,
                 side_sets=None, labeled_sets=None, material_ids=None):
        """
        Creates a 3D mesh from coordinates and connectivity lists.

        coords            : numpy array of shape (NCOORDS, 3)
        face_to_node_conn : list of lists of integer indices into coords specifying an
                            (clockwise OR counterclockwise) ordering of the nodes around
                            the face
        elem_to_face_conn : list of lists of integer indices into face_to_node_conn
                            specifying a list of faces that make up the elem
        """
        assert type(coords) == np.ndarray
        assert len(coords.shape) == 2
        assert coords.shape[1] == 3
            
        self.dim = coords.shape[1]

        self.coords = coords
        self.face_to_node_conn = face_to_node_conn
        self.elem_to_face_conn = elem_to_face_conn

        if labeled_sets is not None:
            self.labeled_sets = labeled_sets
        else:
            self.labeled_sets = []

        if side_sets is not None:
            self.side_sets = side_sets
        else:
            self.side_sets = []
            
        if material_ids is not None:
            self.material_id_list = collections.Counter(material_ids).keys()
            self.material_ids = material_ids
        else:
            self.material_id_list = [10000,]
            self.material_ids = [10000,]*len(self.elem_to_face_conn)

        self.validate()

        
    def validate(self):
        assert self.coords.shape[1] == 3
        assert type(self.face_to_node_conn) is list
        for f in self.face_to_node_conn:
            assert type(f) is list
            assert len(set(f)) == len(f)
            for i in f:
                assert i < self.coords.shape[0]

        assert type(self.elem_to_face_conn) is list
        for e in self.elem_to_face_conn:
            assert type(e) is list
            assert len(set(e)) == len(e)
            for i in e:
                assert i < len(self.face_to_node_conn)

        for ls in self.labeled_sets:
            if ls.entity == "NODE":
                size = self.num_nodes()
            if ls.entity == "FACE":
                size = self.num_faces()
            elif ls.entity == "CELL":
                size = self.num_cells()

            for i in ls.ent_ids:
                assert i < size

        for ss in self.side_sets:
            for j,i in zip(ss.elem_list, ss.side_list):
                assert j < self.num_cells()
                assert i < len(self.elem_to_face_conn[j])



    def num_cells(self):
        return len(self.elem_to_face_conn)

    def num_faces(self):
        return len(self.face_to_node_conn)

    def num_nodes(self):
        return self.coords.shape[0]


    def write_exodus(self, filename, face_block_mode="one block"):
        """Write the 3D mesh to ExodusII using arbitrary polyhedra spec"""

        # put cells in with blocks, which renumbers the cells, so we have to track sidesets.
        # Therefore we keep a map of old cell to new cell ordering
        #
        # also, though not required by the spec, paraview and visit
        # seem to crash if num_face_blocks != num_elem_blocks.  So
        # make face blocks here too, which requires renumbering the faces.

        # -- first pass, form all elem blocks and make the map from old-to-new
        new_to_old_elems = []
        elem_blks = []
        for i_m,m_id in enumerate(self.material_id_list):
            # split out elems of this material, save new_to_old map
            elems_tuple = [(i,c) for (i,c) in enumerate(self.elem_to_face_conn) if self.material_ids[i] == m_id]
            new_to_old_elems.extend([i for (i,c) in elems_tuple])
            elems = [c for (i,c) in elems_tuple]
            elem_blks.append(elems)

        old_to_new_elems = sorted([(old,i) for (i,old) in enumerate(new_to_old_elems)], lambda a,b: int.__cmp__(a[0],b[0]))

        # -- deal with faces, form all face blocks and make the map from old-to-new
        face_blks = []
        if face_block_mode == "one block":
            # no reordering of faces needed
            face_blks.append(self.face_to_node_conn)
            
        elif face_block_mode == "n blocks, not duplicated":
            used_faces = np.zeros((len(self.face_to_node_conn),),'bool')
            new_to_old_faces = []
            for i_m,m_id in enumerate(self.material_id_list):
                # split out faces of this material, save new_to_old map
                def used(f):
                    result = used_faces[f]
                    used_faces[f] = True
                    return result

                elem_blk = elem_blks[i_m]
                faces_tuple = [(f,self.face_to_node_conn[f]) for c in elem_blk for (j,f) in enumerate(c) if not used(f)]
                new_to_old_faces.extend([j for (j,f) in faces_tuple])
                faces = [f for (j,f) in faces_tuple]
                face_blks.append(faces)

            # get the renumbering in the elems
            old_to_new_faces = sorted([(old,j) for (j,old) in enumerate(new_to_old_faces)], lambda a,b: int.__cmp__(a[0],b[0]))
            elem_blks = [[[old_to_new_faces[f][1] for f in c] for c in elem_blk] for elem_blk in elem_blks]

        elif face_block_mode == "n blocks, duplicated":
            elem_blks_new = []
            offset = 0
            for i_m, m_id in enumerate(self.material_id_list):
                used_faces = np.zeros((len(self.face_to_node_conn),),'bool')
                def used(f):
                    result = used_faces[f]
                    used_faces[f] = True
                    return result

                elem_blk = elem_blks[i_m]

                tuple_old_f = [(f,self.face_to_node_conn[f]) for c in elem_blk for f in c if not used(f)]
                tuple_new_old_f = [(new,old,f) for (new,(old,f)) in enumerate(tuple_old_f)]

                old_to_new_blk = np.zeros((len(self.face_to_node_conn),),'i')-1
                for new,old,f in tuple_new_old_f:
                    old_to_new_blk[old] = new + offset

                elem_blk_new = [[old_to_new_blk[f] for f in c] for c in elem_blk]
                #offset = offset + len(ftuple_new)

                elem_blks_new.append(elem_blk_new)
                face_blks.append([f for i,j,f in tuple_new_old_f])
            elem_blks = elem_blks_new
        elif face_block_mode == "one block, repeated":
            # no reordering of faces needed, just repeat
            for eblock in elem_blks:
                face_blks.append(self.face_to_node_conn)
        else:
            raise RuntimeError("Invalid face_block_mode: '%s', valid='one block', 'n blocks, duplicated', 'n blocks, not duplicated'"%face_block_mode)
                

        # open the mesh file
        num_elems = sum(len(elem_blk) for elem_blk in elem_blks)
        num_faces = sum(len(face_blk) for face_blk in face_blks)
        ep = exodus.ex_init_params(title=filename,
                                   num_dim=3,
                                   num_nodes=self.num_nodes(),
                                   num_face=num_faces,
                                   num_face_blk=len(face_blks),
                                   num_elem=num_elems,
                                   num_elem_blk=len(elem_blks),
                                   num_side_sets=len(self.side_sets))
        e = exodus.exodus(filename, mode='w', array_type='numpy', init_params=ep)

        # put the coordinates
        e.put_coord_names(['coordX', 'coordY', 'coordZ'])
        e.put_coords(self.coords[:,0], self.coords[:,1], self.coords[:,2])

        # put the face blocks
        for i_blk, face_blk in enumerate(face_blks):
            face_raveled = [n for f in face_blk for n in f]
            e.put_polyhedra_face_blk(i_blk+1, len(face_blk), len(face_raveled), 0)
            e.put_node_count_per_face(i_blk+1, np.array([len(f) for f in face_blk]))
            e.put_face_node_conn(i_blk+1, np.array(face_raveled)+1)

        # put the elem blocks
        assert len(elem_blks) == len(self.material_id_list)
        for i_blk, (m_id, elem_blk) in enumerate(zip(self.material_id_list, elem_blks)):
            elems_raveled = [f for c in elem_blk for f in c]

            e.put_polyhedra_elem_blk(m_id, len(elem_blk), len(elems_raveled), 0)
            e.put_elem_blk_name(m_id, "MATERIAL_ID_%d"%m_id)
            e.put_face_count_per_polyhedra(m_id, np.array([len(c) for c in elem_blk]))
            e.put_elem_face_conn(m_id, np.array(elems_raveled)+1)

        # add sidesets
        e.put_side_set_names([ss.name for ss in self.side_sets])
        for ss in self.side_sets:
            for elem in ss.elem_list:
                assert old_to_new_elems[elem][0] == elem
            new_elem_list = [old_to_new_elems[elem][1] for elem in ss.elem_list]                
            e.put_side_set_params(ss.setid, len(ss.elem_list), 0)
            e.put_side_set(ss.setid, np.array(new_elem_list)+1, np.array(ss.side_list)+1)

        # finish and close
        e.close()
        

    @classmethod
    def extruded_Mesh2D(cls, mesh2D, layer_types, layer_data, ncells_per_layer, mat_ids):
        """
        Regularly extrude a 2D mesh to make a 3D mesh.

        mesh2D              : a Mesh2D object
        layer_types         : either a string (type) or list of strings (types)
        layer_data          : array of data needed (specific to the type)
        ncells_per_layer    : either a single integer (same number of cells in all
                            : layers) or a list of number of cells in the layer
        mat_ids             : either a single integer (one mat_id for all layers)
                            : or a list of integers (mat_id for each layer)
                            : or a 2D numpy array of integers (mat_id for each layer and
                              each column: [layer_id, surface_cell_id])

        types:
          - 'constant'      : (data=float thickness) uniform thickness
          - 'function'      : (data=function or functor) thickness as a function
                            : of (x,y)
          - 'snapped'       : (data=float z coordinate) snap the layer to
                            : provided z coordinate, telescoping as needed
          - 'node'          : thickness provided on each node of the surface domain
          - 'cell'          : thickness provided on each cell of the surface domain,
                            : interpolate to nodes
    
        NOTE: dz is uniform through the layer in all but the 'snapped' case
        NOTE: 2D mesh is always labeled 'surface', extrusion is always downwards
        """

        # make the data all lists
        # ---------------------------------
        def is_list(data):
            if type(data) is str:
                return False
            try:
                len(data)
            except TypeError:
                return False
            else:
                return True
        
        if is_list(layer_types):
            if not is_list(layer_data):
                layer_data = [layer_data,]*len(layer_types)
            else:
                assert len(layer_data) == len(layer_types)

            if not is_list(ncells_per_layer):
                ncells_per_layer = [ncells_per_layer,]*len(layer_types)
            else:
                assert len(ncells_per_layer) == len(layer_types)

        elif is_list(layer_data):
            layer_type = [layer_type,]*len(layer_data)

            if not is_list(ncells_per_layer):
                ncells_per_layer = [ncells_per_layer,]*len(layer_data)
            else:
                assert len(ncells_per_layer) == len(layer_data)

        elif is_list(ncells_per_layer):
            layer_type = [layer_type,]*len(ncells_per_layer)
            layer_data = [layer_data,]*len(ncells_per_layer)
        else:
            layer_type = [layer_type,]
            layer_data = [layer_data,]
            ncells_per_layer = [ncells_per_layer,]
                
        # helper data and functions for mapping indices from 2D to 3D
        # ------------------------------------------------------------------
        ncells_tall = sum(ncells_per_layer)
        ncells_total = ncells_tall * mesh2D.num_cells()
        nfaces_total = (ncells_tall+1) * mesh2D.num_cells() + ncells_tall * mesh2D.num_edges()
        nnodes_total = (ncells_tall+1) * mesh2D.num_nodes()

        np_mat_ids = np.array(mat_ids, dtype=int)
        if np_mat_ids.size == np.size(np_mat_ids, 0):
            if np_mat_ids.size == 1:
                np_mat_ids = np.full((len(ncells_per_layer), mesh2D.num_cells()), mat_ids[0], dtype=int)
            else:
                np_mat_ids = np.empty((len(ncells_per_layer), mesh2D.num_cells()), dtype=int)
                for ilay in range(len(ncells_per_layer)):
                    np_mat_ids[ilay, :] = np.full(mesh2D.num_cells(), mat_ids[ilay], dtype=int)


        def col_to_id(column, z_cell):
            return z_cell + column * ncells_tall

        def node_to_id(node, z_node):
            return z_node + node * (ncells_tall+1)

        def edge_to_id(edge, z_cell):
            return (ncells_tall + 1) * mesh2D.num_cells() + z_cell + edge * ncells_tall

        # create coordinates
        # ---------------------------------
        coords = np.zeros((mesh2D.coords.shape[0],ncells_tall+1, 3),'d')
        coords[:,:,0:2] = np.expand_dims(mesh2D.coords[:,0:2],1)

        if mesh2D.dim == 3:
            coords[:,0,2] = mesh2D.coords[:,2]
        # else the surface is at 0 depth

        cell_layer_start = 0
        for layer_type, layer_datum, ncells in zip(layer_types, layer_data, ncells_per_layer):
            if layer_type.lower() == 'constant':
                dz = float(layer_datum) / ncells
                for i in range(1,ncells+1):
                    coords[:,cell_layer_start+i,2] = coords[:,cell_layer_start,2] - i * dz

            else:
                # allocate an array of coordinates for the bottom of the layer
                layer_bottom = np.zeros((mesh2D.coords.shape[0],),'d')

                if layer_type.lower() == 'snapped':
                    # layer bottom is uniform
                    layer_bottom[:] = layer_datum

                elif layer_type.lower() == 'function':
                    # layer thickness is given by a function evaluation of x,y
                    for node_col in range(mesh2D.coords.shape[0]):
                        layer_bottom[node_col] = coords[node_col,cell_layer_start,2] - layer_datum(coords[node_col,0,0], coords[node_col,0,1])

                elif layer_type.lower() == 'node':
                    # layer bottom specifically provided through thickness
                    layer_bottom[:] = coords[:,cell_layer_start,2] - layer_datum

                elif layer_type.lower() == 'cell':
                    # interpolate cell thicknesses to node thicknesses
                    import scipy.interpolate
                    centroids = mesh2D.cell_centroids()
                    interp = scipy.interpolate.interp2d(centroids[:,0], centroids[:,1], layer_datum, kind='linear')
                    layer_bottom[:] = coords[:,cell_layer_start,2] - interp(mesh2D.coords[:,0], mesh2D.coords[:,1])

                else:
                    raise RuntimeError("Unrecognized layer_type '%s'"%layer_type)

                # linspace from bottom of previous layer to bottom of this layer
                for node_col in range(mesh2D.coords.shape[0]):
                    coords[node_col,cell_layer_start:cell_layer_start+ncells+1,2] = np.linspace(coords[node_col,cell_layer_start,2], layer_bottom[node_col], ncells+1)
                
            cell_layer_start = cell_layer_start + ncells

        # create faces, face sets, cells
        bottom = []
        surface = []
        faces = []
        cells = [list() for c in range(ncells_total)]

        # -- loop over the columns, adding the horizontal faces
        for col in range(mesh2D.num_cells()):
            nodes_2 = mesh2D.conn[col]
            surface.append(col_to_id(col,0))
            for z_face in range(ncells_tall + 1):
                i_f = len(faces)
                f = [node_to_id(n, z_face) for n in nodes_2]

                if z_face != ncells_tall:
                    cells[col_to_id(col, z_face)].append(i_f)
                if z_face != 0:
                    cells[col_to_id(col, z_face-1)].append(i_f)

                faces.append(f)
            bottom.append(col_to_id(col,ncells_tall-1))

        # -- loop over the columns, adding the vertical faces
        added = dict()
        for col in range(mesh2D.num_cells()):
            nodes_2 = mesh2D.conn[col]
            for i in range(len(nodes_2)):
                edge = mesh2D.edge_hash(nodes_2[i], nodes_2[(i+1)%len(nodes_2)])
                try:
                    i_e = added[edge]
                except KeyError:
                    # faces not yet added to facelist
                    i_e = len(added.keys())
                    added[edge] = i_e
                    
                    for z_face in range(ncells_tall):
                        i_f = len(faces)
                        assert i_f == edge_to_id(i_e, z_face)
                        f = [node_to_id(edge[0], z_face),
                             node_to_id(edge[1], z_face),
                             node_to_id(edge[1], z_face+1),
                             node_to_id(edge[0], z_face+1)]
                        faces.append(f)
                        cells[col_to_id(col, z_face)].append(i_f)

                else:
                    # faces already added from previous column
                    for z_face in range(ncells_tall):
                        i_f = edge_to_id(i_e, z_face)
                        cells[col_to_id(col, z_face)].append(i_f)

        # Do some idiot checking
        assert len(faces) == nfaces_total
        for c in cells:
            assert len(c) > 4
        assert len(surface) == mesh2D.num_cells()
        assert len(bottom) == mesh2D.num_cells()

        # make the material ids
        material_ids = np.zeros((len(cells),),'i')
        for col in range(mesh2D.num_cells()):
            z_cell = 0
            for ilay in range(len(ncells_per_layer)):
                ncells = ncells_per_layer[ilay]
                for i in range(z_cell, z_cell+ncells):
                    material_ids[col_to_id(col, i)] = np_mat_ids[ilay, col]
                z_cell = z_cell + ncells

        # make the side sets
        side_sets = []
        side_sets.append(SideSet("bottom", 1, bottom, [1,]*len(bottom)))
        side_sets.append(SideSet("surface", 2, surface, [0,]*len(surface)))

        # reshape coords
        coords = coords.reshape(nnodes_total, 3)        
        
        for e,s in zip(side_sets[0].elem_list, side_sets[0].side_list):
            face = cells[e][s]
            fz_coords = np.array([coords[n] for n in faces[face]])
            #print "bottom centroid = ", np.mean(fz_coords, axis=0)

        for e,s in zip(side_sets[1].elem_list, side_sets[1].side_list):
            face = cells[e][s]
            fz_coords = np.array([coords[n] for n in faces[face]])
            #print "surface centroid = ", np.mean(fz_coords, axis=0)
        
        # instantiate the mesh
        return cls(coords, faces, cells, side_sets=side_sets, material_ids=material_ids)



def commandline_options():
    parser = argparse.ArgumentParser(description='Extrude a 2D mesh to make a 3D mesh')
    parser.add_argument("-n", "--num-cells", default=10, type=int,
                        help="number of cells to extrude")
    parser.add_argument("-d", "--depth", default=40.0, type=float,
                        help="depth to extrude")
    parser.add_argument("-o", "--outfile", default=None, type=str,
                        help="output filename")
    parser.add_argument("-p", "--plot", default=False, action="store_true",
                        help="plot the 2D mesh")
    parser.add_argument("infile",metavar="INFILE", type=str,
                        help="input filename of surface mesh")

    options = parser.parse_args()

    if options.outfile is None:
        options.outfile = ".".join(options.infile.split(".")[:-1])+".exo"
    

    if os.path.isfile(options.outfile):
        print 'Output file "%s" exists, cowardly not overwriting.'%options.outfile
        sys.exit(1)

    if not os.path.isfile(options.infile):
        print 'No input file provided'        
        parser.print_usage()
        sys.exit(1)

    return options
    


if __name__ == "__main__":
    options = commandline_options()
        
    m2 = Mesh2D.read_VTK(options.infile)
    if options.plot:
        m2.plot()
    m3 = Mesh3D.extruded_Mesh2D(m2, [options.depth,], [options.num_cells,], [10000,])
    m3.write_exodus(options.outfile)

