"""
Extrudes a VTK 2D mesh to generate an ExodusII 3D mesh.

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



Requires building Ethan's hacked exodus python wrappers:
------------------------------------------------------------

first, edit $ATS_SRC_DIR/tools/meshing_ats/cmake-script to fit your
system

then, download exodus-6.09
unzip/untar

$> cd exodus-6.09
$> patch -p1 < $ATS_SRC_DIR/tools/meshing_ats/exodus-6.09.patch

$> cd exodus
$> mkdir build
$> cd build
$> . $ATS_SRC_DIR/tools/meshing_ats/cmake-script
$> make
$> make install


"""

import sys,os
import numpy as np
import collections
sys.path.append("/Users/ecoon/research/coastal/exodus-6.09/exodus/install/python/")
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

    def plot(self):
        import colors
        cm = colors.cm_mapper(0,self.num_cells()-1)
        colors = [cm(i) for i in range(self.num_cells())]

        verts = [[self.coords[i,0:2] for i in f] for f in self.conn]
        from matplotlib import collections
        gons = collections.PolyCollection(verts, facecolors=colors)
        from matplotlib import pyplot as plt
        fig,ax = plt.subplots(1,1)
        ax.add_collection(gons)
        ax.autoscale_view()
        plt.show()
        


    @classmethod
    def read_VTK(cls, filename):
        with open(filename,'r') as fid:
            line = fid.readline()
            while not line.startswith("POINTS"):
                line = fid.readline()
            ncoords = int(line.strip().split()[1])
            coords = np.zeros([ncoords,3],'d')
            for i in range(ncoords):
                coords[i,:] = np.array([float(p) for p in fid.readline().strip().split()])

            line = fid.readline()
            while not line.startswith("POLYGONS"):
                line = fid.readline()
            ngons = int(line.strip().split()[1])
            gons = []
            for i in range(ngons):
                line = [int(n) for n in fid.readline().strip().split()[1:]]
                gon = [int(n) for n in line]

                # check handedness -- need normals to point up!
                d31 = (coords[gon[-1]] - coords[gon[0]])
                d21 = (coords[gon[1]] - coords[gon[0]])
                d31[2] = 0
                d21[2] = 0
                up = np.cross(d21,d31)
                if (up[2] < 0):
                    gon.reverse()
                    
                gons.append(gon)

        return cls(coords, gons)
            
            
                

    

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


    def write_exodus(self, filename):
        """Write the 3D mesh to ExodusII using arbitrary polyhedra spec"""
        # open the mesh file
        ep = exodus.ex_init_params(title=filename,
                                   num_dim=3,
                                   num_nodes=self.num_nodes(),
                                   num_face=self.num_faces(),
                                   num_face_blk=1,
                                   num_elem=self.num_cells(),
                                   num_elem_blk=len(self.material_id_list),
                                   num_side_sets=len(self.side_sets))
        e = exodus.exodus(filename, mode='w', array_type='numpy', init_params=ep)

        # put the coordinates
        e.put_coord_names(['coordX', 'coordY', 'coordZ'])
        e.put_coords(self.coords[:,0], self.coords[:,1], self.coords[:,2])

        # put the face block
        face_raveled = [i for f in self.face_to_node_conn for i in f]
        e.put_polyhedra_face_blk(1, self.num_faces(), len(face_raveled), 0)
        e.put_node_count_per_face(1, np.array([len(f) for f in self.face_to_node_conn]))
        e.put_face_node_conn(1, np.array(face_raveled)+1)

        # put cells in with blocks, which renumbers the cells, so we have to track sidesets.
        # Therefore we keep a map of old cell to new cell ordering
        new_to_old = []
        for i_m,m_id in enumerate(self.material_id_list):
            elems = [(i,c) for (i,c) in enumerate(self.elem_to_face_conn) if self.material_ids[i] == m_id]
            elems_raveled = [f for (i,c) in elems for f in c]
            new_to_old.extend([i for (i,c) in elems])

            e.put_polyhedra_elem_blk(m_id, len(elems), len(elems_raveled), 0)
            e.put_elem_blk_name(m_id, "MATERIAL_ID_%d"%m_id)
            e.put_face_count_per_polyhedra(m_id, np.array([len(c) for i,c in elems]))
            e.put_elem_face_conn(m_id, np.array(elems_raveled)+1)

        # add sidesets
        old_to_new = sorted([(old,i) for (i,old) in enumerate(new_to_old)], lambda a,b: int.__cmp__(a[0],b[0]))
        
        e.put_side_set_names([ss.name for ss in self.side_sets])
        for ss in self.side_sets:

            for elem in ss.elem_list:
                assert old_to_new[elem][0] == elem
            new_elem_list = [old_to_new[elem][1] for elem in ss.elem_list]                
            e.put_side_set_params(ss.setid, len(ss.elem_list), 0)
            e.put_side_set(ss.setid, np.array(new_elem_list)+1, np.array(ss.side_list)+1)

        # finish and close
        e.close()
        

    @classmethod
    def extruded_Mesh2D(cls, mesh2D, layer_thicknesses, ncells_per_layer, mat_ids, snap_bottom=False):
        """
        Regularly extrude a 2D mesh to make a 3D mesh.

        mesh2D              : a Mesh2D object
        layer_thicknesses   : an array/list of layer thicknesses
        ncells_per_layer    : an array/list of number of cells in the layer
        mat_ids             : an array/list of IDs for each layer

        NOTE: dz = layer_thickness[i] / ncells_per_layer[i]
        NOTE: 2D mesh is always labeled 'surface', extrusion is always downwards
        """
        # DBC
        assert len(layer_thicknesses) == len(ncells_per_layer)

        # helpers
        ncells_tall = sum(ncells_per_layer)
        dzs = [thickness/ncells for (thickness,ncells) in zip(layer_thicknesses, ncells_per_layer)]
        ncells_total = ncells_tall * mesh2D.num_cells()
        nfaces_total = (ncells_tall+1) * mesh2D.num_cells() + ncells_tall * mesh2D.num_edges()
        nnodes_total = (ncells_tall+1) * mesh2D.num_nodes()

        def col_to_id(column, z_cell):
            return z_cell + column * ncells_tall

        def node_to_id(node, z_node):
            return z_node + node * (ncells_tall+1)

        def edge_to_id(edge, z_cell):
            return (ncells_tall + 1) * mesh2D.num_cells() + z_cell + edge * ncells_tall

        # create coordinates
        coords = np.zeros((mesh2D.coords.shape[0],ncells_tall+1, 3),'d')
        coords[:,:,0:2] = np.expand_dims(mesh2D.coords[:,0:2],1)

        if mesh2D.dim == 3:
            coords[:,0,2] = mesh2D.coords[:,2]

        if snap_bottom:
            # do all but the bottom layer
            cell = 0
            for ncells, dz in zip(ncells_per_layer[:-1], dzs[:-1]):
                for i in range(ncells):
                    coords[:,cell+1,2] = coords[:,cell,2] - dz
                    cell = cell + 1

            # do the bottom layer as a variable spacing to match the bottom
            for node_col in range(mesh2D.coords.shape[0]):
                coords[node_col,cell:,2] = np.linspace(coords[node_col,cell,2], layer_thicknesses[-1], ncells_per_layer[-1]+1)

        else:
            cell = 0
            for ncells, dz in zip(ncells_per_layer[:-1], dzs[:-1]):
                for i in range(ncells):
                    coords[:,cell+1,2] = coords[:,cell,2] - dz
                    cell = cell + 1
                    


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
            for ncells, m_id in zip(ncells_per_layer, mat_ids):
                for i in range(z_cell, z_cell+ncells):
                    material_ids[col_to_id(col, i)] = m_id
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
            print "bottom centroid = ", np.mean(fz_coords, axis=0)

        for e,s in zip(side_sets[1].elem_list, side_sets[1].side_list):
            face = cells[e][s]
            fz_coords = np.array([coords[n] for n in faces[face]])
            print "surface centroid = ", np.mean(fz_coords, axis=0)
        
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

