/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef _GMVWRITEH_
#define _GMVWRITEH_
#endif


void gmvwrite_openfile(char filenam[]);

void gmvwrite_openfile_ascii(char filenam[]);

void gmvwrite_openfile_ir(char filenam[], int isize, int rsize);

void gmvwrite_openfile_cxir(char filenam[], int isize, int rsize);

void gmvwrite_openfile_ir_ascii(char filenam[], int isize, int rsize);

void gmvwrite_closefile(void);

void gmvwrite_nodes_fromfile(char *filename, long nndes);

void gmvwrite_node_data(void *nndes, void *x, void *y, void *z);

void gmvwrite_node_data_struct(void *nxv, void *nyv, void *nzv,
			       void *x,  void *y, void *z);

void gmvwrite_node_data_lstruct(void *nxv, void *nyv, void *nzv,
			        void *x, void *y, void *z);

void gmvwrite_node_data_amr(int nxc, int nyc, int nzc, void *x0, void *y0,
		            void *z0, void *dx, void *dy, void *dz);

void gmvwrite_nodev_fromfile(char *filename, long nndes);

void gmvwrite_nodev_data(void *nndes, void *x, void *y, void *z);

void gmvwrite_nodev_data_lstruct(void *nxv, void *nyv, void *nzv,
			        void *x, void *y, void *z);

void gmvwrite_cells_amr(void *numcells, void *numtop, void *daughters);

void gmvwrite_cells_fromfile(char *filename, long nclls);

void gmvwrite_cell_header(void *ncells);

void gmvwrite_cell_type(const char cell_type[], int nverts, void *nodes);

void gmvwrite_general_cell_type(char cell_type[], int nverts[], int nfaces,
			        void *nodeids);

void gmvwrite_faces_fromfile(char *filename, long nclls);

void gmvwrite_face_header(void *nfaces, void *ncells);

void gmvwrite_face_data(int nverts, void *nodeids, void *cellid1,
                   void *cellid2);

void gmvwrite_vfaces_fromfile(char *filename, long nclls);

void gmvwrite_vface_header(void *nfaces);

void gmvwrite_vface_data(int nverts, int facepe, void *oppface,
                         int oppfacepe, void *cellid, void *nodeids);

void gmvwrite_xfaces_fromfile(char *filename, long nfces, long nclls);

void gmvwrite_xface_header(void *nfaces);

void gmvwrite_xface_data(long totverts, void *nverts, void *nodeids,
                         void *cellid, void *oppface, void *facepe,
                         void *oppfacepe);

void gmvwrite_material_fromfile(char *filename);

void gmvwrite_material_header(int nmats, int data_type);

void gmvwrite_material_name(char matname[]);

void gmvwrite_material_ids(int matids[], int data_type);

void gmvwrite_velocity_data(int data_type, void *u, void *v, void *w);

void gmvwrite_variable_header(void);

void gmvwrite_variable_name_data(int data_type, char varname[], void *vids);

void gmvwrite_variable_endvars(void);

void gmvwrite_flag_fromfile(char *filename);

void gmvwrite_flag_header(void);

void gmvwrite_flag_name(char flagname[], int numtypes, int data_type);

void gmvwrite_flag_subname(const char subname[]);

void gmvwrite_flag_data(int data_type, int flag_data[]);

void gmvwrite_flag_endflag(void);

void gmvwrite_polygons_fromfile(char *filename);

void gmvwrite_polygons_header(void);

void gmvwrite_polygons_data(int nverts, int matnum, void *x, void *y, void *z);

void gmvwrite_polygons_endpoly(void);

void gmvwrite_tracers_header(int ntracers, void *x, void *y, void *z);

void gmvwrite_tracers_name_data(int ntracers, char tracername[], void *data);

void gmvwrite_tracers_endtrace(void);

void gmvwrite_probtime(double ptime);

void gmvwrite_cycleno(int cyclenum);

void gmvwrite_nodeids(void *nodeids);

void gmvwrite_nodeids_fromfile(char *filename);

void gmvwrite_cellids(void *cellids);

void gmvwrite_cellids_fromfile(char *filename);

void gmvwrite_surface_header(void *nsurf);

void gmvwrite_surface_data(int nverts, void *nodeids);

void gmvwrite_surface_fromfile(char *filename, long nsrf);

void gmvwrite_surfmats(int matids[]);

void gmvwrite_surfvel( void *u, void *v, void *w);

void gmvwrite_surfvars_header(void);

void gmvwrite_surfvars_name_data(char varname[], void *vids);

void gmvwrite_surfvars_endsvar(void);

void gmvwrite_surfflag_header(void);

void gmvwrite_surfflag_name(const char flagname[], int numtypes);

void gmvwrite_surfflag_subname(const char subname[]);

void gmvwrite_surfflag_data(int flag_data[]);

void gmvwrite_surfflag_endflag(void);

void gmvwrite_units_fromfile(char *filename);

void gmvwrite_units_header(void);

void gmvwrite_units_typehdr(int data_type, int numtypes);

void gmvwrite_units_name(char fldname[], char unitname[]);

void gmvwrite_units_endunit(void);

void gmvwrite_vinfo_header(void);

void gmvwrite_vinfo_name_data(int nelem, int nlines, char varname[], void *vids);

void gmvwrite_vinfo_endvinfo(void);

void gmvwrite_traceids(int ntracers, void *traceids);

void gmvwrite_traceids_fromfile(char *filename);

void gmvwrite_faceids(void *faceids);

void gmvwrite_faceids_fromfile(char *filename);

void gmvwrite_group_fromfile(char *filename);

void gmvwrite_group_header(void);

void gmvwrite_group_data(char groupname[], int data_type, int numgrp,
     void *group_data);

void gmvwrite_group_endgroup(void);

void gmvwrite_surfids(void *surfids);

void gmvwrite_codename(char codename[]);

void gmvwrite_codever(char codever[]);

void gmvwrite_simdate(char simdate[]);

void gmvwrite_subvars_header(void);

void gmvwrite_suvvars_name_data(int data_type, int numelem, char varname[],
                              void *vids, void *vdata);

void gmvwrite_subvars_endsubv(void);

void gmvwrite_ghosts(int data_type, int numghst, void *ghost_data);
