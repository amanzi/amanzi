/*
 * unscramble.cpp - Read in an Amanzi Restart file, returns an unpermutated
 *                  restart file
 *
 * Usage: unscramble meshfile.hdf5 vizfile.hdf5 new_filename.hdf5
 *
 * Created:  Erin Iesulauro Barker  Thursday May 3, 2012
 * Modified:
 *
 */

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include "hdf5.h"

#define FALSE 0

std::vector<std::string> datasetList;
std::vector<std::string> groupList;
int num_nodes;
int num_elems;

herr_t group_info(hid_t loc_id, const char *name, void *empty)
{

  H5G_stat_t buf;
  
  H5Gget_objinfo(loc_id, name, FALSE, &buf);
  switch (buf.type) {
    case H5G_GROUP:
      groupList.push_back(name);
      break;
    default:
      break;
  }
  return 0;
  
}

herr_t dataset_info(hid_t loc_id, const char *name, void *empty)
{
  
  H5G_stat_t buf;
  
  H5Gget_objinfo(loc_id, name, FALSE, &buf);
  switch (buf.type) {
    case H5G_DATASET:
      datasetList.push_back(name);
      break;
    default:
      break;
  }
  return 0;
  
}

herr_t unpermute(const char *name, hid_t file_id, hid_t new_fileid, int *nodemap, int *elemmap)
{
  hid_t dataset_id, dataspace;
  herr_t status;
  int rank, num;
  //hsize_t *cdims, *mdims;
  
  // open dataset
  // get dimensions
  // read data
  std::cout << "  E>> in upermute" << std::endl;
  dataset_id = H5Dopen(file_id, name, H5P_DEFAULT) ;
  dataspace = H5Dget_space(dataset_id) ;
  rank = H5Sget_simple_extent_ndims(dataspace) ;
  hsize_t *cdims = new hsize_t[rank];
  hsize_t *mdims = new hsize_t[rank];
  H5Sget_simple_extent_dims(dataspace, cdims, mdims) ;
  num = cdims[0] ;
  hid_t ds_type = H5Dget_type(dataset_id);
  hid_t ds_ntype = H5Tget_native_type(ds_type,H5T_DIR_DEFAULT);
  hid_t	ds_class = H5Tget_class(ds_type);
    
  // comparing class: options = H5T_INTEGER, H5T_FLOAT, H5T_STRING, 
  // or H5T_BITFIELD, H5T_OPAQUE, H5T_COMPOUND, H5T_REFERENCE, 
  // H5T_ENUM, H5T_VLEN, H5T_ARRAY
  if (ds_class == H5T_FLOAT) {
    double *data = new double[cdims[0]] ;
    double *new_data = new double[cdims[0]];
    std::cout << "    E>> doing FLOAT read" << std::endl;
    status = H5Dread(dataset_id, ds_ntype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) ;
    status = H5Dclose(dataset_id) ;
    
    // unpermute data
    std::cout << "    E>> compare "<<cdims[0]<<" to #nodes="<<num_nodes<<" or #elems="<<num_elems<< std::endl;
    if (cdims[0] == num_nodes) {
      std::cout << "    E>> unpermute node data" << std::endl;
      for (int i=0; i<num; ++i) {
        new_data[nodemap[i]] = data[i];
        std::cout << "      E>> "<<i<<" to " << nodemap[i]<< std::endl;
      }
    }
    else if (cdims[0] == num_elems) {
      std::cout << "    E>> unpermute elem data" << std::endl;
      for (int i=0; i<num; ++i) {
        new_data[elemmap[i]] = data[i];
        std::cout << "      E>> "<<i<<" to " << elemmap[i]<< std::endl;
      }
    }
    // write data
    dataspace = H5Screate_simple(2, cdims, NULL);
    dataset_id = H5Dcreate(new_fileid, name, ds_ntype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (cdims[0] == num_nodes || cdims[0] == num_elems) {
      std::cout << "    E>> write permuted data" << std::endl;
      status = H5Dwrite(dataset_id, ds_ntype, H5S_ALL, H5S_ALL, H5P_DEFAULT, new_data);
    } else {    
      std::cout << "    E>> write straight data" << std::endl;
      status = H5Dwrite(dataset_id, ds_ntype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    }
    delete[] data;
    delete[] new_data;
  } 

  else if (ds_class == H5T_INTEGER) {
    int *data = new int[cdims[0]] ;
    int *new_data = new int[cdims[0]];
    std::cout << "    E>> doing INT read" << std::endl;
    status = H5Dread(dataset_id, ds_ntype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data) ;
    status = H5Dclose(dataset_id) ;
    
    // unpermute data
    
    if (cdims[0] == num_nodes) {
      std::cout << "    E>> unpermute node data" << std::endl;
      for (int i=0; i<num; ++i) {
        new_data[nodemap[i]] = data[i];
      }
    }
    else if (cdims[0] == num_elems) {
      std::cout << "    E>> unpermute elem data" << std::endl;
      for (int i=0; i<num; ++i) {
        new_data[elemmap[i]] = data[i];
      }
    }
    // write data
    dataspace = H5Screate_simple(2, cdims, NULL);
    dataset_id = H5Dcreate(new_fileid, name, ds_ntype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (cdims[0] == num_nodes || cdims[0] == num_elems) {
      std::cout << "    E>> write permuted data" << std::endl;
      status = H5Dwrite(dataset_id, ds_ntype, H5S_ALL, H5S_ALL, H5P_DEFAULT, new_data);
    } else {    
      std::cout << "    E>> write straight data" << std::endl;
      status = H5Dwrite(dataset_id, ds_ntype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    }
    delete[] data;
    delete[] new_data;
  } 
  else if (ds_class == H5T_STRING) {
    // make string array thing
  } else {
    // error/exit because we can't handle type
  }
  
  H5Dclose(dataset_id);
  H5Sclose(dataspace);
  H5Tclose(ds_type);
  std::cout << "  E>> done in upermute" << std::endl;
  return 0;
}


int main (int argc, char *argv[])
{
  FILE *fp;
  hid_t mesh_file, dataset_id, dataspace, new_file, dt, data_file = 0;
  hsize_t *cdims, *mdims, dimsf[2];
  herr_t status;
  int rank;
  
  // open mesh file
  if ((fp = fopen(argv[1],"r")) == NULL)
  {
    printf("AMANZI:: Unable to open mesh file %s\n",argv[1]);
    // exit
  } else {
    fclose(fp);
    mesh_file = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
    
    // create new hdf5 file
    if ((fp = fopen(argv[3],"r")) == NULL)
    {
      new_file = H5Fcreate(argv[3], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } else {
      fclose(fp);
      printf("AMANZI:: Opening existing H5 file %s - will overwrite contents\n",argv[3]);
      new_file = H5Fopen(argv[3], H5F_ACC_RDWR, H5P_DEFAULT);
    }
    
    // read Mesh/NodeMap, store num_nodes
    std::cout << "E>> read NodeMap" << std::endl;
    dataset_id = H5Dopen(mesh_file, "/Mesh/NodeMap", H5P_DEFAULT) ;
    dataspace = H5Dget_space(dataset_id) ;
    rank = H5Sget_simple_extent_ndims(dataspace) ;
    cdims = (hsize_t*)malloc(rank*sizeof(hsize_t)) ;
    mdims = (hsize_t*)malloc(rank*sizeof(hsize_t)) ;
    H5Sget_simple_extent_dims(dataspace, cdims, mdims) ;
    num_nodes = cdims[0] ;
    int nodemap[num_nodes] ;
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nodemap) ;
    status = H5Dclose(dataset_id) ;
    
    // read Mesh/ElementMap, store num_elems
    std::cout << "E>> read ElementMap" << std::endl;
    dataset_id = H5Dopen(mesh_file, "/Mesh/ElementMap", H5P_DEFAULT) ;
    dataspace = H5Dget_space(dataset_id) ;
    rank = H5Sget_simple_extent_ndims(dataspace) ;
    cdims = (hsize_t*)malloc(rank*sizeof(hsize_t)) ;
    mdims = (hsize_t*)malloc(rank*sizeof(hsize_t)) ;
    H5Sget_simple_extent_dims(dataspace, cdims, mdims) ;
    num_elems = cdims[0] ;
    int elemmap[num_elems] ;
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, elemmap) ;
    status = H5Dclose(dataset_id) ;
    
    // write out unpermuted mesh
    
    // read Mesh
    std::cout << "E>> read Nodes" << std::endl;
    dataset_id = H5Dopen(mesh_file, "/Mesh/Nodes", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset_id);
    rank = H5Sget_simple_extent_ndims(dataspace) ;
    cdims = (hsize_t*)malloc(rank*sizeof(hsize_t)) ;
    mdims = (hsize_t*)malloc(rank*sizeof(hsize_t)) ;
    H5Sget_simple_extent_dims(dataspace, cdims, mdims) ;
    if (cdims[0] != num_nodes) {
      // need to throw an error here, something is wrong!
    }
    if (cdims[1] != 3) {
      // write a warning to the screen, but keep going
    }
    float nodes[cdims[0]][cdims[1]] ;
    std::cout << "  E>> read dims: " << cdims[0] << " x " << cdims[1] << std::endl;
    status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, nodes) ;
    status = H5Dclose(dataset_id) ;
    
    std::cout << "E>> read mixedelements" << std::endl;
    dataset_id = H5Dopen(mesh_file, "/Mesh/MixedElements", H5P_DEFAULT);
    dataspace = H5Dget_space(dataset_id);
    rank = H5Sget_simple_extent_ndims(dataspace) ;
    cdims = (hsize_t*)malloc(rank*sizeof(hsize_t)) ;
    mdims = (hsize_t*)malloc(rank*sizeof(hsize_t)) ;
    H5Sget_simple_extent_dims(dataspace, cdims, mdims) ;
    std::cout << "  E>> read dims: " << cdims[0] << " x " << cdims[1] << std::endl;
    int elem_len = cdims[0] ;
    int elems[elem_len] ;
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, elems) ;
    status = H5Dclose(dataset_id) ;
    
    // unpermute nodes
    std::cout << "E>> unpermute nodes (" <<num_nodes<<")"<< std::endl;
    float mapnodes[num_nodes][3];
    std::cout << "  E>> working node: ";
    for (int i=0; i<num_nodes; i++) {
      std::cout << "  " << i ;
      mapnodes[nodemap[i]][0] = nodes[i][0];
      mapnodes[nodemap[i]][1] = nodes[i][1];
      mapnodes[nodemap[i]][2] = nodes[i][2];
    }
    std::cout << std::endl;
    
    // unpermute mixed elements
    // loop through elements to get type list
    std::cout << "E>> collect element information" << std::endl;
    int elem_types[num_elems][3]; // 0)element type id 1)connectivity length 2)offset is element array
    int elem_cnt = 0;
    for (int i=0; i<num_elems; i++) {
      elem_types[i][0] = elems[elem_cnt];
      int conn_len;
      switch (elem_types[i][0]) {
        case 4:
          conn_len = 3; //TRI
          break;
        case 5:
          conn_len = 4; //QUAD
          break;
        case 6:
          conn_len = 4; //TET
          break;
        case 7:
          conn_len = 5; //PYRAMID
          break;
        case 8:
          conn_len = 6; //PRISM or WEDGE
          break;
        case 9:
          conn_len = 8; // HEX
          break;
        //case 3:
        //  conn_len = ????; //POLYGON or POLYHED
        //  break;
        //TODO(barker): deal with polygons and polyhedra
        default:
          conn_len = -1;
          break;
      }
      elem_types[i][1] = conn_len;
      elem_types[i][2] = elem_cnt;
      elem_cnt += conn_len + 1;      
      std::cout << "  E>> id = "<<i<<" type = "<<elem_types[i][0]<<" conn_len = "<<elem_types[i][1]<<" offset = "<<elem_types[i][2]<<std::endl;
    }
    // do element unpermute
    int map_offset = 0;
    int org_offset = 0;
    int mapelems[elem_len];
    std::cout << "E>> create reverse elemmap" << std::endl;
    int rev_elemmap[num_elems];
    for (int i=0; i<num_elems; i++) {
      rev_elemmap[elemmap[i]] = i;
      std::cout << "  E>> rev["<<elemmap[i]<<"] = " <<i<< std::endl;
    }
    std::cout << "E>> now unpermute the element connectivities" << std::endl;
    for (int i=0; i<num_elems; i++) {
      int id = rev_elemmap[i];
      org_offset = elem_types[id][2];
      mapelems[map_offset] = elem_types[id][0];
      std::cout << "  E>> elem "<<i<<" = "<<id<<" with conn_len = "<<elem_types[i][1]<<" and offset = "<<elem_types[i][2]<< std::endl;
      for (int j=1; j<elem_types[id][1]+1; j++) {
        std::cout << "    E>> getting elems[" << org_offset<<"+"<<j<<"] = "<<elems[org_offset+j] << "  to ["<<map_offset+j<<"]"<< std::endl;
        mapelems[map_offset+j] = nodemap[elems[org_offset+j]];
      }
      map_offset += elem_types[id][1]+1;
    }
    std::cout << "E>> new ordered element connectivities:";
    for (int i=0; i<elem_len; i++) {
      std::cout << " " << mapelems[i];
    }
    std::cout << std::endl;
    
    // create Mesh group
    hid_t group = H5Gcreate(new_file, "/Mesh",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // write out mesh
    std::cout << "E>> write out new nodes" << std::endl;
    dimsf[0] = num_nodes;
    dimsf[1] = 3;
    dt = H5Tcopy(H5T_NATIVE_FLOAT);
    dataspace = H5Screate_simple(rank, dimsf, NULL);
    dataset_id = H5Dcreate(new_file,"/Mesh/Nodes", dt, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mapnodes);
    
    std::cout << "E>> write out new elements" << std::endl;
    dimsf[0] = elem_len;
    dimsf[1] = 1;
    dt = H5Tcopy(H5T_NATIVE_INT);
    dataspace = H5Screate_simple(rank, dimsf, NULL);
    dataset_id = H5Dcreate(new_file,"/Mesh/MixedElements", dt, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mapelems);
    
  
    // open data file
    if ((fp = fopen(argv[2],"r")) == NULL)
    {
      printf("AMANZI:: Unable to open data file %s\n",argv[2]);
      // exit
    } else {
      fclose(fp);
      data_file = H5Fopen(argv[2], H5F_ACC_RDWR, H5P_DEFAULT);
      
      // iterate over groups, fill in list of names
      std::cout << "E>> iterate to get field names" << std::endl;
      H5Giterate(data_file, "/", NULL, group_info, NULL);
      for (int i=0; i<groupList.size(); ++i) {
        std::cout << " " << groupList[i];
      }
      std::cout << std::endl;
      for (int i=0; i<groupList.size(); ++i) {
        std::stringstream group_name;
        group_name << "/" << groupList[i];
        std::cout << "E>> creating " << group_name.str() << std::endl;

        group = H5Gcreate(new_file, group_name.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        datasetList.erase(datasetList.begin(), datasetList.end());
        H5Giterate(data_file, group_name.str().c_str(), NULL, dataset_info, NULL);
        for (int j=0; j<datasetList.size(); j++) {
          std::stringstream ds_name;
          ds_name << group_name.str() << "/" << datasetList[j];
          std::cout << "  E>> unpermute " << ds_name.str() << std::endl;
          status = unpermute(ds_name.str().c_str(), data_file, new_file, nodemap, elemmap);
        }
      }
    }
    
    // close restart file, new file
    std::cout << "E>> closing files" << std::endl;
    status = H5Fclose(mesh_file);
    status = H5Fclose(data_file);
    status = H5Fclose(new_file);
  }
}

