/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for a Field.

Field also stores some basic metadata for Vis, checkpointing, etc.
------------------------------------------------------------------------- */

#include <string>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "exodusII.h" 

#include "dbc.hh"
#include "errors.hh"
#include "CompositeVector.hh"
#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "Field.hh"
#include "Field_CompositeVector.hh"

namespace Amanzi {

Field_CompositeVector::Field_CompositeVector(std::string fieldname, std::string owner) :
    Field::Field(fieldname, owner), data_() {
  type_ = COMPOSITE_VECTOR_FIELD;
};

Field_CompositeVector::Field_CompositeVector(std::string fieldname, std::string owner,
        const std::vector<std::vector<std::string> >& subfield_names) :
    Field::Field(fieldname, owner),
    data_(),
    subfield_names_(subfield_names) {
  type_ = COMPOSITE_VECTOR_FIELD;
};

Field_CompositeVector::Field_CompositeVector(std::string fieldname, std::string owner,
                                           Teuchos::RCP<CompositeVector>& data) :
    Field::Field(fieldname, owner), data_(data) {
  type_ = COMPOSITE_VECTOR_FIELD;
};

// copy constructor:
Field_CompositeVector::Field_CompositeVector(const Field_CompositeVector& other) :
    Field::Field(other),
    subfield_names_(other.subfield_names_) {
  data_ = Teuchos::rcp(new CompositeVector(*other.data_));
};

// Virtual copy constructor
Teuchos::RCP<Field> Field_CompositeVector::Clone() const {
  return Teuchos::rcp(new Field_CompositeVector(*this));
}

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> Field_CompositeVector::Clone(std::string fieldname) const {
  Teuchos::RCP<Field_CompositeVector> other = Teuchos::rcp(new Field_CompositeVector(*this));
  other->fieldname_ = fieldname;
  return other;
};

// Virtual copy constructor with non-empty name
Teuchos::RCP<Field> Field_CompositeVector::Clone(std::string fieldname, std::string owner) const {
  Teuchos::RCP<Field_CompositeVector> other = Teuchos::rcp(new Field_CompositeVector(*this));
  other->fieldname_ = fieldname;
  other->owner_ = owner;
  return other;
};

// Create the data
void Field_CompositeVector::CreateData() {}

// write-access to the data
Teuchos::RCP<CompositeVector> Field_CompositeVector::GetFieldData() {
  return data_;
};

// Overwrite data by pointer, not copy
void Field_CompositeVector::SetData(const Teuchos::RCP<CompositeVector>& data) {
  data_ = data;
};

void Field_CompositeVector::SetData(const CompositeVector& data) {
  *data_ = data;
};

void Field_CompositeVector::Initialize(Teuchos::ParameterList& plist) {
  // Protect against unset names
  EnsureSubfieldNames_();

  bool checkpoint_io = plist.get<bool>("write checkpoint", true);
  set_io_checkpoint(checkpoint_io);

  bool vis_io = plist.get<bool>("write vis", true);
  set_io_vis(vis_io);


  // First try all initialization method which set the entire data structure.
  // ------ Try to set values from a restart file -----
  if (plist.isParameter("restart file")) {
    std::string filename = plist.get<std::string>("restart file");
    bool read = ReadCheckpoint_(filename);
    if (!read) {
      Errors::Message message;
      message << "Field: \"" << fieldname() << "\" is not available in restart file \"" << filename << "\"";
      Exceptions::amanzi_throw(message);
    }      
    set_initialized();
    return;
  }


  // ------ Try to set values from an file -----
  if (plist.isSublist("exodus file initialization")) {
    // data must be pre-initialized to zero in case Exodus file does not
    // provide all values.
    data_->PutScalar(0.0);

    Teuchos::ParameterList file_list = plist.sublist("exodus file initialization");
    ReadVariableFromExodusII_(file_list);
    set_initialized();
  }

  // Next try all partial initialization methods -- typically cells.
  // ------ Try to set cell values from a restart file -----
  if (plist.isParameter("cells from file")) {
    std::string filename = plist.get<std::string>("cells from file");
    ReadCellsFromCheckpoint_(filename);
    set_initialized();
  }

  // ------ Set values using a constant -----
  if (plist.isParameter("constant")) {
    double value = plist.get<double>("constant");
    data_->PutScalar(value);
    set_initialized();
  }
  if (plist.isParameter("value")) {
    double value = plist.get<double>("value");
    data_->PutScalar(value);
    set_initialized();
  }

  // ------ Set values from 1D solution -----
  if (plist.isSublist("initialize from 1D column")) {
    Teuchos::ParameterList& init_plist = plist.sublist("initialize from 1D column");
    InitializeFromColumn_(init_plist);
    set_initialized();
  }

  // ------ Set values using a function -----
  if (plist.isSublist("function")) {
    Teuchos::ParameterList func_plist = plist.sublist("function");
 
    // -- potential use of a mapping operator first -- 
    bool map_normal = plist.get<bool>("dot with normal", false); 
    if (map_normal) { 
      // map_normal take a vector and dots it with face normals 
      AMANZI_ASSERT(data_->NumComponents() == 1); // one comp 
      AMANZI_ASSERT(data_->HasComponent("face")); // is named face 
      AMANZI_ASSERT(data_->Location("face") == AmanziMesh::FACE); // is on face 
      AMANZI_ASSERT(data_->NumVectors("face") == 1);  // and is scalar 
 
      // create a vector on faces of the appropriate dimension 
      int dim = data_->Mesh()->space_dimension(); 

      CompositeVectorSpace cvs;
      cvs.SetMesh(data_->Mesh());
      cvs.SetComponent("face", AmanziMesh::FACE, dim);
      Teuchos::RCP<CompositeVector> vel_vec = Teuchos::rcp(new CompositeVector(cvs));

      // evaluate the full xD velocity function 
      Teuchos::RCP<Functions::CompositeVectorFunction> func = 
          Functions::CreateCompositeVectorFunction(func_plist, vel_vec->Map()); 
      func->Compute(0.0, vel_vec.ptr()); 
 
      // CV's map may differ from the regular mesh face map 
      const auto& fmap = *data_->Map().Map("face", true);
      unsigned int nfaces_owned = data_->Mesh() 
          ->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED); 
 
      Epetra_MultiVector& dat_f = *data_->ViewComponent("face",false); 
      const Epetra_MultiVector& vel_f = *vel_vec->ViewComponent("face",false); 
 
      AmanziGeometry::Point vel(dim); 
      for (unsigned int f=0; f!=nfaces_owned; ++f) { 
        AmanziGeometry::Point normal = data_->Mesh()->face_normal(f); 
        for (int i = 0; i < dim; ++i) vel[i] = vel_f[i][f]; 

        int g = fmap.FirstPointInElement(f);
        for (int i = 0; i < fmap.ElementSize(f); ++i) {
          dat_f[0][g + i] = vel * normal; 
        }
      } 
      set_initialized(); 
 
    } else { 
      // no map, just evaluate the function 
      Teuchos::RCP<Functions::CompositeVectorFunction> func = 
          Functions::CreateCompositeVectorFunction(func_plist, data_->Map()); 
      func->Compute(0.0, data_.ptr()); 
      set_initialized(); 
    }
  }

  // ------ Set face values by interpolation -----
  if ((data_->HasComponent("face") || data_->HasComponent("boundary_face"))  && data_->HasComponent("cell") &&
      plist.get<bool>("initialize faces from cells", false)) {
    DeriveFaceValuesFromCellValues(*data_);
  }

  return;
};


void Field_CompositeVector::WriteVis(Visualization& vis) {
  Key name = vis.name();
  if (name == "domain") name = "";
  if (io_vis_ &&
      ((name == Keys::getDomain(fieldname_)) ||
       (Keys::isDomainSet(fieldname_) && name == Keys::getDomainSetName(fieldname_)))) {
    EnsureSubfieldNames_();

    // loop over the components and dump them to the vis file if possible
    int i = 0;
    for (CompositeVector::name_iterator compname=data_->begin();
         compname!=data_->end(); ++compname) {
      // check that this vector is a cell vector (currently this is the only
      // type of vector we can visualize
      if (data_->Location(*compname) == AmanziMesh::CELL) {
        // get the MultiVector that should be dumped
        Teuchos::RCP<Epetra_MultiVector> v = data_->ViewComponent(*compname, false);

        // construct the name for vis
        std::vector< std::string > vis_names(subfield_names_[i]);
        for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
          vis_names[j] = fieldname_ + std::string(".") + *compname
            + std::string(".") + subfield_names_[i][j];
        }
        vis.WriteVector(*v, vis_names);
      }
      i++;
    }
  }
};


void Field_CompositeVector::WriteCheckpoint(Checkpoint& chk) {
  if (io_checkpoint_) {
    EnsureSubfieldNames_();

    // loop over the components and dump them to the checkpoint file if possible
    int i = 0;
    for (CompositeVector::name_iterator compname=data_->begin();
         compname!=data_->end(); ++compname) {
      // get the MultiVector that should be dumped
      Teuchos::RCP<Epetra_MultiVector> v = data_->ViewComponent(*compname, false);

      // construct name for the field in the checkpoint
      std::vector<std::string> chkp_names(subfield_names_[i]);
      for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
        chkp_names[j] = fieldname_ + "." + *compname + "." + subfield_names_[i][j];
      }
      chk.WriteVector(*v, chkp_names);
      i++;
    }
  }
};


void Field_CompositeVector::ReadCellsFromCheckpoint_(const std::string& filename) {
  Teuchos::RCP<Amanzi::HDF5_MPI> file_input =
      Teuchos::rcp(new Amanzi::HDF5_MPI(data_->Comm(), filename));
  file_input->open_h5file();
  EnsureSubfieldNames_();

  int i = 0;
  for (CompositeVector::name_iterator compname=data_->begin();
       compname!=data_->end(); ++compname) {
    if (*compname == std::string("cell")) {

      // get the MultiVector that should be read
      Teuchos::RCP<Epetra_MultiVector> vec = data_->ViewComponent(*compname, false);
      AMANZI_ASSERT(vec != Teuchos::null);
      
      // construct name for the field in the checkpoint file
      for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
        std::string chkp_name = fieldname_ + "." + *compname + "." + subfield_names_[i][j];
        file_input->readData(*(*vec)(j), chkp_name);
      }
    }
    ++i;
  }

  file_input->close_h5file();
}


bool Field_CompositeVector::ReadCheckpoint_(std::string filename) {
  Amanzi::HDF5_MPI file_input(data_->Comm(), filename);
  file_input.open_h5file();
  bool read_complete = ReadCheckpoint(file_input);
  file_input.close_h5file();
  return read_complete;
}


// modify methods
// -- set data from file
bool Field_CompositeVector::ReadCheckpoint(HDF5_MPI& file_input) {
  EnsureSubfieldNames_();

  // loop over the components and read them from the checkpoint file if possible
  int i = 0;
  for (CompositeVector::name_iterator compname=data_->begin();
       compname!=data_->end(); ++compname) {
    // get the MultiVector that should be read
    Teuchos::RCP<Epetra_MultiVector> vec = data_->ViewComponent(*compname, false);

    // construct name for the field in the checkpoint file
    std::vector<std::string> chkp_names(subfield_names_[i]);
    for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
      chkp_names[j] = fieldname_ + "." + *compname + "." + subfield_names_[i][j];
    }
    for (unsigned int j = 0; j!=subfield_names_[i].size(); ++j) {
      bool read = file_input.readData(*(*vec)(j), chkp_names[j]);
      if (!read) 
        return false;     
    }
    i++;
  }
  return true;
}


void Field_CompositeVector::InitializeFromColumn_(Teuchos::ParameterList& plist) {
  // get filename, data names
  if (!plist.isParameter("file")) {
    Errors::Message message("Missing InitializeFromColumn parameter \"file\"");
    Exceptions::amanzi_throw(message);
  }
  std::string filename = plist.get<std::string>("file");
  std::string z_str = plist.get<std::string>("z header", "/z");
  std::string f_str = std::string("/")+fieldname_;
  f_str = plist.get<std::string>("f header", f_str);

  // Create the function
  Teuchos::ParameterList func_list;
  Teuchos::ParameterList& func_sublist = func_list.sublist("function-tabular");
  func_sublist.set("file", filename);
  func_sublist.set("x header", z_str);
  func_sublist.set("y header", f_str);
  FunctionFactory fac;
  Teuchos::RCP<Function> func = Teuchos::rcp(fac.Create(func_list));
      
  // orientation
  std::string orientation = plist.get<std::string>("coordinate orientation", "standard");
  if (orientation != "standard" && orientation != "depth") {
    Errors::Message message("InitializeFromColumn parameter \"orientation\" must be either \"standard\" (bottom up) or \"depth\" (top down)");
    Exceptions::amanzi_throw(message);
  }

  // starting surface
  Teuchos::Array<std::string> sidesets;
  if (plist.isParameter("surface sideset")) {
    sidesets.push_back(plist.get<std::string>("surface sideset"));
  } else if (plist.isParameter("surface sidesets")) {
    sidesets = plist.get<Teuchos::Array<std::string> >("surface sidesets");
  } else {
    Errors::Message message("Missing InitializeFromColumn parameter \"surface sideset\" or \"surface sidesets\"");
    Exceptions::amanzi_throw(message);
  }

  data_->Mesh()->build_columns();
  // evaluate
  Epetra_MultiVector& vec = *data_->ViewComponent("cell",false);
  if (orientation == "depth") {
    double z0;
    std::vector<double> z(1);

    AmanziMesh::Entity_ID_List surf_faces; 
    for (Teuchos::Array<std::string>::const_iterator setname=sidesets.begin();
         setname!=sidesets.end(); ++setname) {
      data_->Mesh()->get_set_entities(*setname,AmanziMesh::FACE,
              AmanziMesh::Parallel_type::OWNED, &surf_faces);

      for (AmanziMesh::Entity_ID_List::const_iterator f=surf_faces.begin();
           f!=surf_faces.end(); ++f) {
        // Collect the reference coordinate z0
        AmanziGeometry::Point x0 = data_->Mesh()->face_centroid(*f);
        z0 = x0[x0.dim()-1];

        // Iterate down the column
        AmanziMesh::Entity_ID_List cells;
        data_->Mesh()->face_get_cells(*f, AmanziMesh::Parallel_type::OWNED, &cells);
        AMANZI_ASSERT(cells.size() == 1);
        AmanziMesh::Entity_ID c = cells[0];

        while (c >= 0) {
          AmanziGeometry::Point x1 = data_->Mesh()->cell_centroid(c);
          z[0] = z0 - x1[x1.dim()-1];
          vec[0][c] = (*func)(z);
          c = data_->Mesh()->cell_get_cell_below(c);
        }
      }
    }

  } else {
    double z0;
    std::vector<double> z(1);

    AmanziMesh::Entity_ID_List surf_faces;
    for (Teuchos::Array<std::string>::const_iterator setname=sidesets.begin();
         setname!=sidesets.end(); ++setname) {
      data_->Mesh()->get_set_entities(*setname,AmanziMesh::FACE,
              AmanziMesh::Parallel_type::OWNED, &surf_faces);

      for (AmanziMesh::Entity_ID_List::const_iterator f=surf_faces.begin();
           f!=surf_faces.end(); ++f) {
        // Collect the reference coordinate z0
        AmanziGeometry::Point x0 = data_->Mesh()->face_centroid(*f);
        z0 = x0[x0.dim()-1];

        // Iterate down the column
        AmanziMesh::Entity_ID_List cells;
        data_->Mesh()->face_get_cells(*f, AmanziMesh::Parallel_type::OWNED, &cells);
        AMANZI_ASSERT(cells.size() == 1);
        AmanziMesh::Entity_ID c = cells[0];

        while (c >= 0) {
          AmanziGeometry::Point x1 = data_->Mesh()->cell_centroid(c);
          z[0] = x1[x1.dim()-1] - z0;
          vec[0][c] = (*func)(z);
          c = data_->Mesh()->cell_get_cell_above(c);
        }
      }
    }
  }
}


void Field_CompositeVector::EnsureSubfieldNames_() {
  // set default values for subfield names, ensuring they are unique
  if (subfield_names_.size() == 0) {
    subfield_names_.resize(data_->NumComponents());

    unsigned int i = 0;
    for (CompositeVector::name_iterator compname=data_->begin();
         compname!=data_->end(); ++compname) {
      subfield_names_[i].resize(data_->NumVectors(*compname));

      for (unsigned int j=0; j!=subfield_names_[i].size(); ++j) {
        std::stringstream s;
        s << j;
        subfield_names_[i][j] = s.str();
      }
      ++i;
    }
  } else {

    if (subfield_names_.size() != data_->NumComponents()) {
      subfield_names_.resize(data_->NumComponents());
    }

    unsigned int i = 0;
    for (CompositeVector::name_iterator compname=data_->begin();
         compname!=data_->end(); ++compname) {
      
      if (subfield_names_[i].size() == 0) {
        subfield_names_[i].resize(data_->NumVectors(*compname));
      }

      for (int j=0; j!=subfield_names_[i].size(); ++j) {
        if (subfield_names_[i][j].length() == 0) {
          if (i==0) {
            std::stringstream s;
            s << j;
            subfield_names_[i][j] = s.str();
          }else{
            subfield_names_[i][j] = subfield_names_[0][j];
          }
        }
      }
      i++;      
    }
  }
};


// -----------------------------------------------------------------------------
// New routine for reading cell-based varibles as attributes.
// It recongnizes parallel and serial inputs.
// -----------------------------------------------------------------------------
void Field_CompositeVector::ReadVariableFromExodusII_(Teuchos::ParameterList& file_list) 
{ 
  Epetra_MultiVector& dat_f = *data_->ViewComponent("cell", false); 
  int nvectors = dat_f.NumVectors(); 
 
  std::string file_name = file_list.get<std::string>("file"); 
  std::vector<std::string> attributes = 
      file_list.get<Teuchos::Array<std::string> >("attributes").toVector(); 
 
  // open ExodusII file 
  auto comm = data_->Comm();
 
  if (comm->NumProc() > 1) { 
    int ndigits = (int)floor(log10(comm->NumProc())) + 1;
    std::string fmt = boost::str(boost::format("%%s.%%d.%%0%dd") % ndigits);
    file_name = boost::str(boost::format(fmt) % file_name % comm->NumProc() % comm->MyPID());
  } 
 
  int CPU_word_size(8), IO_word_size(0), ierr; 
  float version; 
  int exoid = ex_open(file_name.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version); 
  if (comm->MyPID() == 0) {
    printf("Trying file: %s ws=%d %d  id=%d\n", file_name.c_str(), CPU_word_size, IO_word_size, exoid); 
  }

  // check if we have to use serial file
  int fail = (exoid < 0) ? 1 : 0;
  int fail_tmp(fail);
  bool distributed_data(true);

  comm->SumAll(&fail_tmp, &fail, 1);
  if (fail == comm->NumProc()) {
    Errors::Message msg("Rao is working on new data layout which we need to proceed.");
    Exceptions::amanzi_throw(msg);

    file_name = file_list.get<std::string>("file"); 
    distributed_data = false;
    if (comm->MyPID() == 0) {
      exoid = ex_open(file_name.c_str(), EX_READ, &CPU_word_size, &IO_word_size, &version); 
      printf("Opening file: %s ws=%d %d  id=%d\n", file_name.c_str(), CPU_word_size, IO_word_size, exoid); 
    }
  } else if (fail > 0) {
    Errors::Message msg("A few parallel Exodus files are missing, but not all.");
    Exceptions::amanzi_throw(msg);
  }
 
  // read database parameters 
  if (comm->MyPID() == 0 || distributed_data) {  
    int dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets; 
    char title[MAX_LINE_LENGTH + 1]; 
    ierr = ex_get_init(exoid, title, &dim, &num_nodes, &num_elem, 
                       &num_elem_blk, &num_node_sets, &num_side_sets); 
 
    int* ids = (int*) calloc(num_elem_blk, sizeof(int)); 
    ierr = ex_get_ids(exoid, EX_ELEM_BLOCK, ids); 
 
    // read number of variables 
    int num_vars;
    auto obj_type = ex_var_type_to_ex_entity_type('e');
    ierr = ex_get_variable_param(exoid, obj_type, &num_vars);
    if (ierr < 0) printf("Exodus file has no variables.\n");

    char* var_names[num_vars];
    for (int i = 0; i < num_vars; i++) {
      var_names[i] = (char*) calloc ((MAX_STR_LENGTH+1), sizeof(char));
    }

    obj_type = ex_var_type_to_ex_entity_type('e');
    ierr = ex_get_variable_names(exoid, obj_type, num_vars, var_names);
    if (ierr < 0) printf("Exodus file cannot read variable names.\n");

    int var_index(-1), ncells;
    for (int k = 0; k < nvectors; ++k) {
      for (int i = 0; i < num_vars; i++) {
        std::string tmp(var_names[i]);
        if (tmp == attributes[k]) var_index = i + 1;
      }
      if (var_index < 0) printf("Exodus file has no variable \"%s\".\n", attributes[k].c_str());
      printf("Variable \"%s\" has index %d.\n", attributes[k].c_str(), var_index);

      // read variable with the k-th attribute
      int offset = 0; 
      char elem_type[MAX_LINE_LENGTH + 1]; 
      for (int i = 0; i < num_elem_blk; i++) { 
        int num_elem_this_blk, num_attr, num_nodes_elem; 
        ierr = ex_get_block(exoid, EX_ELEM_BLOCK, ids[i], elem_type, &num_elem_this_blk, 
                            &num_nodes_elem, 0, 0, &num_attr); 
 
        double* var_values = (double*) calloc(num_elem_this_blk, sizeof(double)); 
        ierr = ex_get_var(exoid, 1, EX_ELEM_BLOCK, var_index, ids[i], num_elem_this_blk, var_values); 
 
        for (int n = 0; n < num_elem_this_blk; n++) { 
          int c = n + offset; 
          dat_f[k][c] = var_values[n]; 
        } 
        free(var_values); 
        printf("MyPID=%d  ierr=%d  id=%d  ncells=%d\n", comm->MyPID(), ierr, ids[i], num_elem_this_blk); 
 
        offset += num_elem_this_blk; 
      } 
      ncells = offset; 
    }

    for (int i = 0; i < num_vars; i++) {
      free(var_names[i]);
    }
 
    ierr = ex_close(exoid); 
    printf("Closing file: %s ncells=%d error=%d\n", file_name.c_str(), ncells, ierr); 
  }
} 


long int Field_CompositeVector::GetLocalElementCount() {
  long int count(0);
  for (CompositeVector::name_iterator compname=data_->begin();
       compname!=data_->end(); ++compname) {
    // get the MultiVector that should be dumped
    Teuchos::RCP<Epetra_MultiVector> v = data_->ViewComponent(*compname, false);
    count += v->NumVectors() * v->MyLength();
  }
  return count;
}

// void Field_CompositeVector::CopyFace2BndFace(const Epetra_MultiVector& vec_face, Epetra_MultiVector& vec_bndface) {
//   for (int bf=0   
// }
  
} // namespace Amanzi
