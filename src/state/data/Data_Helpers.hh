/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License: BSD, see $AMANZI_DIR/COPYRIGHT
Author: Ethan Coon

Helpers that know how to read/write/etc data.

------------------------------------------------------------------------- */

#ifndef AMANZI_DATA_HELPERS_HH_
#define AMANZI_DATA_HELPERS_HH_

#include "Teuchos_ParameterList.hpp"

#include "errors.hh"

#include "CompositeVector.hh"
#include "Op.hh"

#include "FunctionFactory.hh"

#include "StateDefs.hh"
#include "Visualization.hh"
#include "Checkpoint.hh"



namespace Amanzi {
namespace Helpers {



//
// Visualization
// ======================================================================

// Default simply dispatches to Vis.  This fails to compile! if not either
// specifically implemented in Visualization class or specialized below.
template<typename T>
void
WriteVis(const Visualization& vis, const Key& fieldname,
         const std::vector<std::string>& subfieldnames,
         const T& t)
{
  UserWriteVis(vis, fieldname, subfieldnames, t);
}


//
// Checkpoint
// ======================================================================

// WriteCheckpoint writes data to file
//
// Default simply dispatches to Checkpoint.  This fails to compile! if not either
// specifically implemented in Checkpoint class or specialized below.
template<typename T>
void
WriteCheckpoint(const Checkpoint& chkp, const Key& fieldname,
                const T& t)
{
  UserWriteCheckpoint(chkp, fieldname, t);
}


// ReadCheckpoint reads data from file
//
// Default simply dispatches to Checkpoint.  This fails to compile! if not either
// specifically implemented in Checkpoint class or specialized below.
template<typename T>
void
ReadCheckpoint(const Checkpoint& chkp, const Key& fieldname, T& t)
{
  UserReadCheckpoint(chkp, fieldname, t);
}



//
// Initialize
// ======================================================================

// This fails to compile if not provided by user.
template<typename T>
bool
Initialize(Teuchos::ParameterList& plist, T& t,
           const Key& fieldname,
           const std::vector<std::string>& subfieldnames) {
  return UserInitialize(plist, t, fieldname, subfieldnames);
}



//
// Specializations for simple data types
// ======================================================================
template<>
void
WriteVis<double>(const Visualization& vis, const Key& fieldname,
                 const std::vector<std::string>& subfieldnames,
                 const double& t);

template<>
void
WriteCheckpoint<double>(const Checkpoint& chkp, const Key& fieldname,
                        const double& t);

template<>
void
ReadCheckpoint<double>(const Checkpoint& chkp, const Key& fieldname,
                       double& t);

template<>
bool
Initialize<double>(Teuchos::ParameterList& plist, double& t,
                   const Key& fieldname,
                   const std::vector<std::string>& subfieldnames);


template<>
void
WriteVis<int>(const Visualization& vis, const Key& fieldname,
                 const std::vector<std::string>& subfieldnames,
                 const int& t);

template<>
void
WriteCheckpoint<int>(const Checkpoint& chkp, const Key& fieldname,
                        const int& t);

template<>
void
ReadCheckpoint<int>(const Checkpoint& chkp, const Key& fieldname,
                       int& t);

template<>
bool
Initialize<int>(Teuchos::ParameterList& plist, int& t,
                   const Key& fieldname,
                   const std::vector<std::string>& subfieldnames);



//
// Specializations for CompositeVector
// ======================================================================

template<>
void
WriteVis<CompositeVector>(const Visualization& vis, const Key& fieldname,
                          const std::vector<std::string>& subfieldnames,
                          const CompositeVector& vec);

template<>
void
WriteCheckpoint<CompositeVector>(const Checkpoint& chkp, const Key& fieldname,
        const CompositeVector& vec);

template<>
void
ReadCheckpoint<CompositeVector>(const Checkpoint& chkp, const Key& fieldname,
        CompositeVector& vec);

template<>
bool
Initialize<CompositeVector>(Teuchos::ParameterList& plist, CompositeVector& t,
                            const Key& fieldname,
                            const std::vector<std::string>& subfieldnames);




//
// Specializations for Op
// ======================================================================

template<>
inline
void
WriteVis<Operators::Op>(const Visualization& vis, const Key& fieldname,
             const std::vector<std::string>& subfieldnames,
             const Operators::Op& vec) {}

template<>
inline
void
WriteCheckpoint<Operators::Op>(const Checkpoint& chkp, const Key& fieldname,
                    const Operators::Op& vec) {}

template<>
inline
void
ReadCheckpoint<Operators::Op>(const Checkpoint& chkp, const Key& fieldname,
                   Operators::Op& vec){}

template<>
inline
bool
Initialize<Operators::Op>(Teuchos::ParameterList& plist, Operators::Op& t,
               const Key& fieldname,
               const std::vector<std::string>& subfieldnames) { return true; }




template<class Op_t>
class Op_Factory {
 public:
  Op_Factory() {}

  void set_mesh(const Teuchos::RCP<AmanziMesh::Mesh>& mesh) { mesh_ = mesh; }
  void set_name(const std::string& name) { name_ = name; }
  Teuchos::RCP<Op_t> Create() const {
    return Teuchos::rcp(new Op_t(name_, mesh_));
  }

 private:
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  std::string name_;
};



} // namespace
} // namespace
    
#endif
