/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins, others
*/

//! Exceptions for MeshFactory
#ifndef AMANZI_MESH_EXCEPTION_HH_
#define AMANZI_MESH_EXCEPTION_HH_

#include "errors.hh"

namespace Amanzi {
namespace AmanziMesh {

/// A special exception type for MeshFactory errors
class Message : public Errors::Message {
 public:
  explicit Message(void) : Errors::Message(){};
  explicit Message(const char* message) : Errors::Message(message){};
  virtual ~Message(void) noexcept {};
};


// A special exception to identify file problems
class FileMessage : public Message {
 public:
  explicit FileMessage(void) : Message(){};
  explicit FileMessage(const char* message) : Message(message){};
  virtual ~FileMessage(void) noexcept {};
};

// A special exception to identify frameworks problems
class FrameworkMessage : public Message {
 public:
  explicit FrameworkMessage(void) : Message(){};
  explicit FrameworkMessage(const char* message) : Message(message){};
  virtual ~FrameworkMessage(void) noexcept {};
};

} // namespace AmanziMesh
} // namespace Amanzi

#endif
