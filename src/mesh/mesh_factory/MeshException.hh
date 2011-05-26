// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   MeshException.hh
 * @author William A. Perkins
 * @date Mon Mar 14 09:02:33 2011
 * 
 * @brief declaration of several exceptions used by MeshFactory and
 * related code
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 12, 2011 by William A. Perkins
// Last Change: Mon Mar 14 09:02:33 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

// SCCS ID: $Id$ Battelle PNL

#ifndef _MeshException_hh_
#define _MeshException_hh_

#include "errors.hh"
#include "MeshFramework.hh"

namespace Amanzi {
namespace AmanziMesh {

  /// A special exception type for MeshFactory errors
  class Message : public Errors::Message {
  public:
    explicit Message(void) : Errors::Message() {}
    explicit Message(const char* message) : Errors::Message(message) {}
    virtual ~Message(void) throw() {};
  };

  /// A special exception to identify framework problems
  /**
   * This is only thrown if a specific framework is requested but
   * not available
   * 
   */
  class FrameworkMessage : public Message {
  protected:
    Framework my_framework;
    void initialize(void)
    {
      add_data(framework_name(my_framework).c_str());
      add_data(": mesh framework not available");
    }

  public:
    explicit FrameworkMessage(const Framework& fw) : Message(), my_framework(fw) 
    {
      this->initialize();
    }
    explicit FrameworkMessage(const Framework& fw, const char* message) : 
      Message(), my_framework(fw) 
    { 
      this->initialize();
      add_data(": ");
      add_data(message);
    }
    ~FrameworkMessage(void) throw() {};
  };
    
  /// A special exception to identify file problems
  class FileMessage : public Message {
  public:
    explicit FileMessage(void) : Message() {}
    explicit FileMessage(const char* message) : Message(message) {}
    ~FileMessage(void) throw() {};
  };    
    
} // namespace AmanziMesh
} // namespace Amanzi

#endif
