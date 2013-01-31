/*
This is the transport component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef __TRANSPORT_CONSTANTS_HH__
#define __TRANSPORT_CONSTANTS_HH__


namespace Amanzi {
namespace AmanziTransport {

const int TRANSPORT_NULL = 0;
const int TRANSPORT_FLOW_AVAILABLE = 1;
const int TRANSPORT_STATE_BEGIN = 2;
const int TRANSPORT_STATE_COMPLETE = 3;

const int TRANSPORT_INTERNAL_ERROR = 911;  // contact (lipnikov@lanl.gov)

const double TRANSPORT_LARGE_TIME_STEP = 1e+99;
const double TRANSPORT_SMALL_TIME_STEP = 1e-12;

const int TRANSPORT_BC_CONSTANT_TCC = 1;
const int TRANSPORT_BC_DISPERSION_FLUX = 2;
const int TRANSPORT_BC_NULL = 3;

const int TRANSPORT_FLOW_STEADYSTATE = 1;
const int TRANSPORT_FLOW_TRANSIENT = 2;

const double TRANSPORT_CONCENTRATION_OVERSHOOT = 1e-6;
const double TRANSPORT_CONCENTRATION_INFINITY = 1e+99;

const int TRANSPORT_MAX_FACES = 14;  // Kelvin's tetrakaidecahedron
const int TRANSPORT_MAX_NODES = 47;  // These olyhedron parameters must
const int TRANSPORT_MAX_EDGES = 60;  // be calculated in Init().

const int TRANSPORT_DISPERSIVITY_MODEL_NULL = 1;
const int TRANSPORT_DISPERSIVITY_MODEL_ISOTROPIC = 2;
const int TRANSPORT_DISPERSIVITY_MODEL_BEAR = 3;
const int TRANSPORT_DISPERSIVITY_MODEL_LICHTNER = 4;

const int TRANSPORT_LIMITER_BARTH_JESPERSEN = 1; 
const int TRANSPORT_LIMITER_TENSORIAL = 2;
const int TRANSPORT_LIMITER_KUZMIN = 3;
const double TRANSPORT_LIMITER_TOLERANCE = 1e-14;

const int TRANSPORT_VERBOSITY_NONE = 0;
const int TRANSPORT_VERBOSITY_LOW = 1;
const int TRANSPORT_VERBOSITY_MEDIUM = 2;
const int TRANSPORT_VERBOSITY_HIGH = 3;
const int TRANSPORT_VERBOSITY_EXTREME = 4;

const int TRANSPORT_AMANZI_VERSION = 2;  

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

