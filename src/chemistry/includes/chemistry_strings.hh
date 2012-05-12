/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_STRINGS_HH_
#define AMANZI_CHEMISTRY_STRINGS_HH_

#include <string>

namespace amanzi {
namespace chemistry {
namespace strings {

/*
static const std::string ("");
*/

/*
**  Verbosity Levels
*/
static const std::string kVerbositySilent("silent");
static const std::string kVerbosityTerse("terse");
static const std::string kVerbosityVerbose("verbose");
static const std::string kVerbosityWarning("warning");
static const std::string kVerbosityError("error");

static const std::string kVerbosityDebugBeaker("Debug Input File");
static const std::string kVerbosityDebugMineralKinetics("Debug Mineral Kinetics");
static const std::string kVerbosityDebugInputFile("debug Input File");
static const std::string kVerbosityDebugDatabase("debug Database");
static const std::string kVerbosityDebugActivityModel("debug ActivityModel");
static const std::string kVerbosityDebugSpeciation("debug Speciation");
static const std::string kVerbosityDebugLinearSolver("debug LinearSolver");
static const std::string kVerbosityDebugChemistryProcessKernel("debug ChemistryProcessKernel");



/*
**  Activity Coefficient strings
*/
static const std::string kActivityModel("Activity Model");
static const std::string kDebyeHuckel("Debye-Huckel");
static const std::string kDebyeHuckelBdot("Debye-Huckel B-dot");
static const std::string kPitzerHWM("pitzer-hwm");
static const std::string kUnit("unit");

/*
**  Species Names
*/
static const std::string kSpeciesWater("H2O");

/*
**  Kinetic Rates
*/
static const std::string kTST("TST");


/*
**  Evaluator Strings
*/
static const std::string kCoordinator("Coordinator");
static const std::string kSpeciation("speciation");
static const std::string kOperatorSplittingNR("operator splitting newton-raphson");
static const std::string kOperatorSplittingODE("operator splitting ode");
static const std::string kGlobalImplicit("global implicit");

/*
**  Database Strings
*/
static const std::string kCrunchFlow("CrunchFlow");
static const std::string kPFloTran("PFloTran");
static const std::string kPFLOTRAN_preprocessed("PFLOTRAN_preprocessed");

/*
**  InputLoader Strings
*/
static const std::string kLoaderChemText("Chem Text");
static const std::string kSuffixChemText("bgd");

/*
**  Interpolation Strings
*/
const std::string kLinearInterpolation = "Linear Interpolation";

/*
**  Linear Solver Strings
*/
const std::string kLU = "LU";

/*
**  Nonlinear Solver Strings
*/
const std::string kNewtonRaphson = "Newton-Raphson";
const std::string kNewtonRaphsonUnderRelaxation = "Newton-Raphson-Under-Relaxation";

/*
**  ODE Solver Strings
*/
const std::string kBackwardEuler = "Backward Euler";
const std::string kWeightedEuler = "Weighted Euler";


}  // namespace strings
}  // namespace chemistry
}  // namespace amanzi
#endif     /* AMANZI_CHEMISTRY_STRINGS_HH_ */
