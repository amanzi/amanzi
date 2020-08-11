/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Source term evaluator for the lake temperature model.

  Authors: Svetlana Tokareva (tokareva@lanl.gov)
*/

#include "lake_source_terms_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// constructor format for all derived classes
  LakeThermoSourceEvaluator::LakeThermoSourceEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  InitializeFromPlist_();
}

  LakeThermoSourceEvaluator::LakeThermoSourceEvaluator(const LakeThermoSourceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    internal_enthalpy_key_(other.internal_enthalpy_key_),
    external_enthalpy_key_(other.external_enthalpy_key_),
    mass_source_key_(other.mass_source_key_),
    internal_density_key_(other.internal_density_key_),
    external_density_key_(other.external_density_key_),
    cell_vol_key_(other.cell_vol_key_),
    conducted_source_key_(other.conducted_source_key_),
    include_conduction_(other.include_conduction_),
    source_units_(other.source_units_)
{}

Teuchos::RCP<FieldEvaluator>
LakeThermoSourceEvaluator::Clone() const {
  return Teuchos::rcp(new LakeThermoSourceEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
LakeThermoSourceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  const Epetra_MultiVector& temp = *S->GetFieldData(temperature_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cell_vol_key_)
      ->ViewComponent("cell",false);

  Epetra_MultiVector& res = *result->ViewComponent("cell",false);
  
  //S->GetFieldEvaluator(density_key_)->HasFieldChanged(S.ptr(), name_);

  // evaluate density
  const Epetra_MultiVector& rho =
  *S->GetFieldData(density_key_)->ViewComponent("cell",false);

  // precipitation rate
  Epetra_Vector& rvec = *S->GetConstantVectorData("precipitation", "state");
  double r_ = std::fabs(rvec[1]);

  // precipitation rate
  Epetra_Vector& Evec = *S->GetConstantVectorData("evaporation", "state");
  double E_ = std::fabs(Evec[1]);

  // surface runoff
  Epetra_Vector& Rsvec = *S->GetConstantVectorData("surface runoff", "state");
  double R_s_ = std::fabs(Rsvec[1]);

  // bottom runoff
  Epetra_Vector& Rbvec = *S->GetConstantVectorData("bottom runoff", "state");
  double R_b_ = std::fabs(Rbvec[1]);

  // extinction coefficient
  Epetra_Vector& alpha_e_vec = *S->GetConstantVectorData("extinction coefficient", "state");
  double alpha_e_ = std::fabs(alpha_e_vec[1]);

  // heat capacity of water
  Epetra_Vector& cp_vec = *S->GetConstantVectorData("heat capacity", "state");
  double cp_ = std::fabs(cp_vec[1]);
  cp_ = 4184.;

  double dhdt = r_ - E_ - R_s_ - R_b_;
  double S0 = 1.;

  // THIS WON"T WORK BECAUSE EVALUATOR DOESN"T KNOW THE MESH

  unsigned int ncells = res.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    const AmanziGeometry::Point& xc; // = mesh_->cell_centroid(c);
    res[0][c] = S0*exp(-alpha_e_*h_*xc[0])*(-alpha_e_*h_) + cp_*rho[0][c]*temp[0][c]*dhdt/h_;
  }
}

void
LakeThermoSourceEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
  if (include_conduction_ && wrt_key == conducted_source_key_) {
    *result->ViewComponent("cell",false) = *S->GetFieldData(cell_vol_key_)
        ->ViewComponent("cell",false);
  } else {
    result->PutScalar(0.);
  }
}

void
LakeThermoSourceEvaluator::InitializeFromPlist_() {

  if (my_key_.empty()) {
    if (include_conduction_) {
      my_key_ = plist_.get<std::string>("energy source key",
              "total_energy_source");
    } else {
      my_key_ = plist_.get<std::string>("energy source key",
              "advected_energy_source");
    }
  }
  std::string domain = Keys::getDomain(my_key_);

  internal_enthalpy_key_ = Keys::readKey(plist_, domain, "internal enthalpy", "enthalpy");
  external_enthalpy_key_ = Keys::readKey(plist_, domain, "external enthalpy", "mass_source_enthalpy");
  mass_source_key_ = Keys::readKey(plist_, domain, "mass source", "mass_source");

  dependencies_.insert(internal_enthalpy_key_);
  dependencies_.insert(external_enthalpy_key_);
  dependencies_.insert(mass_source_key_);

  // this handles both surface fluxes (in m/s) and subsurface fluxes (in mol/s)
  std::string source_units = plist_.get<std::string>("mass source units");
  if (source_units == "mol s^-1") {
    source_units_ = SOURCE_UNITS_MOLS_PER_SECOND;
  } else if (source_units == "m s^-1") {
    source_units_ = SOURCE_UNITS_METERS_PER_SECOND;
  } else if (source_units == "mol m^-2 s^-1" ||
             source_units == "mol m^-3 s^-1") {
    source_units_ = SOURCE_UNITS_MOLS_PER_SECOND_PER_METERSD;
  } else {
    Errors::Message message;
    message << "LakeThermoSourceEvaluator: "
            << my_key_
            << ": invalid units \""
            << source_units
            << "\" for \"mass source units\", valid are \"mol s^-1\", \"m s^-1\", \"mol m^-2 s^-1\","
            << " and \"mol m^-3 s^-1\".";
    Exceptions::amanzi_throw(message);
  }

  if (source_units_ == SOURCE_UNITS_METERS_PER_SECOND) {
    internal_density_key_ = plist_.get<std::string>("internal density key",
            Keys::getKey(domain, "molar_density_liquid"));
    external_density_key_ = plist_.get<std::string>("external density key",
            Keys::getKey(domain, "source_molar_density"));

    dependencies_.insert(internal_density_key_);
    dependencies_.insert(external_density_key_);
  }

  // this enables the addition of provided diffusive fluxes as well
  include_conduction_ = plist_.get<bool>("include conduction");
  if (include_conduction_) {
    conducted_source_key_ = plist_.get<std::string>("conducted energy source key",
            Keys::getKey(domain, "conducted_energy_source"));
    dependencies_.insert(conducted_source_key_);
  }

  cell_vol_key_ = plist_.get<std::string>("cell volume key",
          Keys::getKey(domain, "cell_volume"));

}


} //namespace
} //namespace

