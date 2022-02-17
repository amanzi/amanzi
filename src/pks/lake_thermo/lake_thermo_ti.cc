/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "boost/math/special_functions/fpclassify.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "Debugger.hh"
#include "BoundaryFunction.hh"
#include "FieldEvaluator.hh"
#include "lake_thermo_pk.hh"
#include "Op.hh"

namespace Amanzi {
namespace LakeThermo {

#define DEBUG_FLAG 1
#define MORE_DEBUG_FLAG 0

// Lake_Thermo_PK is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void Lake_Thermo_PK::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  Teuchos::OSTab tab = vo_->getOSTab();

  bool ice_cover_ = false; // first always assume that there is no ice

  // get temperature
  Teuchos::RCP<const CompositeVector> temp = S_inter_->GetFieldData(temperature_key_);

  Epetra_MultiVector& g_c = *g->Data()->ViewComponent("cell",false);
  unsigned int ncomp = g_c.MyLength();

  std::cout << "ncomp in ti = " << ncomp << std::endl;

	const Epetra_MultiVector& temp_v = *temp->ViewComponent("cell",false);

//	int ncomp = temp->size(*comp, false);

	int i_ice_max = 0;
	int i_ice_min = ncomp;

	for (int i=0; i!=ncomp; ++i) {
	  if (temp_v[0][i] < 273.15) { // check if there is ice cover
//		std::cout << "temp_v[0][" << i << "] = " << temp_v[0][i] << std::endl;
		ice_cover_ = true;
		i_ice_max = i;
	  }
	} // i

	for (int i=ncomp-2; i!=0; --i) {
	  if (temp_v[0][i] < 273.15) { // check if there is ice cover
//		std::cout << "temp_v[0][" << i << "] = " << temp_v[0][i] << std::endl;
		ice_cover_ = true;
		i_ice_min = i;
	  }
	} // i

	std::cout << "ncomp = " << ncomp << std::endl;
	std::cout << "i_ice_max = " << i_ice_max << ", i_ice_min = " << i_ice_min << std::endl;

//	// if melting occured at the top, swap cells
//	for (int i=0; i!=ncomp; ++i) {
//	  if (ice_cover_ && i < i_ice_max && temp_v[0][i] >= 273.15 ) {
//		temp_v[0][i] = temp_v[0][i+1];
//		temp_v[0][i+i_ice_max] = temp_v[0][i+i_ice_max+1];
//	  }
//	} // i

//	// if melting occured at the top, swap cells
//	for (int i=ncomp-1; i!=1; --i) {
//	  if (ice_cover_ && i > i_ice_max && temp_v[0][i] >= 273.15 ) {
//		std::cout << "swapping cells at i = " << i << std::endl;
//		temp_v[0][i] = temp_v[0][i-1];
//		temp_v[0][i_ice_min] = temp_v[0][i_ice_min-1];
//	  }
//	} // i


	int d_thawed = ncomp-1-i_ice_max;
	int d_ice = i_ice_max-i_ice_min+1;

	std::cout << "d_thawed = " << d_thawed << std::endl;
	std::cout << "d_ice = " << d_ice << std::endl;

	std::vector<double> temp_new(ncomp);

	for (int i=0; i!=ncomp; ++i) {
		temp_new[i] = temp_v[0][i];
	}

	if (ice_cover_ && d_thawed > 0) {
		// if thawing occured at the top, swap cells
		for (int i=0; i < std::min(d_thawed,d_ice); ++i) {
			std::cout << "swapping cells at i = " << i << std::endl;
//			std::cout << "temp_v[0][" << i_ice_max+1+i << "] = temp_v[0][" << std::max(i_ice_min+1+i,i_ice_max-d_thawed+1+i) << "]" << std::endl;
//			std::cout << "temp_v[0][" << i_ice_min+i << "] = temp_v[0][" << i_ice_max+1+i << "]" << std::endl;
//			temp_new[i_ice_max+1+i] = temp_v[0][std::max(i_ice_min+1+i,i_ice_max-d_thawed+1+i)];
//			temp_new[i_ice_min+i] = temp_v[0][i_ice_max+1+i];
			std::cout << "temp_v[0][" << ncomp-1-i << "] = temp_v[0][" << i_ice_max-i << "]" << std::endl;
			std::cout << "temp_v[0][" << ncomp-d_ice-1-i << "] = temp_v[0][" << ncomp-1-i << "]" << std::endl;
			temp_new[ncomp-1-i] = temp_v[0][i_ice_max-i];
			temp_new[ncomp-d_ice-1-i] = temp_v[0][ncomp-1-i];
		} // i

		std::cout << "Temperature after swapping cells" << std::endl;
		for (int i=0; i!=ncomp; ++i) {
			temp_v[0][i] = temp_new[i];
			std::cout << "temp_v[0][" << i << "] = " << temp_v[0][i] << std::endl;
		}

	}

  // increment, get timestep
  niter_++;
  double h = t_new - t_old;

  // pointer-copy temperature into states and update any auxilary data
  Solution_to_State(*u_new, S_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << h << std::endl;

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("T_old"); vnames.push_back("T_new");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_inter_->GetFieldData(key_).ptr()); vecs.push_back(u.ptr());
  db_->WriteVectors(vnames, vecs, true);

  // vnames[0] = "sl"; vnames[1] = "si";
  // vecs[0] = S_next_->GetFieldData("saturation_liquid").ptr();
  // vecs[1] = S_next_->GetFieldData("saturation_ice").ptr();
  // db_->WriteVectors(vnames, vecs, false);
#endif

  // update boundary conditions
  bc_temperature_->Compute(t_new);
  bc_diff_flux_->Compute(t_new);
//  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  Teuchos::ParameterList& param_list = plist_->sublist("parameters");
  FunctionFactory fac;
  Teuchos::RCP<Function> r_func_ = Teuchos::rcp(fac.Create(param_list.sublist("precipitation")));
  Teuchos::RCP<Function> E_func_ = Teuchos::rcp(fac.Create(param_list.sublist("evaporation")));

  std::vector<double> args(1);
  args[0] = S_inter_->time();
  r_ = (*r_func_)(args);
  E_ = (*E_func_)(args);

  // Compute evaporartion rate
    Teuchos::ParameterList& param_list_eval = plist_->sublist("surface flux evaluator");
    Teuchos::ParameterList& param_list_param = param_list_eval.sublist("parameters");
    // read parameters from the met data
    Teuchos::RCP<Amanzi::Function> SS_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("solar radiation")));
    Teuchos::RCP<Amanzi::Function> E_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("atmospheric downward radiation")));
    Teuchos::RCP<Amanzi::Function> T_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("air temperature")));
    Teuchos::RCP<Amanzi::Function> H_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("air humidity")));
    Teuchos::RCP<Amanzi::Function> P_a_func_ = Teuchos::rcp(fac.Create(param_list_param.sublist("atmospheric pressure")));

    args[0] = S_inter_->time();
    double SS = (*SS_func_)(args);
    double E_a = (*E_a_func_)(args);
    double T_a = (*T_a_func_)(args);
    double q_a = (*H_a_func_)(args);
    double P_a = (*P_a_func_)(args);

//    // get temperature
//    const Epetra_MultiVector& temp_v = *S_inter_->GetFieldData(temperature_key_)
//		->ViewComponent("cell",false);

    double T_s = temp_v[0][ncomp-1];

	double b1_vap   = 610.78;        // Coefficient [N m^{-2} = kg m^{-1} s^{-2}]
    double b3_vap   = 273.16;        // Triple point [K]
	double b2w_vap  = 17.2693882;    // Coefficient (water)
	double b2i_vap  = 21.8745584;    // Coefficient (ice)
	double b4w_vap  = 35.86;         // Coefficient (temperature) [K]
	double b4i_vap  = 7.66;          // Coefficient (temperature) [K]

	// Saturation water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
	double wvpres_s;

	double h_ice = 0.;
	double h_Ice_min_flk = 1.e-9;
	if (h_ice < h_Ice_min_flk) {  // Water surface
		wvpres_s = b1_vap*std::exp(b2w_vap*(T_s-b3_vap)/(T_s-b4w_vap));
	} else {                      // Ice surface
		wvpres_s = b1_vap*std::exp(b2i_vap*(T_s-b3_vap)/(T_s-b4i_vap));
	}

	// Saturation specific humidity at T=T_s
	double tpsf_R_dryair    = 2.8705e2;  // Gas constant for dry air [J kg^{-1} K^{-1}]
	double tpsf_R_watvap    = 4.6151e2;  // Gas constant for water vapour [J kg^{-1} K^{-1}]
	double tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap;
	double q_s = tpsf_Rd_o_Rv*wvpres_s/(P_a-(1.-tpsf_Rd_o_Rv)*wvpres_s);

	double height_tq = 2;
	double tpsf_kappa_q_a   = 2.4e-05; // Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

	double LE = -tpsf_kappa_q_a*(q_a-q_s)/height_tq;

	double rho_a = P_a/tpsf_R_dryair/T_s/(1.+(1./tpsf_Rd_o_Rv-1.)*q_s);

	double tpsf_c_a_p  = 1.005e3; // Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
	double tpsf_L_evap = 2.501e6; // Specific heat of evaporation [J kg^{-1}]
	double tpl_L_f     = 3.3e5;   // Latent heat of fusion [J kg^{-1}]

	double Q_watvap   = LE*rho_a;
	LE = tpsf_L_evap;
	if (h_ice >= h_Ice_min_flk) LE = LE + tpl_L_f;   // Add latent heat of fusion over ice
	LE = Q_watvap*LE;

	double row0 = 1.e+3;
	double evap_rate = LE/(row0*tpsf_L_evap)*1000.;
	evap_rate = abs(evap_rate);

	std::cout << "evap_rate = " << evap_rate << std::endl;

	E_ = evap_rate;

	std::cout << "precip rate = " << r_ << std::endl;

  // update depth
  double dt = S_next_->time() - S_inter_->time();
  double dhdt = r_ - E_ - R_s_ - R_b_;
  h_ += dhdt*dt;

  std::cout << "h_ = " << h_ << std::endl;

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

//  std::cout << "Starting res" << std::endl;
//  res->Print(std::cout);

  // diffusion term, implicit
  ApplyDiffusion_(S_next_.ptr(), res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("K",S_next_->GetFieldData(conductivity_key_).ptr(),true);
  db_->WriteVector("res (diff)", res.ptr(), true);
#endif

//  std::cout << "After ApplyDiffusion_" << std::endl;
//  res->Print(std::cout);

  // source terms
  AddSources_(S_next_.ptr(), res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("res (src)", res.ptr());
#endif

//  std::cout << "After AddSources_" << std::endl;
//  res->Print(std::cout);

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  vnames[0] = "T_old";
  vnames[1] = "T_new";
  vecs[0] = S_inter_->GetFieldData(temperature_key_).ptr();
  vecs[1] = S_next_->GetFieldData(temperature_key_).ptr();
  db_->WriteVectors(vnames, vecs, true);
  db_->WriteVector("res (acc)", res.ptr());
#endif

//  std::cout << "After AddAccumulation_" << std::endl;
//  res->Print(std::cout);

  // advection term
  if (implicit_advection_) {
    AddAdvection_(S_next_.ptr(), res.ptr(), true);
  } else {
    AddAdvection_(S_inter_.ptr(), res.ptr(), true);
  }
#if DEBUG_FLAG
  db_->WriteVector("res (adv)", res.ptr());
#endif

//  std::cout << "After AddAdvection_" << std::endl;
//  res->Print(std::cout);

  // Dump residual to state for visual debugging.
#if MORE_DEBUG_FLAG
  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << domain_prefix_ << "energy_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << domain_prefix_ << "energy_solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *u;
  }
#endif

  g->Data()->Print(std::cout);

//  exit(0);

//  if (ice_cover_) exit(0);

  // save file with temperature distribution
  std::cout << "niter_ = " << niter_ << std::endl;
  if (int(t_new) % 86400 == 0) {
	  std::ofstream tempfile;
	  tempfile.open ("temperature.txt", std::ios::app);
	  for (int i=0; i!=ncomp; ++i) {
		  const AmanziGeometry::Point& xc = mesh_->cell_centroid(i);
		  tempfile << temp_v[0][i] << " ";
	  }
	  tempfile << "\n";
  }

};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int Lake_Thermo_PK::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
#if DEBUG_FLAG
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon application:" << std::endl;
  db_->WriteVector("T_res", u->Data().ptr(), true);
#endif

//  preconditioner_->PrintDiagnostics();

  // apply the preconditioner
//  std::cout << "Printing data" << std::endl;
//  u->Data()->Print(std::cout);
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());
//  std::cout << "ApplyPreconditioner ierr = " << ierr << std::endl;

#if DEBUG_FLAG
  db_->WriteVector("PC*T_res", Pu->Data().ptr(), true);
#endif
  
  Pu->Data()->ViewComponent("boundary_face")->PutScalar(0.0); // correction 01/22/21

//  std::cout << "Lake ApplyPreconditioner" << std::endl;
//  Pu->Data()->Print(std::cout);

//  std::cout << *Pu->Data()->ViewComponent("cell") << std::endl;

//  std::cout << "ApplyPreconditioner ierr = " << ierr << std::endl;
//  std::cout << *preconditioner_->A() << std::endl;

//  int ierr = 0;
//  *Pu = *u;

  return (ierr > 0) ? 0 : 1;
};

// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void Lake_Thermo_PK::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.

  AMANZI_ASSERT(std::abs(S_next_->time() - t) <= 1.e-4*t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, S_next_);

  Teuchos::RCP<const CompositeVector> temp = S_next_ -> GetFieldData(key_);

  // update boundary conditions
  bc_temperature_->Compute(S_next_->time());
  bc_diff_flux_->Compute(S_next_->time());
//  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_(S_next_.ptr());

  // div K_e grad u
  UpdateConductivityData_(S_next_.ptr());
  if (jacobian_) UpdateConductivityDerivativeData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> conductivity =
      S_next_->GetFieldData(uw_conductivity_key_);

  // jacobian term
  Teuchos::RCP<const CompositeVector> dKdT = Teuchos::null;
  if (jacobian_) {
    if (!duw_conductivity_key_.empty()) {
      dKdT = S_next_->GetFieldData(duw_conductivity_key_);
    } else {
      dKdT = S_next_->GetFieldData(dconductivity_key_);
    }
  }

  // create local matrices
  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(conductivity, dKdT);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
  preconditioner_diff_->ApplyBCs(true, true, true);

  if (jacobian_) {
    Teuchos::RCP<CompositeVector> flux = Teuchos::null;

    flux = S_next_->GetFieldData(energy_flux_key_, name_);
    preconditioner_diff_->UpdateFlux(up->Data().ptr(), flux.ptr());

    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), up->Data().ptr());
  }

  // update with accumulation terms
  // -- update the accumulation derivatives, de/dT
  S_next_->GetFieldEvaluator(energy_key_)
      ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
  const Epetra_MultiVector& de_dT = *S_next_->GetFieldData(Keys::getDerivKey(energy_key_, key_))
      ->ViewComponent("cell",false);
  unsigned int ncells = de_dT.MyLength();

  CompositeVector acc(S_next_->GetFieldData(energy_key_)->Map());
  auto& acc_c = *acc.ViewComponent("cell", false);

#if DEBUG_FLAG
//  db_->WriteVector("    de_dT", S_next_->GetFieldData(Keys::getDerivKey(temperature_key_, key_)).ptr());
#endif

  if (coupled_to_subsurface_via_temp_ || coupled_to_subsurface_via_flux_) {
    // do not add in de/dT if the height is 0
    const Epetra_MultiVector& pres = *S_next_->GetFieldData(Keys::getKey(domain_,"pressure"))
        ->ViewComponent("cell",false);
    const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

    for (unsigned int c=0; c!=ncells; ++c) {
      acc_c[0][c] = pres[0][c] >= patm ? de_dT[0][c] / h : 0.;
    }
  } else {
    if (precon_used_) {
      for (unsigned int c=0; c!=ncells; ++c) {
        //      AMANZI_ASSERT(de_dT[0][c] > 1.e-10);
        // ?? Not using e_bar anymore apparently, though I didn't think we were ever.  Need a nonzero here to ensure not singlar.
        acc_c[0][c] = std::max(de_dT[0][c], 1.e-12) / h;
      }
    } else {
      for (unsigned int c=0; c!=ncells; ++c) {
        //      AMANZI_ASSERT(de_dT[0][c] > 1.e-10);
        // ?? Not using e_bar anymore apparently, though I didn't think we were ever.  Need a nonzero here to ensure not singlar.
        // apply a diagonal shift manually for coupled problems
        acc_c[0][c] = de_dT[0][c] / h + 1.e-6;
      }
    }
  }
  preconditioner_acc_->AddAccumulationTerm(acc, "cell");

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(S_next_.ptr(), h);

//  // update with advection terms
//  if (implicit_advection_ && implicit_advection_in_pc_) {
//    Teuchos::RCP<const CompositeVector> mass_flux = S_next_->GetFieldData(flux_key_);
//    S_next_->GetFieldEvaluator(enthalpy_key_)
//        ->HasFieldDerivativeChanged(S_next_.ptr(), name_, key_);
////    Teuchos::RCP<const CompositeVector> dhdT = S_next_->GetFieldData(Keys::getDerivKey(enthalpy_key_, key_));
//    Teuchos::RCP<const CompositeVector> dhdT = S_next_->GetFieldData(Keys::getDerivKey(temperature_key_, key_));
//    preconditioner_adv_->Setup(*mass_flux);
//    preconditioner_adv_->SetBCs(bc_adv_, bc_adv_);
//    preconditioner_adv_->UpdateMatrices(mass_flux.ptr(), dhdT.ptr());
//    ApplyDirichletBCsToEnthalpy_(S_next_.ptr());  //!!!!!!!!!!!!! IMPORTANT !!!!!!
//    preconditioner_adv_->ApplyBCs(false, true, false);
//
//  }

  // Apply boundary conditions.
  preconditioner_diff_->ApplyBCs(true, true, true);
};

//double Lake_Thermo_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
//        Teuchos::RCP<const TreeVector> res) {
//
//  const Epetra_MultiVector& temp = *S_inter_->GetFieldData(temperature_key_)
//        ->ViewComponent("cell",true);
//  const Epetra_MultiVector& cv = *S_inter_->GetFieldData(cell_vol_key_)
//        ->ViewComponent("cell",true);
//
//  int ncells_owned = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
//
//  double err_max = 0.;
//  double err_L1 = 0.;
//
//  for (int c = 0; c < ncells_owned; c++) {
//    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
//    double T0 = 280., T1 = 400.;
//    double coef = log(T1/T0);
//    double temp_ex_c = T0*exp(coef*xc[2]);
//    err_max = std::max(err_max,std::abs(temp_ex_c-temp[0][c]));
//    err_L1 += std::abs(temp_ex_c-temp[0][c])*cv[0][c];
//  }
//
//  double err_max_tmp(err_max);
//  double err_L1_tmp(err_L1);
//
//  mesh_->get_comm()->MaxAll(&err_max_tmp, &err_max, 1);
//  mesh_->get_comm()->SumAll(&err_L1_tmp, &err_L1, 1);
//
////  if (mesh_->get_comm()->MyPID() == 0) {
////    std::cout << "err_max = " << err_max << std::endl;
////    std::cout << "err_L1  = " << err_L1 << std::endl;
////  }
//
//  return PK_PhysicalBDF_Default::ErrorNorm(u,res);
//
//}

double Lake_Thermo_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
                                    Teuchos::RCP<const TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // Relative error in cell-centered temperature
  const Epetra_MultiVector& uc = *u->Data()->ViewComponent("cell", false);
  const Epetra_MultiVector& duc = *du->Data()->ViewComponent("cell", false);

  int ncells_owned = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  double error_t(0.0);
  double ref_temp(273.0);
  for (int c = 0; c < ncells_owned; c++) {
    double tmp = fabs(duc[0][c]) / (fabs(uc[0][c] - ref_temp) + ref_temp);
    if (tmp > error_t) {
      error_t = tmp;
    }
  }

  // Cell error is based upon error in energy conservation relative to
  // a characteristic energy
  double error_e(0.0);

  double error = std::max(error_t, error_e);

#ifdef HAVE_MPI
  double buf = error;
  du->Data()->Comm()->MaxAll(&buf, &error, 1);  // find the global maximum
#endif

  return error;
}


/*
// -----------------------------------------------------------------------------
// Default enorm that uses an abs and rel tolerance to monitor convergence.
// -----------------------------------------------------------------------------
double Lake_Thermo_PK::ErrorNorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> res) {
  // Abs tol based on old conserved quantity -- we know these have been vetted
  // at some level whereas the new quantity is some iterate, and may be
  // anything from negative to overflow.

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  const Epetra_MultiVector& temp_c = *S_inter_->GetFieldData(temperature_key_)->ViewComponent("cell");

  std::cout << "PRINTING Temperature: " << std::endl;
  for (int c = 0; c < ncells_owned; c++) {
    std::cout << temp_c[0][c] << std::endl;
  }
  std::cout << "end" << std::endl;



  int cycle = S_next_->cycle();
  if (cycle == 32) {
    std::cout << "we is here" << std::endl;
  }

  S_inter_->GetFieldEvaluator(energy_key_)->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& energy = *S_inter_->GetFieldData(energy_key_)
      ->ViewComponent("cell",true);

//  S_inter_->GetFieldEvaluator(wc_key_)->HasFieldChanged(S_inter_.ptr(), name_);
//  const Epetra_MultiVector& wc = *S_inter_->GetFieldData(wc_key_)
//      ->ViewComponent("cell",true);

  const Epetra_MultiVector& wc = *S_inter_->GetFieldData(wc_key_)
          ->ViewComponent("cell",true);

  const Epetra_MultiVector& cv = *S_inter_->GetFieldData(cell_vol_key_)
      ->ViewComponent("cell",true);

//  std::cout << "mass_atol_ = " << mass_atol_ << std::endl;
//  std::cout << "soil_atol_ = " << soil_atol_ << std::endl;
//  std::cout << "atol_ = " << atol_ << std::endl;




  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "ENorm (Infnorm) of: " << conserved_key_ << ": " << std::endl;

  Teuchos::RCP<const CompositeVector> dvec = res->Data();
  double h = S_next_->time() - S_inter_->time();

  Teuchos::RCP<const Comm_type> comm_p = mesh_->get_comm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp=dvec->begin();
       comp!=dvec->end(); ++comp) {
    double enorm_comp = 0.0;
    int enorm_loc = -1;
    const Epetra_MultiVector& dvec_v = *dvec->ViewComponent(*comp, false);

    if (*comp == std::string("cell")) {
      // error done in two parts, relative to mass but absolute in
      // energy since it doesn't make much sense to be relative to
      // energy
      int ncells = dvec->size(*comp,false);
      for (unsigned int c=0; c!=ncells; ++c) {
        double mass = std::max(mass_atol_, wc[0][c] / cv[0][c]);
        double energy = mass * atol_ + soil_atol_;
        double enorm_c = std::abs(h * dvec_v[0][c]) / (energy*cv[0][c]);

//        std::cout << "c = " << c << ", enorm_c = " << enorm_c << std::endl;
//        std::cout << "mass = " << mass << std::endl;
//        std::cout << "energy = " << energy << std::endl;
//        std::cout << "h = " << h << std::endl;
//        std::cout << "dvec_v[0][c] = " << dvec_v[0][c] << std::endl;
//        std::cout << "cv[0][c] = " << cv[0][c] << std::endl;
//        std::cout << "wc[0][c] = " << wc[0][c] << std::endl;

        if (enorm_c > enorm_comp) {
          enorm_comp = enorm_c;
          enorm_loc = c;
        }
      }

    } else if (*comp == std::string("face")) {
      // error in flux -- relative to cell's extensive conserved quantity
      int nfaces = dvec->size(*comp, false);

      for (unsigned int f=0; f!=nfaces; ++f) {
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);
        double cv_min = cells.size() == 1 ? cv[0][cells[0]]
            : std::min(cv[0][cells[0]],cv[0][cells[1]]);
        double mass_min = cells.size() == 1 ? wc[0][cells[0]]/cv[0][cells[0]]
          : std::min(wc[0][cells[0]]/cv[0][cells[0]], wc[0][cells[1]]/cv[0][cells[1]]);
        mass_min = std::max(mass_min, mass_atol_);

        double energy = mass_min * atol_ + soil_atol_;
        double enorm_f = fluxtol_ * h * std::abs(dvec_v[0][f])
          / (energy * cv_min);

//        std::cout << "f = " << f << ", enorm_f = " << enorm_f << std::endl;

        if (enorm_f > enorm_comp) {
          enorm_comp = enorm_f;
          enorm_loc = f;
        }
      }

    } else {
      // boundary face components had better be effectively identically 0
      double norm;
      dvec_v.Norm2(&norm);
      AMANZI_ASSERT(norm < 1.e-15);
    }


    // Write out Inf norms too.
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      double infnorm(0.);
      dvec_v.NormInf(&infnorm);

      ENorm_t err;
      ENorm_t l_err;
      l_err.value = enorm_comp;
      l_err.gid = dvec_v.Map().GID(enorm_loc);

      int ierr;
      ierr = MPI_Allreduce(&l_err, &err, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
      AMANZI_ASSERT(!ierr);
      *vo_->os() << "  ENorm (" << *comp << ") = " << err.value << "[" << err.gid << "] (" << infnorm << ")" << std::endl;
    }

    enorm_val = std::max(enorm_val, enorm_comp);
  }

  double enorm_val_l = enorm_val;

  int ierr;
  ierr = MPI_Allreduce(&enorm_val_l, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, comm);
  AMANZI_ASSERT(!ierr);
  return enorm_val;
};
*/


} // namespace LakeThermo
} // namespace Amanzi
