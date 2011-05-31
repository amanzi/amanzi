#include "Transient_Richards_PK.hpp"

#include "RichardsProblem.hpp"


Transient_Richards_PK::Transient_Richards_PK(Teuchos::ParameterList &plist, const Teuchos::RCP<const Flow_State> FS_) : FS(FS_)
{
  // Create the flow boundary conditions object.
  Teuchos::ParameterList bc_plist = plist.sublist("Flow BC");
  bc = Teuchos::rcp<FlowBC>(new FlowBC(bc_plist, FS->mesh()));

  // Create the Richards flow problem.
  problem = new RichardsProblem(FS->mesh(), plist.sublist("Richards Problem"), bc);

  // Create the solution vectors.
  solution = new Epetra_Vector(problem->Map());
  pressure = problem->CreateCellView(*solution);
  richards_flux = new Epetra_Vector(problem->FaceMap());

  

};

Transient_Richards_PK::~Transient_Richards_PK()
{
  delete richards_flux;
  delete pressure;
  delete solution;
  delete problem;
};


int Transient_Richards_PK::advance()
{
  // Set problem parameters.
  problem->SetFluidDensity(FS->fluid_density());
  problem->SetFluidViscosity(FS->fluid_viscosity());
  problem->SetPermeability(FS->permeability());
  problem->SetGravity(FS->gravity());
  problem->SetFlowState(FS);


  // Derive the Richards fluxes on faces
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *richards_flux, l1_error);
  std::cout << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;

}

void Transient_Richards_PK::GetSaturation(Epetra_Vector &s) const
{
  //for (int i = 0; i < s.MyLength(); ++i) s[i] = 1.0;

  problem->DeriveVanGenuchtenSaturation(*pressure, s);

}
