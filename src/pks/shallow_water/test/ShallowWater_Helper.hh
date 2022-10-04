#include "CompositeVector.hh"
#include "OutputXDMF.hh"
#include "State.hh"

void IO_Fields(double t_out, int iter, int MyPID,
               Amanzi::OutputXDMF& io, const Amanzi::State& S,
               Epetra_MultiVector* hh_ex,
               Epetra_MultiVector* vel_ex)
{
  const auto& hh = *S.Get<Amanzi::CompositeVector>("surface-ponded_depth").ViewComponent("cell");
  const auto& ht = *S.Get<Amanzi::CompositeVector>("surface-total_depth").ViewComponent("cell");
  const auto& vel = *S.Get<Amanzi::CompositeVector>("surface-velocity").ViewComponent("cell");
  const auto& q = *S.Get<Amanzi::CompositeVector>("surface-discharge").ViewComponent("cell");
  const auto& B = *S.Get<Amanzi::CompositeVector>("surface-bathymetry").ViewComponent("cell");

  // create pid vector
  Epetra_MultiVector pid(hh);
  for (int c = 0; c < pid.MyLength(); c++) pid[0][c] = MyPID;

  auto comp = Amanzi::AmanziMesh::CELL;
  io.InitializeCycle(t_out, iter, "");
  io.WriteVector(*hh(0), "depth", comp);
  io.WriteVector(*ht(0), "total_depth", comp);
  io.WriteVector(*vel(0), "vx", comp);
  io.WriteVector(*vel(1), "vy", comp);
  io.WriteVector(*q(0), "qx", comp);
  io.WriteVector(*q(1), "qy", comp);
  io.WriteVector(*B(0), "B", comp);
  io.WriteVector(*pid(0), "pid", comp);

  if (hh_ex) {
    io.WriteVector(*(*hh_ex)(0), "hh_ex", comp);
    io.WriteVector(*(*vel_ex)(0), "vx_ex", comp);
    io.WriteVector(*(*vel_ex)(1), "vy_ex", comp);
  }
  io.FinalizeCycle();
}
