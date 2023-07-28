#include "libmesh/face_quad4.h" // for meshing
#include "SystemBase.h" // for access solution
#include "SimpleHeatStructure2.h"

/**
 * This SimpleHeatStructure2 (option 2) follows the traiditional RELAP/TRACE way of finite volume method
 * in structured mesh.
 * As depicted, unlike the true finite volume method, for the left/right surfaces, the surface cell
 * is half the size of a regular internal cells.
 * The positive side of this arrangement is that the surface temperature can be used directly
 * for handling of convective heat transfer and thermal radiation.
 * The drawback is when a Dirichlet boundary condition is applied on the left/right surface temperature,
 * the numerical response becomes unphysical:
 * 1) when the Dirichlet boundary temperature changes with time, the surface cell temperature (a volume average)
 *    responses infinitely fast
 * 2) there is no easy while still reliable way to estimate the heat flux on the left/right surface while still
 *    maintaining local (and therefore global) heat balance.
 *
 *                   -----------------------------------------
 *                   |    |         |         |         |    |
 *                   x    |    x    |    x    |    x    |    x
 *                   |    |         |         |         |    |
 *                   -----------------------------------------
 *                   |    |         |         |         |    |
 *                   x    |    x    |    x    |    x    |    x
 *                   |    |         |         |         |    |
 *                   -----------------------------------------
 *                   |    |         |         |         |    |
 *                   x    |    x    |    x    |    x    |    x
 *                   |    |         |         |         |    |
 *                   -----------------------------------------
 *
 * Another drawback of this approach is when modeling a multi-layer structure, as dipicted,
 * the interface cell (represented by #) is a composite cell with both materials 1 and 2.
 *
 *                   <------------- material 1 -------------> <-------- material 2 -----
 *
 *                   ---------------------------------------- --------------------------
 *                   |    |         |         |         |    |       |              |
 *                   x    |    x    |    x    |    x    | 1  #   2   |      x       |
 *                   |    |         |         |         |    |       |              |
 *                   ---------------------------------------- --------------------------
 *                   |    |         |         |         |    |       |              |
 *                   x    |    x    |    x    |    x    |    #       |      x       |
 *                   |    |         |         |         |    |       |              |
 *                   ---------------------------------------- --------------------------
 *                   |    |         |         |         |    |       |              |
 *                   x    |    x    |    x    |    x    |    #       |      x       |
 *                   |    |         |         |         |    |       |              |
 *                   ---------------------------------------- --------------------------
 */

registerMooseObject("delphiApp", SimpleHeatStructure2);

InputParameters
SimpleHeatStructure2::validParams()
{
  InputParameters params = SimpleHeatStructure::validParams();
  return params;
}

SimpleHeatStructure2::SimpleHeatStructure2(const InputParameters & parameters)
  : SimpleHeatStructure(parameters)
{
}

void
SimpleHeatStructure2::buildMesh()
{
  const std::vector<Real> & pos = getParam<std::vector<Real>>("position");
  const std::vector<Real> & dir = getParam<std::vector<Real>>("orientation");

  Point position = Point(pos[0], pos[1], pos[2]);
  RealVectorValue dir_L = VectorValue<Real>(dir[0], dir[1], dir[2]);
  if (dir_L.norm() < 1e-16)
    mooseError("'orientation' cannot be a zero vector.");
  dir_L = dir_L.unit();

  RealVectorValue x_axis = VectorValue<Real>(1.0, 0.0, 0.0);
  RealVectorValue y_axis = VectorValue<Real>(0.0, 1.0, 0.0);
  RealVectorValue dir_W = dir_L.cross(y_axis);
  if (dir_W.norm() < 1e-16)
    dir_W = dir_L.cross(x_axis);
  dir_W = dir_W.unit();

  _subdomain_id = getNextSubdomainId();
  _subdomain_name = Moose::stringify(_subdomain_id);
  _mesh.setSubdomainName(_subdomain_id, _subdomain_name);

  // points/nodes
  _nodes.resize(_NL + 1);
  for (unsigned j = 0; j < _NL + 1; j++) // loop on row in axial direction
  {
    Point p = position + _dL * j * dir_L; // left most point
    for (unsigned i = 0; i < _NW + 2; i++) // loop on column in width/radial direction
    {
      Real W_step = 0;
      if (i == 0) W_step = 0.0; // don't move
      else if (i == 1 || i == _NW + 1) W_step = 0.5 * _dW; // first/last cell move have dW
      else W_step = _dW; // interior cells move full dW

      p += W_step * dir_W;
      Node * nd = _mesh.getMesh().add_point(p);
      _nodes[j].push_back(nd);
    }
  }

  // elems
  _elems.resize(_NL);
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
    {
      Elem * elem = _mesh.getMesh().add_elem(new Quad4);
      elem->subdomain_id() = _subdomain_id;
      elem->set_node(0) = _nodes[j][i];
      elem->set_node(1) = _nodes[j][i+1];
      elem->set_node(2) = _nodes[j+1][i+1];
      elem->set_node(3) = _nodes[j+1][i];
      _elems[j].push_back(elem);
    }
}

void
SimpleHeatStructure2::setupIC(double * u)
{
  Real Ts_init = getParam<Real>("Ts_init");
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
    {
      _Ts[j][i] = Ts_init;
      u[j * (_NW + 1) + i] = Ts_init;
    }

  _Ts_old = _Ts;
  _Ts_oo = _Ts;
}

void
SimpleHeatStructure2::addPhysicalModel()
{
  // _Ts is arranged as _Ts[Y/axial][X/radial], so when using one-dimensional assumption, the matrix is denser
  _Ts.resize(_NL);
  _Ts_DOFs.resize(_NL);
  for (unsigned j = 0; j < _Ts.size(); j++)
  {
    _Ts[j].resize(_NW + 1, 0.0);

    for (unsigned i = 0; i < _NW + 1; i++)
      _Ts_DOFs[j].push_back(j * (_NW + 1) + i + _DOF_offset);
  }

  _Ts_old = _Ts;
  _Ts_oo = _Ts;

  _volume = _Ts;
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
      _volume[j][i] = (i == 0 || i == _NW) ? 0.5 * _dW * _dL : _dW * _dL;

  _n_DOFs = (_NW + 1) * _NL;

  _solid = &(_sim.getUserObject<ThermalSolidProperties>(getParam<UserObjectName>("solid")));

  // coupling
  if (isParamValid("name_comp_left"))
  {
    std::string comp_name = getParam<std::string>("name_comp_left");
    _pipe_left = dynamic_cast<OneDFlowChannel *>(_sim.getComponentByName(comp_name));
    if (_pipe_left)
    {
      Real aw_left = getParam<Real>("HT_surface_area_density_left");
      _hw_left = getParam<Real>("Hw_left");
      _depth = aw_left * _pipe_left->getArea();
      _pipe_left->addWallHeating(this, _hw_left, aw_left, 0);
    }
    else
      mooseError(comp_name + " is not a OneDFlowChannel.");
  }
  if (isParamValid("name_comp_right"))
  {
    std::string comp_name = getParam<std::string>("name_comp_right");
    _pipe_right = dynamic_cast<OneDFlowChannel *>(_sim.getComponentByName(comp_name));
    if (_pipe_right)
    {
      Real aw_right = getParam<Real>("HT_surface_area_density_right");
      _hw_right = getParam<Real>("Hw_right");
      _depth = aw_right * _pipe_right->getArea(); // TODO: check consistence between this depth and the other one
      _pipe_right->addWallHeating(this, _hw_right, aw_right, 1);
    }
    else
      mooseError(comp_name + " is not a OneDFlowChannel.");
  }
}

void
SimpleHeatStructure2::addExternalVariables()
{
  _sim.addMooseAuxVar("T", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
}

void
SimpleHeatStructure2::updateSolution(double * u)
{
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
      _Ts[j][i] = u[j * (_NW + 1) + i];
}

void
SimpleHeatStructure2::computeTranRes(double * res)
{
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
    {
      Real rho = _solid->rho_from_T(_Ts[j][i]);
      Real cp = _solid->cp_from_T(_Ts[j][i]);

      Real dT_dt = 0.0;
      if ((_sim.TS() == Moose::TI_IMPLICIT_EULER) || (_sim.timeStep() == 1)) // bdf1 or step 1 of bdf2
        dT_dt = (_Ts[j][i] - _Ts_old[j][i]) / _sim.dt();
      else //bdf2 and step > 1
        dT_dt = DELPHI::BDF2_ddt(_Ts[j][i], _Ts_old[j][i], _Ts_oo[j][i], _sim.dt(), _sim.dtOld());

      res[j * (_NW + 1) + i] = _volume[j][i] * rho * cp * dT_dt;
    }

  // (optional, if applicable) west and east boundary DirichletBC, zero transient residual
  if (isParamValid("T_bc_left"))
    for (unsigned j = 0; j < _NL; j++)
      res[j * (_NW + 1) + 0] = 0;

  if (isParamValid("T_bc_right"))
    for (unsigned j = 0; j < _NL; j++)
      res[j * (_NW + 1) + _NW] = 0;
}

void
SimpleHeatStructure2::computeSpatialRes(double * res)
{
  _Q_out_left = 0.0;
  _Q_out_right = 0.0;
  // loop on vertical interior faces
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW; i++)
    {
      Real TW = _Ts[j][i];
      Real TE = _Ts[j][i+1];
      Real T_avg = 0.5 * (TW + TE);
      Real k = _solid->k_from_T(T_avg);
      Real q_flux = -k * (TE - TW) / _dW;

      // West
      res[j * (_NW + 1) + i] += q_flux * _dL;
      // East
      res[j * (_NW + 1) + i + 1] -= q_flux * _dL;

      // use these values to approximate q_wall (left/right) when Dirichlet BC applied on (left/right)
      // this may cause artificial energy imbalance in HS, but there is no other better way to approximate wall heat flux
      if (i == 0) _Q_out_left += q_flux * _dL * _depth;
      if (i == _NW-1) _Q_out_right += (-q_flux * _dL * _depth);
    }

  // heat
  _Q_internal = 0.0;
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
    {
      Real x = _dW * i;
      Real q_vol = M_PI * M_PI / (_width * _width) * sin(M_PI * x);
      _Q_internal += q_vol * _volume[j][i] * _depth;
      res[j * (_NW + 1) + i] -= q_vol * _volume[j][i];
    }

  // convective heat transfer
  if (_pipe_left)
  {
    _Q_out_left = 0.0;
    std::vector<OneDCell*> & fluid_cells = _pipe_left->getCells();
    for (unsigned j = 0; j < _NL; j++)
    {
      Real T_fluid = fluid_cells[j]->T();
      Real q_flux_out = _hw_left * (_Ts[j][0] - T_fluid) * _dL;
      _Q_out_left -= q_flux_out * _depth;
      res[j * (_NW + 1) + 0] += q_flux_out;
    }
  }
  if (_pipe_right)
  {
    _Q_out_right = 0.0;
    std::vector<OneDCell*> & fluid_cells = _pipe_right->getCells();
    for (unsigned j = 0; j < _NL; j++)
    {
      Real T_fluid = fluid_cells[j]->T();
      Real q_flux_out = _hw_right * (_Ts[j][_NW] - T_fluid) * _dL;
      _Q_out_right += q_flux_out * _depth;
      res[j * (_NW + 1) + _NW] += q_flux_out;
    }
  }

  // (optional, if applicable) apply west and east boundary DirichletBC
  if (isParamValid("T_bc_left"))
    for (unsigned j = 0; j < _NL; j++)
      res[j * (_NW + 1) + 0] = _Ts[j][0] - _T_bc_left;

  if (isParamValid("T_bc_right"))
    for (unsigned j = 0; j < _NL; j++)
      res[j * (_NW + 1) + _NW] = _Ts[j][_NW] - _T_bc_right;
}

void
SimpleHeatStructure2::onTimestepBegin()
{
}

void
SimpleHeatStructure2::onTimestepEnd()
{
  _Ts_oo = _Ts_old;
  _Ts_old = _Ts;

  // output
  MooseVariableFieldBase & T_var = _sim.getVariable(0, "T");
  NumericVector<Number> & T_sln = T_var.sys().solution();
  for (unsigned j = 0; j < _elems.size(); j++)
    for (unsigned i = 0; i < _elems[j].size(); i++)
    {
      dof_id_type dof = _elems[j][i]->dof_number(T_var.sys().number(), T_var.number(), 0);
      T_sln.set(dof, _Ts[j][i]);
    }

  T_sln.close();
}

void
SimpleHeatStructure2::writeTextOutput()
{
  FILE * file = _sim.getTextOutputFile();

  fprintf(file, "Component = %s\n", name().c_str());
  fprintf(file, "Energy balance:\n");
  fprintf(file, "  Internal heating  [W]: %20.8e\n", _Q_internal);
  fprintf(file, "  Heat loss (left side)  [W]: %20.8e\n", _Q_out_left);
  fprintf(file, "  Heat loss (right side)  [W]: %20.8e\n", _Q_out_right);
  fprintf(file, "  Energy balance  [W]: %20.8e\n", _Q_internal + _Q_out_left + _Q_out_right);
  fprintf(file, "\n");
}

void
SimpleHeatStructure2::FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp)
{
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
    {
      int row = Ts_DOF(j, i);
      // SELF
      mnzp->addEntry(row, row);
      if (i > 0) // WEST
        mnzp->addEntry(row, Ts_DOF(j, i-1));
      if (i < _NW) // EAST
        mnzp->addEntry(row, Ts_DOF(j, i+1));
    }

  if (_pipe_left)
  {
    std::vector<OneDCell*> & fluid_cells = _pipe_left->getCells();
    for (unsigned j = 0; j < _NL; j++)
      mnzp->addEntry(Ts_DOF(j, 0), fluid_cells[j]->TDOF());
  }
  if (_pipe_right)
  {
    std::vector<OneDCell*> & fluid_cells = _pipe_right->getCells();
    for (unsigned j = 0; j < _NL; j++)
      mnzp->addEntry(Ts_DOF(j, _NW), fluid_cells[j]->TDOF());
  }
}
