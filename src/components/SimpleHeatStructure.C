#include "libmesh/face_quad4.h" // for meshing
#include "SystemBase.h" // for access solution
#include "SimpleHeatStructure.h"

/**
 * This SimpleHeatStructure follows the cell-center finite volume method for heat conduction equation
 * in structured mesh.
 * As depicted, discretized temperature is arranged in the cell center (x).
 * At the same time, left/right surface temperatures are also added for easier handling of convective
 * heat transfer and thermal radiation.
 *
 *                   -------------------------------
 *                   |         |         |         |
 *                   o    x    |    x    |    x    o
 *                   |         |         |         |
 *                   -------------------------------
 *                   |         |         |         |
 *                   o    x    |    x    |    x    o
 *                   |         |         |         |
 *                   -------------------------------
 *                   |         |         |         |
 *                   o    x    |    x    |    x    o
 *                   |         |         |         |
 *                   -------------------------------
 */

registerMooseObject("delphiApp", SimpleHeatStructure);

InputParameters
SimpleHeatStructure::validParams()
{
  InputParameters params = OneDComponent::validParams();

  params.addRequiredParam<std::vector<Real>>("position", "Origin (start) of the component");
  params.addRequiredParam<std::vector<Real>>("orientation", "Orientation vector of the component");
  params.addRequiredParam<Real>("length", "Length of the SimpleHeatStructure");
  params.addRequiredParam<Real>("width", "Width of the SimpleHeatStructure");
  params.addRequiredParam<unsigned>("elem_number_axial", "number of element in the axial direction");
  params.addRequiredParam<unsigned>("elem_number_width", "number of element in the width/radial direction");
  params.addRequiredParam<UserObjectName>("solid", "solid property UserObject name");

  // initial condition
  params.addRequiredParam<Real>("Ts_init", "Initial Ts");

  // heat source
  params.addParam<Real>("heat_source_solid", 0.0, "volumetric heat source in solid");

  // boundary conditions
  params.addParam<Real>("T_bc_left", "User-specified temperature at left boundary");
  params.addParam<Real>("T_bc_right", "User-specified temperature at right boundary");

  params.addParam<std::string>("name_comp_left", "The name of the left flow component");
  params.addParam<Real>("HT_surface_area_density_left", "Aw at the left surface");
  params.addParam<Real>("Hw_left", "Hw at the left surface");
  params.addParam<std::string>("name_comp_right", "The name of the right flow component");
  params.addParam<Real>("HT_surface_area_density_right", "Aw at the right surface");
  params.addParam<Real>("Hw_right", "Hw at the right surface");
  return params;
}

SimpleHeatStructure::SimpleHeatStructure(const InputParameters & parameters)
  : OneDComponent(parameters),
    _length(getParam<Real>("length")),
    _width(getParam<Real>("width")),
    _NL(getParam<unsigned>("elem_number_axial")),
    _NW(getParam<unsigned>("elem_number_width")),
    _dL(_length / _NL),
    _dW(_width / _NW),
    _qv(getParam<Real>("heat_source_solid"))
{
  _T_bc_left = isParamValid("T_bc_left") ? getParam<Real>("T_bc_left") : 0.0;
  _T_bc_right = isParamValid("T_bc_right") ? getParam<Real>("T_bc_right") : 0.0;
  _depth = 1.0;
}

void
SimpleHeatStructure::buildMesh()
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
    for (unsigned i = 0; i < _NW + 1; i++) // loop on column in width/radial direction
    {
      Node * nd = _mesh.getMesh().add_point(p);
      _nodes[j].push_back(nd);
      p += _dW * dir_W;
    }
  }

  // elems
  _elems.resize(_NL);
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW; i++)
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
SimpleHeatStructure::setupIC(double * u)
{
  Real Ts_init = getParam<Real>("Ts_init");
  for (unsigned j = 0; j < _NL; j++)
  {
    _Tw_left[j] = Ts_init;
    u[_Tw_left_DOFs[j] - _DOF_offset] = Ts_init;
  
    for (unsigned i = 0; i < _NW; i++)
    {
      _Ts[j][i] = Ts_init;
      u[_Ts_DOFs[j][i] - _DOF_offset] = Ts_init;
    }

    _Tw_right[j] = Ts_init;
    u[_Tw_right_DOFs[j] - _DOF_offset] = Ts_init;
  }

  _Ts_old = _Ts;
  _Ts_oo = _Ts;

  _rho_ref = _solid->rho_from_T(Ts_init);
  _cp_ref = _solid->cp_from_T(Ts_init);
}

void
SimpleHeatStructure::addPhysicalModel()
{
  // _Ts is arranged as _Ts[Y/axial][X/radial], so when using one-dimensional assumption, the matrix is denser
  _Ts.resize(_NL);
  _Ts_DOFs.resize(_NL);
  for (unsigned j = 0; j < _Ts.size(); j++)
  {
    _Tw_left.push_back(0.0);
    _Ts[j].resize(_NW, 0.0);
    _Tw_right.push_back(0.0);

    _Tw_left_DOFs.push_back(j * (_NW + 2) + 0 + _DOF_offset);
    for (unsigned i = 0; i < _NW; i++)
      _Ts_DOFs[j].push_back(j * (_NW + 2) + i + 1 + _DOF_offset);
    _Tw_right_DOFs.push_back(j * (_NW + 2) + _NW + 1 + _DOF_offset);
  }

  _Ts_old = _Ts;
  _Ts_oo = _Ts;

  _n_DOFs = (_NW + 2) * _NL;

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

  _volume = _Ts;
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW; i++)
      _volume[j][i] = _dW * _dL * _depth;
}

void
SimpleHeatStructure::addExternalVariables()
{
  _sim.addMooseAuxVar("T", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
}

void
SimpleHeatStructure::updateSolution(double * u)
{
  for (unsigned j = 0; j < _NL; j++)
  {
    _Tw_left[j] = u[_Tw_left_DOFs[j] - _DOF_offset];
    for (unsigned i = 0; i < _NW ; i++)
      _Ts[j][i] = u[_Ts_DOFs[j][i] - _DOF_offset];
    _Tw_right[j] = u[_Tw_right_DOFs[j] - _DOF_offset];
  }
}

void
SimpleHeatStructure::computeTranRes(double * res)
{
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW; i++)
    {
      // Real rho = _solid->rho_from_T(_Ts[j][i]);
      // Real cp = _solid->cp_from_T(_Ts[j][i]);

      Real dT_dt = 0.0;
      if ((_sim.TS() == Moose::TI_IMPLICIT_EULER) || (_sim.timeStep() == 1)) // bdf1 or step 1 of bdf2
        dT_dt = (_Ts[j][i] - _Ts_old[j][i]) / _sim.dt();
      else //bdf2 and step > 1
        dT_dt = DELPHI::BDF2_ddt(_Ts[j][i], _Ts_old[j][i], _Ts_oo[j][i], _sim.dt(), _sim.dtOld());

      // option 1:
      // res[_Ts_DOFs[j][i] - _DOF_offset] = rho * cp * dT_dt / (_rho_ref * _cp_ref);
      // option 2:
      res[_Ts_DOFs[j][i] - _DOF_offset] = dT_dt;
    }
}

void
SimpleHeatStructure::computeSpatialRes(double * res)
{
  // loop on vertical interior faces
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW - 1; i++)
    {
      Real TW = _Ts[j][i];
      Real TE = _Ts[j][i+1];
      Real T_avg = 0.5 * (TW + TE);
      Real k = _solid->k_from_T(T_avg);
      Real q_flux = -k * (TE - TW) / _dW;

      // West
      res[_Ts_DOFs[j][i] - _DOF_offset] += q_flux * _dL * _depth;
      // East
      res[_Ts_DOFs[j][i+1] - _DOF_offset] -= q_flux * _dL * _depth;
    }

  // heat
  _Q_internal = 0.0;
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW; i++)
    {
      // Real x = _dW * (i + 0.5);
      //Real q_vol = M_PI * M_PI / (_width * _width) * sin(M_PI * x);
      _Q_internal += _qv * _volume[j][i];
      res[_Ts_DOFs[j][i] - _DOF_offset] -= _qv * _volume[j][i];
    }

  // Apply BCs
  _Q_out_left = 0.0;
  _Q_out_right = 0.0;
  if (isParamValid("T_bc_left"))
  {
    for (unsigned j = 0; j < _NL; j++)
    {
      Real q = -2.0 * _solid->k_from_T(_T_bc_left) * (_Ts[j][0] - _T_bc_left) / _dW * _dL * _depth;
      _Q_out_left += q;
      res[_Ts_DOFs[j].front() - _DOF_offset] -= q;
      res[_Tw_left_DOFs[j] - _DOF_offset] = _Tw_left[j] - _T_bc_left;
    }
  }
  else if (_pipe_left)
  {
    std::vector<OneDCell*> & fluid_cells = _pipe_left->getCells();
    for (unsigned j = 0; j < _NL; j++)
    {
      Real T_fluid = fluid_cells[j]->T();
      Real q_flux_out = _hw_left * (_Tw_left[j] - T_fluid) * _dL * _depth;
      _Q_out_left -= q_flux_out;
      res[_Ts_DOFs[j].front() - _DOF_offset] += q_flux_out;

      Real q = -2.0 * _solid->k_from_T(_Tw_left[j]) * (_Ts[j][0] - _Tw_left[j]) / _dW * _dL * _depth;
      res[_Tw_left_DOFs[j] - _DOF_offset] = (q_flux_out + q) / (_hw_left * _dL * _depth);
    }
  }

  if (isParamValid("T_bc_right"))
    for (unsigned j = 0; j < _NL; j++)
    {
      Real q = -2.0 * _solid->k_from_T(_T_bc_right) * (_T_bc_right - _Ts[j][_NW-1]) / _dW * _dL * _depth;
      _Q_out_right -= q;
      res[_Ts_DOFs[j].back() - _DOF_offset] += q;
      res[_Tw_right_DOFs[j] - _DOF_offset] = _Tw_right[j] - _T_bc_right;
    }

  // scaling
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW; i++)
    {
      // option 1:
      // res[_Ts_DOFs[j][i] - _DOF_offset] /= (_volume[j][i] * _rho_ref * _cp_ref);

      // option 2:
      Real rho = _solid->rho_from_T(_Ts[j][i]);
      Real cp = _solid->cp_from_T(_Ts[j][i]);
      res[_Ts_DOFs[j][i] - _DOF_offset] /= (_volume[j][i] * rho * cp);
    }
}

void
SimpleHeatStructure::onTimestepBegin()
{
}

void
SimpleHeatStructure::onTimestepEnd()
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
SimpleHeatStructure::writeTextOutput()
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
SimpleHeatStructure::FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp)
{
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW; i++)
    {
      int row = _Ts_DOFs[j][i]; // Ts_DOF(j, i);
      // SELF
      mnzp->addEntry(row, row);
      
      if (i > 0) // WEST
        mnzp->addEntry(row, _Ts_DOFs[j][i-1]);

      if (i < _NW-1) // EAST
        mnzp->addEntry(row, _Ts_DOFs[j][i+1]);
    }

  for (unsigned j = 0; j < _NL; j++)
  {
    mnzp->addEntry(_Tw_left_DOFs[j], _Tw_left_DOFs[j]);
    mnzp->addEntry(_Tw_left_DOFs[j], _Ts_DOFs[j].front());
    mnzp->addEntry(_Ts_DOFs[j].front(), _Tw_left_DOFs[j]);

    mnzp->addEntry(_Tw_right_DOFs[j], _Tw_right_DOFs[j]);
    mnzp->addEntry(_Tw_right_DOFs[j], _Ts_DOFs[j].back());
    mnzp->addEntry(_Ts_DOFs[j].back(), _Tw_right_DOFs[j]);
  }

  if (_pipe_left)
  {
    std::vector<OneDCell*> & fluid_cells = _pipe_left->getCells();
    for (unsigned j = 0; j < _NL; j++)
    {
      mnzp->addEntry(_Tw_left_DOFs[j], fluid_cells[j]->TDOF());
      mnzp->addEntry(_Ts_DOFs[j].front(), fluid_cells[j]->TDOF());
    }
  }
  if (_pipe_right)
  {
    std::vector<OneDCell*> & fluid_cells = _pipe_right->getCells();
    for (unsigned j = 0; j < _NL; j++)
    {
      mnzp->addEntry(_Tw_right_DOFs[j], fluid_cells[j]->TDOF());
      mnzp->addEntry(_Ts_DOFs[j].back(), fluid_cells[j]->TDOF());
    }
  }
}
