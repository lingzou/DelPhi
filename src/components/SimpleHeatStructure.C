#include "libmesh/face_quad4.h" // for meshing
#include "SystemBase.h" // for access solution
#include "SimpleHeatStructure.h"

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
  return params;
}

SimpleHeatStructure::SimpleHeatStructure(const InputParameters & parameters)
  : OneDComponent(parameters),
    _length(getParam<Real>("length")),
    _width(getParam<Real>("width")),
    _NL(getParam<unsigned>("elem_number_axial")),
    _NW(getParam<unsigned>("elem_number_width")),
    _dL(_length / _NL),
    _dW(_width / _NW)
{
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
SimpleHeatStructure::setupIC(double * u)
{
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
    {
      _Ts[j][i] = 0.0;
      u[j * (_NW + 1) + i] = 0.0;
    }

  _Ts_old = _Ts;
  _Ts_oo = _Ts;
}

void
SimpleHeatStructure::addPhysicalModel()
{
  // _Ts is arranged as _Ts[Y/axial][X/radial], so when using one-dimensional assumption, the matrix is denser
  _Ts.resize(_NL);
  for (unsigned i = 0; i < _Ts.size(); i++)
    _Ts[i].resize(_NW + 1, 0.0);

  _Ts_old = _Ts;
  _Ts_oo = _Ts;

  _volume = _Ts;
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
      _volume[j][i] = (i == 0 || i == _NW) ? 0.5 * _dW * _dL : _dW * _dL;

  _n_DOFs = (_NW + 1) * _NL;

  _solid = &(_sim.getUserObject<ThermalSolidProperties>(getParam<UserObjectName>("solid")));
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
    for (unsigned i = 0; i < _NW + 1; i++)
      _Ts[j][i] = u[j * (_NW + 1) + i];
}

void
SimpleHeatStructure::computeTranRes(double * res)
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

  // west and east boundary DirichletBC
  for (unsigned j = 0; j < _NL; j++)
  {
    res[j * (_NW + 1) + 0] = 0;
    res[j * (_NW + 1) + _NW] = 0;
  }
}

void
SimpleHeatStructure::computeSpatialRes(double * res)
{
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
    }

  // heat
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
    {
      Real x = _dW * i;
      Real q_vol = M_PI * M_PI / (_width * _width) * sin(M_PI * x);
      res[j * (_NW + 1) + i] -= q_vol * _volume[j][i];
    }

  // west and east boundary DirichletBC
  for (unsigned j = 0; j < _NL; j++)
  {
    res[j * (_NW + 1) + 0] = _Ts[j][0] - 0.0;
    res[j * (_NW + 1) + _NW] = _Ts[j][_NW] - 0.0;
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
SimpleHeatStructure::FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp)
{
  for (unsigned j = 0; j < _NL; j++)
    for (unsigned i = 0; i < _NW + 1; i++)
    {
      int row = j * (_NW + 1) + i;
      // SELF
      mnzp->addEntry(row + _DOF_offset, row + _DOF_offset);
      if (i > 0) // WEST
        mnzp->addEntry(row + _DOF_offset, row - 1 + _DOF_offset);
      if (i < _NW) // EAST
        mnzp->addEntry(row + _DOF_offset, row + 1 + _DOF_offset);
    }
}
