#include "SystemBase.h" // for access solution
#include "TestComponent.h"

registerMooseObject("delphiApp", TestComponent);

InputParameters
TestComponent::validParams()
{
  InputParameters params = DelPhiComponent::validParams();

  params.addRequiredParam<std::vector<Real>>("position", "Origin (start) of the component");
  params.addRequiredParam<std::vector<Real>>("orientation", "Orientation vector of the component");
  params.addRequiredParam<Real>("length", "Length of the OneDComp");
  params.addRequiredParam<unsigned>("n_elems", "number of element in this OneDComp");

  return params;
}

TestComponent::TestComponent(const InputParameters & parameters)
  : DelPhiComponent(parameters),
    _length(getParam<Real>("length")),
    _n_elem(getParam<unsigned>("n_elems")),
    _dL(_length / _n_elem)
{
}

void
TestComponent::buildMesh()
{
  const std::vector<Real> & pos = getParam<std::vector<Real>>("position");
  const std::vector<Real> & dir = getParam<std::vector<Real>>("orientation");

  Point position = Point(pos[0], pos[1], pos[2]);
  RealVectorValue direction = VectorValue<Real>(dir[0], dir[1], dir[2]);
  direction = direction.unit();

  _subdomain_id = getNextSubdomainId();
  _subdomain_name = Moose::stringify(_subdomain_id);
  _mesh.setSubdomainName(_subdomain_id, _subdomain_name);

  // points
  Point p = position;
  for (unsigned i = 0; i < _n_elem + 1; i++)
  {
    Node * nd = _mesh.getMesh().add_point(p);
    _nodes.push_back(nd);
    p += _dL * direction;
  }

  // elems
  for (unsigned i = 0; i < _n_elem; i++)
  {
    Elem * elem = _mesh.getMesh().add_elem(new Edge2);
    elem->subdomain_id() = _subdomain_id;
    elem->set_node(0) = _nodes[i];
    elem->set_node(1) = _nodes[i + 1];
    _elems.push_back(elem);
  }
}

void
TestComponent::addExternalVariables()
{
  _n_DOFs = _n_elem;
  _T.resize(_n_elem, 0);
  _T_old.resize(_n_elem, 0);
  _sim.addMooseAuxVar("T_cell", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("T_analytical", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
}

void
TestComponent::updateSolution(double * u)
{
  for (unsigned i = 0; i < _T.size(); i++)
    _T[i] = u[i];
}

void
TestComponent::computeTranRes(double * res)
{
  // with rho = 1, cp = 1
  for (unsigned i = 0; i < _T.size(); i++)
    res[i] = (_T[i] - _T_old[i]) * _dL / _sim.dt();
}

void
TestComponent::computeSpatialRes(double * res)
{
  // with k = 1, q" = -k dT/dx
  double q0 = -2 * (_T[0] - 0) / _dL;
  res[0] -= q0;
  double qN = -2 * (0 - _T.back()) / _dL;
  res[_n_elem - 1] += qN;

  for (unsigned i = 1; i < _n_elem; i++)
  {
    double qi = -(_T[i] - _T[i - 1]) / _dL;
    res[i - 1] += qi;
    res[i] -= qi;
  }

  // source term
  for (unsigned i = 0; i < _n_elem; i++)
  {
    double x = _dL * (i + 0.5);
    double q_vol = M_PI * M_PI / (_length * _length) * sin(M_PI * x / _length);
    res[i] -= q_vol * _dL;
  }
}

void
TestComponent::onTimestepBegin()
{
}

void
TestComponent::onTimestepEnd()
{
  _T_old = _T;

  // output
  MooseVariableFieldBase & T_var = _sim.getVariable(0, "T_cell");
  NumericVector<Number> & T_sln = T_var.sys().solution();
  MooseVariableFieldBase & T_ana = _sim.getVariable(0, "T_analytical");
  NumericVector<Number> & T_ana_sln = T_ana.sys().solution();
  for (unsigned i = 0; i < _elems.size(); i++)
  {
    dof_id_type dof = _elems[i]->dof_number(T_var.sys().number(), T_var.number(), 0);
    T_sln.set(dof, _T[i]);

    dof = _elems[i]->dof_number(T_ana.sys().number(), T_ana.number(), 0);
    T_ana_sln.set(dof, sin(M_PI / _length * (i + 0.5) * _dL));
  }
}

void
TestComponent::FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp)
{
  // When loop, use int, not unsigned int; 'col' might go below 0, unsigned won't handle correctly
  for (int row = 0; row < int(_n_DOFs); row++)
    for (int col = row - 1; col <= row + 1; col++)
    {
      if (col >= 0 && col < int(_n_DOFs))
        mnzp->addEntry(row + _DOF_offset, col + _DOF_offset);
    }
}
