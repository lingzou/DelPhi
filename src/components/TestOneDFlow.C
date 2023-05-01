#include "SystemBase.h" // for access solution
#include "TestOneDFlow.h"
//#include "SinglePhaseFluidProperties.h"

registerMooseObject("delphiApp", TestOneDFlow);

InputParameters
TestOneDFlow::validParams()
{
  InputParameters params = DelPhiComponent::validParams();

  params.addRequiredParam<std::vector<Real>>("position", "Origin (start) of the component");
  params.addRequiredParam<std::vector<Real>>("orientation", "Orientation vector of the component");
  params.addRequiredParam<Real>("length", "Length of the OneDComp");
  params.addRequiredParam<unsigned>("n_elems", "number of element in this OneDComp");
  params.addRequiredParam<UserObjectName>("eos", "equation of states");

  return params;
}

TestOneDFlow::TestOneDFlow(const InputParameters & parameters)
  : DelPhiComponent(parameters),
    _length(getParam<Real>("length")),
    _n_elem(getParam<unsigned>("n_elems")),
    _dL(_length / _n_elem)
{
}

void
TestOneDFlow::buildMesh()
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
TestOneDFlow::addExternalVariables()
{
  _n_DOFs = (_n_elem) + (_n_elem) + (_n_elem + 1); // p + T + v

  _p.resize(_n_elem);
  _p_old.resize(_n_elem);
  _T.resize(_n_elem);
  _T_old.resize(_n_elem);
  _v.resize(_n_elem + 1);
  _v_old.resize(_n_elem + 1);

  _rho.resize(_n_elem);
  _rho_old.resize(_n_elem);
  _rho_edge.resize(_n_elem + 1);

  _h.resize(_n_elem);
  _h_old.resize(_n_elem);

  _mass_flux.resize(_n_elem + 1, 0.0);
  _enthalpy_flux.resize(_n_elem + 1, 0.0);

  _sim.addMooseAuxVar("p", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("T", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("rho", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("enthalpy", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("v", FEType(FIRST, LAGRANGE), {_subdomain_name});
}

void
TestOneDFlow::setupIC(double * u)
{
  // handle eos first
  const UserObjectName & uo_name = getParam<UserObjectName>("eos");
  const UserObject & uo = _sim.getUserObject<UserObject>(uo_name);
  if (dynamic_cast<const SinglePhaseFluidProperties *>(&uo) == nullptr)
    mooseError("cannot convert: " + uo_name);
  else
    _eos = dynamic_cast<const SinglePhaseFluidProperties *>(&uo);

  // setup initial conditions
  Real v_init = 0.5;
  Real p_init = 1e5;
  Real T_init = 300.0;
  unsigned idx = 0;
  for(unsigned i = 0; i < _n_elem + 1; i++)
  {
    _v[i] = v_init;
    u[idx++] = _v[i];

    if (i < _n_elem)
    {
      _p[i] = p_init;
      u[idx++] = _p[i];

      _T[i] = T_init;
      u[idx++] = _T[i];

      _rho[i] = _eos->rho_from_p_T(_p[i], _T[i]);
      _h[i] = _eos->h_from_p_T(_p[i], _T[i]);
    }
  }

  _p_old = _p;
  _T_old = _T;
  _v_old = _v;
  _rho_old = _rho;
  _h_old = _h;

  _rho_ref = _eos->rho_from_p_T(p_init, T_init);
  _rhoh_ref = _rho_ref * _eos->h_from_p_T(p_init, T_init);
}

void
TestOneDFlow::updateSolution(double * u)
{
  unsigned idx = 0;
  for(unsigned i = 0; i < _n_elem + 1; i++)
  {
    _v[i] = u[idx++];
    if (i < _n_elem)
    {
      _p[i] = u[idx++];
      _T[i] = u[idx++];
      _rho[i] = _eos->rho_from_p_T(_p[i], _T[i]);
      _h[i] = _eos->h_from_p_T(_p[i], _T[i]);
    }
  }
  //update edge values
  _rho_edge[0] = _rho[0];
  _rho_edge[_n_elem] = _rho[_n_elem - 1];
  for(unsigned i = 1; i < _n_elem; i++)
    _rho_edge[i] = 0.5 * (_rho[i-1] + _rho[i]);
}

void
TestOneDFlow::computeTranRes(double * res)
{
  unsigned idx = 0;
  for(unsigned i = 0; i < _n_elem + 1; i++)
  {
    res[idx++] = _rho_edge[i] * (_v[i] - _v_old[i]) / _sim.dt() / _rho_ref;
    if (i < _n_elem)
    {
      res[idx++] = (_rho[i] - _rho_old[i]) / _sim.dt() / _rho_ref;
      res[idx++] = (_rho[i] * _h[i] - _rho_old[i] * _h_old[i]) / _sim.dt() / _rhoh_ref;
    }
  }

  // for testing
  // Zero-out transient residual, to apply Dirichlet BC at the inlet later
  res[0] = 0.0;
}

void
TestOneDFlow::updateFluxes()
{
  // assume v_inlet (>0) and p_outlet condition
  Real rho_in = _eos->rho_from_p_T(_p[0], 300.0);
  Real h_in = _eos->h_from_p_T(_p[0], 300.0);
  Real rho_out = _eos->rho_from_p_T(1e5, 300.0);
  Real h_out = _eos->h_from_p_T(1e5, 300.0);

  // Upwind donor cell method for void fraction and mass balance equations
  _mass_flux[0] = (_v[0] > 0) ? _v[0] * rho_in : _v[0] * _rho[0];
  _enthalpy_flux[0] = (_v[0] > 0) ? _v[0] * rho_in * h_in : _v[0] * _rho[0] * _h[0];

  _mass_flux[_n_elem] = (_v[_n_elem] > 0) ? _v[_n_elem] * _rho[_n_elem-1] : _v[_n_elem] * rho_out;
  _enthalpy_flux[_n_elem] = (_v[_n_elem] > 0) ? _v[_n_elem] * _rho[_n_elem-1] * _h[_n_elem-1]
                         : _v[_n_elem] * rho_out * h_out;

  for(unsigned i = 1; i < _n_elem; i++)
  {
    _mass_flux[i] = (_v[i] > 0) ? _v[i] * _rho[i-1] : _v[i] * _rho[i];
    _enthalpy_flux[i] = (_v[i] > 0) ? _v[i] * _rho[i-1] * _h[i-1] : _v[i] * _rho[i] * _h[i];
  }
}

void
TestOneDFlow::computeSpatialRes(double * res)
{
  updateFluxes();

  // Momentum equations RHS
  Real f = 0.1;
  Real dh = 0.01;
  res[0] = _v[0] - 1.0; // Dirichlet at inlet

  for(unsigned i = 1; i < _n_elem + 1; i++) // loop on the remaining edges
  {
    // east and west velocities
    Real v_east = (i == _n_elem) ? _v[_n_elem] : _v[i+1];
    Real v_west = _v[i-1];
    Real dv_dx = (_v[i] > 0) ? (_v[i] - v_west) / _dL : (v_east - _v[i]) / _dL;

    // dp_dx term
    Real dp_dx = (i == _n_elem) ? (1e5 - _p[_n_elem-1])/_dL*2.0 : (_p[i] - _p[i-1])/_dL;

    // friction term
    Real fric  = 0.5 * f / dh * _rho_edge[i] * _v[i] * std::fabs(_v[i]);

    // assemble spatial terms
    res[3*i] = (_rho_edge[i] * _v[i] * dv_dx + dp_dx + fric) / _rho_ref;
  }

  // spatial terms for mass and energy equations
  for(unsigned i = 0; i < _n_elem; i++) //loop on elements
  {
    res[3*i+1] = (_mass_flux[i+1] - _mass_flux[i]) / _dL / _rho_ref;
    res[3*i+2] = ((_enthalpy_flux[i+1] - _enthalpy_flux[i]) / _dL - 1e5) / _rhoh_ref; // 1e3 source term
  }
}

void
TestOneDFlow::onTimestepBegin()
{
}

void
TestOneDFlow::onTimestepEnd()
{
  // save old solutions
  _p_old = _p;
  _T_old = _T;
  _v_old = _v;

  _rho_old = _rho;
  _h_old = _h;

  // output (element/cell value)
  MooseVariableFieldBase & T_var = _sim.getVariable(0, "T");
  NumericVector<Number> & T_sln = T_var.sys().solution();
  MooseVariableFieldBase & p_var = _sim.getVariable(0, "p");
  NumericVector<Number> & p_sln = p_var.sys().solution();
  MooseVariableFieldBase & rho_var = _sim.getVariable(0, "rho");
  NumericVector<Number> & rho_sln = rho_var.sys().solution();
  MooseVariableFieldBase & h_var = _sim.getVariable(0, "enthalpy");
  NumericVector<Number> & h_sln = h_var.sys().solution();
  for (unsigned i = 0; i < _elems.size(); i++)
  {
    dof_id_type dof = _elems[i]->dof_number(T_var.sys().number(), T_var.number(), 0);
    T_sln.set(dof, _T[i]);

    dof = _elems[i]->dof_number(p_var.sys().number(), p_var.number(), 0);
    p_sln.set(dof, _p[i]);

    dof = _elems[i]->dof_number(rho_var.sys().number(), rho_var.number(), 0);
    rho_sln.set(dof, _rho[i]);

    dof = _elems[i]->dof_number(h_var.sys().number(), h_var.number(), 0);
    h_sln.set(dof, _h[i]);
  }

  // output (node/edge value)
  MooseVariableFieldBase & v_var = _sim.getVariable(0, "v");
  NumericVector<Number> & v_sln = v_var.sys().solution();
  for (unsigned i = 0; i < _nodes.size(); i++)
  {
    dof_id_type dof = _nodes[i]->dof_number(v_var.sys().number(), v_var.number(), 0);
    v_sln.set(dof, _v[i]);
  }
}

void
TestOneDFlow::FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp)
{
  int n_Var = 3;
  for (int i = 0; i < _n_elem + 1; i++)
    for (int var = 0; var < n_Var; var++)
    {
      int i_dof = i * n_Var + var;
      for (int j_dof = (i - 2) * n_Var; j_dof < (i + 3) * n_Var; j_dof++)
      {
        if ((i_dof >= 0) && (i_dof < int(_n_DOFs)) && (j_dof >= 0) && (j_dof < int(_n_DOFs)))
          mnzp->addEntry(i_dof + _DOF_offset, j_dof + _DOF_offset);
      }
    }
}
