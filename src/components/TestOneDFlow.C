#include "SystemBase.h" // for access solution
#include "TestOneDFlow.h"
#include "OneDFlowModel.h"
//#include "SinglePhaseFluidProperties.h"

registerMooseObject("delphiApp", TestOneDFlow);

InputParameters
TestOneDFlow::validParams()
{
  InputParameters params = OneDComponent::validParams();

  params.addRequiredParam<std::vector<Real>>("position", "Origin (start) of the component");
  params.addRequiredParam<std::vector<Real>>("orientation", "Orientation vector of the component");
  params.addRequiredParam<Real>("length", "Length of the OneDComp");
  params.addRequiredParam<unsigned>("n_elems", "number of element in this OneDComp");
  params.addRequiredParam<UserObjectName>("eos", "equation of states");

  return params;
}

TestOneDFlow::TestOneDFlow(const InputParameters & parameters)
  : OneDComponent(parameters),
    _length(getParam<Real>("length")),
    _n_elem(getParam<unsigned>("n_elems")),
    _dL(_length / _n_elem)
{
}

TestOneDFlow::~TestOneDFlow()
{
  for(auto& cell : _cells)   delete cell;
  for(auto& edge : _edges)   delete edge;
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

  // handle eos first
  const UserObjectName & uo_name = getParam<UserObjectName>("eos");
  const UserObject & uo = _sim.getUserObject<UserObject>(uo_name);
  if (dynamic_cast<const SinglePhaseFluidProperties *>(&uo) == nullptr)
    mooseError("cannot convert: " + uo_name);
  else
    _eos = dynamic_cast<const SinglePhaseFluidProperties *>(&uo);

  _cells.resize(_n_elem, NULL);
  for (unsigned i = 0; i < _n_elem; i++)
  {
    InputParameters pars = emptyInputParameters();
    pars.set<std::string>("name") = name() + ":cell_" + std::to_string(i);
    pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
    pars.set<Real>("dL") = _dL;
    _cells[i] = new CellBase(pars);
  }

  _edges.resize(_n_elem + 1, NULL);
  { // inlet vEdge
    InputParameters pars = emptyInputParameters();
    pars.set<std::string>("name") = name() + ":inlet_vEdge";
    pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
    pars.set<CellBase *>("west_cell") = NULL;
    pars.set<CellBase *>("east_cell") = _cells[0];
    pars.set<Real>("v_bc") = 1.0;
    pars.set<Real>("T_bc") = 300.0;
    _edges[0] = new vBCEdge(pars);
  }
  for (unsigned i = 1; i < _n_elem; i++) // int edge only, bc edges will be handled later
  {
    InputParameters pars = emptyInputParameters();
    pars.set<std::string>("name") = name() + ":int_edge_" + std::to_string(i);
    pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
    pars.set<CellBase *>("west_cell") = _cells[i-1];
    pars.set<CellBase *>("east_cell") = _cells[i];
    _edges[i] = new IntEdge(pars);
  }
  { // outlet pEdge
    InputParameters pars = emptyInputParameters();
    pars.set<std::string>("name") = name() + ":outlet_pEdge";
    pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
    pars.set<CellBase *>("west_cell") = _cells[_n_elem-1];
    pars.set<CellBase *>("east_cell") = NULL;
    pars.set<Real>("p_bc") = 1.0e5;
    pars.set<Real>("T_bc") = 300.0;
    _edges[_n_elem] = new pBCEdge(pars);
  }

  for (unsigned i = 0; i < _edges.size(); i++)
    _edges[i]->setDOF(_DOF_offset + 3 * i);
  for (unsigned i = 0; i < _cells.size(); i++)
    _cells[i]->setDOF(_DOF_offset + 3 * i + 1, _DOF_offset + 3 * i + 2);

  _sim.addMooseAuxVar("p", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("T", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("rho", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("enthalpy", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("v", FEType(FIRST, LAGRANGE), {_subdomain_name});
}

void
TestOneDFlow::setExtendedNeighbors()
{
  for(auto& cell : _cells)  cell->setExtendedNeighborCells();
  for(auto& edge : _edges)  edge->setExtendedNeighborEdges();
}

void
TestOneDFlow::setupIC(double * u)
{
  // setup initial conditions
  Real v_init = 0.5;
  Real p_init = 1e5;
  Real T_init = 300.0;
  unsigned idx = 0;
  for(unsigned i = 0; i < _n_elem + 1; i++)
  {
    u[idx++] = v_init;
    _edges[i]->initialize(v_init);

    if (i < _n_elem)
    {
      u[idx++] = p_init;
      u[idx++] = T_init;
      _cells[i]->initialize(p_init, T_init);
    }
  }

  _rho_ref = _eos->rho_from_p_T(p_init, T_init);
  _rhoh_ref = _rho_ref * _eos->h_from_p_T(p_init, T_init);
}

void
TestOneDFlow::updateSolution(double * u)
{
  unsigned idx = 0;
  for(unsigned i = 0; i < _n_elem + 1; i++)
  {
    _edges[i]->updateSolution(u[idx++]);
    if (i < _n_elem)
    {
      Real p = u[idx++];
      Real T = u[idx++];
      _cells[i]->updateSolution(p, T);
    }
  }

  vBCEdge * inlet_edge = dynamic_cast<vBCEdge *>(_edges.front());
  if (inlet_edge)
    inlet_edge->updateGhostPressure(_cells.front()->p()); // we can do linear projection
  vBCEdge * outlet_edge = dynamic_cast<vBCEdge *>(_edges.back());
  if (outlet_edge)
    outlet_edge->updateGhostPressure(_cells.back()->p()); // we can do linear projection
}

void
TestOneDFlow::computeTranRes(double * res)
{
  unsigned idx = 0;
  for(unsigned i = 0; i < _n_elem + 1; i++)
  {
    res[idx++] = _edges[i]->rho_edge() * _edges[i]->dv_dt(_sim.dt()) / _rho_ref;

    if (i < _n_elem)
    {
      res[idx++] = (_cells[i]->rho() - _cells[i]->rho_o()) / _sim.dt() / _rho_ref;
      res[idx++] = (_cells[i]->rhoh() - _cells[i]->rhoh_o()) / _sim.dt() / _rhoh_ref;
    }
  }
}

void
TestOneDFlow::computeSpatialRes(double * res)
{
  // Momentum equations RHS
  Real f = 0.1;
  Real dh = 0.01;
  for(unsigned i = 0; i < _n_elem + 1; i++) // loop on edges
  {
    Real rho_edge = _edges[i]->rho_edge();
    Real v = _edges[i]->v();

    Real dv_dx = _edges[i]->dv_dx();
    Real dp_dx = _edges[i]->dp_dx();
    // friction term
    Real fric = 0.5 * f / dh * rho_edge * v * std::fabs(v);

    // assemble spatial terms
    res[3*i] = (rho_edge * v * dv_dx + dp_dx + fric) / _rho_ref;
  }

  // apply Dirichlet BC for inlet/outlet v (if applicable)
  vBCEdge * inlet_edge = dynamic_cast<vBCEdge *>(_edges.front());
  if (inlet_edge)
    res[0] = inlet_edge->computeDirichletBCResidual();
  vBCEdge * outlet_edge = dynamic_cast<vBCEdge *>(_edges.back());
  if (outlet_edge)
    res[_n_DOFs-1] = outlet_edge->computeDirichletBCResidual();

  // spatial terms for mass and energy equations
  for(unsigned i = 0; i < _n_elem; i++) //loop on elements
  {
    EdgeBase * w_edge = _cells[i]->wEdge();
    EdgeBase * e_edge = _cells[i]->eEdge();

    res[3*i+1] = (e_edge->mass_flux() - w_edge->mass_flux()) / _dL / _rho_ref;
    res[3*i+2] = ((e_edge->enthalpy_flux() - w_edge->enthalpy_flux()) / _dL - 1e5) / _rhoh_ref; // 1e3 source term
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
  for(auto& cell : _cells)  cell->saveOldSlns();
  for(auto& edge : _edges)  edge->saveOldSlns();

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
    T_sln.set(dof, _cells[i]->T());

    dof = _elems[i]->dof_number(p_var.sys().number(), p_var.number(), 0);
    p_sln.set(dof, _cells[i]->p());

    dof = _elems[i]->dof_number(rho_var.sys().number(), rho_var.number(), 0);
    rho_sln.set(dof, _cells[i]->rho());

    dof = _elems[i]->dof_number(h_var.sys().number(), h_var.number(), 0);
    h_sln.set(dof, _cells[i]->h());
  }

  // output (node/edge value)
  MooseVariableFieldBase & v_var = _sim.getVariable(0, "v");
  NumericVector<Number> & v_sln = v_var.sys().solution();
  for (unsigned i = 0; i < _nodes.size(); i++)
  {
    dof_id_type dof = _nodes[i]->dof_number(v_var.sys().number(), v_var.number(), 0);
    v_sln.set(dof, _edges[i]->v());
  }
}

void
TestOneDFlow::FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp)
{
  for(auto& cell : _cells)
  {
    mnzp->addRow(cell->pDOF(), cell->getConnectedDOFs());
    mnzp->addRow(cell->TDOF(), cell->getConnectedDOFs());
  }
  for(auto& edge : _edges)
  {
    mnzp->addRow(edge->vDOF(), edge->getConnectedDOFs());
  }
}
