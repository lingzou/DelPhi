#include "SystemBase.h" // for access solution
#include "TestOneDFlow.h"
#include "OneDFlowModel.h"
#include "DelPhiTypes.h"

registerMooseObject("delphiApp", TestOneDFlow);

InputParameters
TestOneDFlow::validParams()
{
  InputParameters params = OneDComponent::validParams();

  params.addParam<unsigned>("order", 1, "order for spatial discretization");

  params.addRequiredParam<std::vector<Real>>("position", "Origin (start) of the component");
  params.addRequiredParam<std::vector<Real>>("orientation", "Orientation vector of the component");
  params.addRequiredParam<Real>("length", "Length of the OneDComp");
  params.addRequiredParam<unsigned>("n_elems", "number of element in this OneDComp");
  params.addRequiredParam<UserObjectName>("eos", "equation of states");

  params.addRequiredParam<Real>("initial_P", "Initial value for pressure");
  params.addRequiredParam<Real>("initial_V", "Initial value for velocity");
  params.addRequiredParam<Real>("initial_T", "Initial value for temperature");

  params.addRequiredParam<Real>("A", "Flow area");
  params.addRequiredParam<Real>("Dh", "Hydraulic diameter");

  // wall friction input parameters
  params.addParam<Real>("f", "friction factor");
  params.addParam<std::string>("WF_user_option", "User-provided wall friction option");
  params.addParam<Real>("Welander_constant", "Welander constant value used in Welander cases");

  params.addParam<Real>("heat_source", 0, "volumetric heat source");

  return params;
}

TestOneDFlow::TestOneDFlow(const InputParameters & parameters)
  : OneDComponent(parameters),
    _order(getParam<unsigned>("order")),
    _length(getParam<Real>("length")),
    _n_elem(getParam<unsigned>("n_elems")),
    _dL(_length / _n_elem),
    _flow_area(getParam<Real>("A")),
    _dh(getParam<Real>("Dh")),
    _f(getParam<Real>("f")),
    _qv(getParam<Real>("heat_source"))
{
  if (isParamValid("f") && isParamValid("WF_user_option"))
    mooseError("'f' and 'WF_user_option' cannot be both specified.");
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
TestOneDFlow::addPhysicalModel()
{
  _n_DOFs = (_n_elem) + (_n_elem) + (_n_elem + 1); // p + T + v

  // handle eos first
  _eos = _sim.getSinglePhaseEOS(getParam<UserObjectName>("eos"));

  _cells.resize(_n_elem, NULL);
  for (unsigned i = 0; i < _n_elem; i++)
  {
    InputParameters pars = emptyInputParameters();
    pars.set<DelPhiSimulation *>("_sim") = &_sim;
    pars.set<std::string>("name") = name() + ":cell_" + std::to_string(i);
    pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
    pars.set<Real>("Dh") = getParam<Real>("Dh");
    pars.set<Real>("dL") = _dL;
    if (isParamValid("f"))
    {
      pars.set<DELPHI::WallFrictionModel>("WF_option") = DELPHI::CONST_FRICTION;
      pars.set<Real>("f") = getParam<Real>("f");
    }
    else if (isParamValid("WF_user_option"))
    {
      DELPHI::WallFrictionModel wf_option = DELPHI::stringToEnum<DELPHI::WallFrictionModel>(getParam<std::string>("WF_user_option"));
      pars.set<DELPHI::WallFrictionModel>("WF_option") = wf_option;
      if (wf_option == DELPHI::WELANDER)
        pars.set<Real>("Welander_constant") = getParam<Real>("Welander_constant");
    }
    else
      mooseError("To be implemented");
    _cells[i] = new OneDCell(pars);
    _cells[i]->setDOF(_DOF_offset + 3 * i + 1, _DOF_offset + 3 * i + 2);
  }

  _edges.resize(_n_elem + 1, NULL);
  for (unsigned i = 1; i < _n_elem; i++) // int edge only, bc edges will be handled later
  {
    InputParameters pars = emptyInputParameters();
    pars.set<DelPhiSimulation *>("_sim") = &_sim;
    pars.set<std::string>("name") = name() + ":int_edge_" + std::to_string(i);
    pars.set<const SinglePhaseFluidProperties *>("eos") = _eos;
    pars.set<CellBase *>("west_cell") = _cells[i-1];
    pars.set<CellBase *>("east_cell") = _cells[i];
    _edges[i] = new IntEdge(pars);
    _edges[i]->setDOF(_DOF_offset + 3 * i);
  }
}

void
TestOneDFlow::addExternalVariables()
{
  _sim.addMooseAuxVar("p", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("T", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("rho", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("enthalpy", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("v", FEType(FIRST, LAGRANGE), {_subdomain_name});
  // _sim.addMooseAuxVar("T_edge", FEType(FIRST, LAGRANGE), {_subdomain_name});
}

void
TestOneDFlow::setBoundaryEdge(DELPHI::EEndType end, EdgeBase* edge)
{
  if (end == DELPHI::IN)
  {
    _edges.front() = edge;
    (_edges.front())->setDOF(_DOF_offset);
  }
  else if (end == DELPHI::OUT)
  {
    _edges.back() = edge;
    (_edges.back())->setDOF(_DOF_offset + 3 * _n_elem);
  }
  else
    mooseError("error");
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
  Real v_init = getParam<Real>("initial_V");
  Real p_init = getParam<Real>("initial_P");
  Real T_init = getParam<Real>("initial_T");
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
TestOneDFlow::highOrderReconstruction()
{
  if (_order == 1) return;
  if (_n_elem < 2) return;

  boundaryEdge * inletEdge = dynamic_cast<boundaryEdge*>(_edges.front());
  Real T_in = 0.0;
  if (inletEdge)
    T_in = (inletEdge->v() > 0.0) ? inletEdge->T_bc() : _cells.front()->T();
  else
    T_in = _cells.front()->wCell()->T();

  _cells[0]->linearReconstruction(_cells[0]->p(), T_in, _cells[1]->p(), _cells[1]->T());

  boundaryEdge * outletEdge = dynamic_cast<boundaryEdge*>(_edges.back());
  Real T_out = 0.0;
  if (outletEdge)
    T_out = (outletEdge->v() > 0.0) ? _cells.back()->T() : outletEdge->T_bc();
  else
    T_out = _cells.back()->eCell()->T();

  _cells[_n_elem-1]->linearReconstruction(_cells[_n_elem-2]->p(), _cells[_n_elem-2]->T(), _cells[_n_elem-1]->p(), T_out);

  for (unsigned i = 1; i < _n_elem - 1; i++)
    _cells[i]->linearReconstruction(_cells[i-1]->p(), _cells[i-1]->T(), _cells[i+1]->p(), _cells[i+1]->T());
}

void
TestOneDFlow::computeHelperVariables()
{
  for(auto& edge : _edges)
    edge->computeFluxes();

  for(auto& cell : _cells)
    cell->computeHelperVariables();
}

void
TestOneDFlow::computeTranRes(double * res)
{
  unsigned idx = 0;
  for(unsigned i = 0; i < _n_elem + 1; i++)
  {
    res[idx++] = _edges[i]->rho_edge() * _edges[i]->dv_dt() / _rho_ref;

    if (i < _n_elem)
    {
      res[idx++] = _cells[i]->drho_dt() / _rho_ref;
      res[idx++] = _cells[i]->drhoh_dt() / _rhoh_ref;
    }
  }
}

void
TestOneDFlow::computeSpatialRes(double * res)
{
  // Momentum equations RHS
  for(unsigned i = 0; i < _n_elem + 1; i++) // loop on edges
  {
    Real rho_edge = _edges[i]->rho_edge();
    Real v = _edges[i]->v();

    Real dv_dx = _edges[i]->dv_dx();
    Real dp_dx = _edges[i]->dp_dx();
    Real fric = _edges[i]->dp_dx_friction();

    // assemble spatial terms
    res[3*i] = (rho_edge * v * dv_dx + dp_dx + fric) / _rho_ref;
  }

  _edges.front()->applyDirichletBC(res[0]);
  _edges.back()->applyDirichletBC(res[_n_DOFs-1]);

  // spatial terms for mass and energy equations
  for(unsigned i = 0; i < _n_elem; i++) //loop on elements
  {
    EdgeBase * w_edge = _cells[i]->wEdge();
    EdgeBase * e_edge = _cells[i]->eEdge();

    res[3*i+1] = (e_edge->mass_flux() - w_edge->mass_flux()) / _dL / _rho_ref;
    res[3*i+2] = ((e_edge->enthalpy_flux() - w_edge->enthalpy_flux()) / _dL - _qv) / _rhoh_ref; // 1e3 source term
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
  T_sln.close();
  p_sln.close();
  rho_sln.close();
  h_sln.close();

  // output (node/edge value)
  MooseVariableFieldBase & v_var = _sim.getVariable(0, "v");
  NumericVector<Number> & v_sln = v_var.sys().solution();
  // MooseVariableFieldBase & T_edge_var = _sim.getVariable(0, "T_edge");
  // NumericVector<Number> & T_edge_sln = T_edge_var.sys().solution();
  for (unsigned i = 0; i < _nodes.size(); i++)
  {
    dof_id_type dof = _nodes[i]->dof_number(v_var.sys().number(), v_var.number(), 0);
    v_sln.set(dof, _edges[i]->v());

    // dof = _nodes[i]->dof_number(T_edge_var.sys().number(), T_edge_var.number(), 0);
    // T_edge_sln.set(dof, _edges[i]->T_edge());
  }
  v_sln.close();
  // T_edge_sln.close();
}

void
TestOneDFlow::writeTextOutput()
{
  FILE * file = _sim.getTextOutputFile();

  fprintf(file, "Component = %s\n", name().c_str());
  fprintf(file, "%20s%20s%20s%20s%20s%20s\n", "Cell", "x", "p", "T", "rho", "h");
  for (unsigned i = 0; i < _cells.size(); i++)
    fprintf(file, "%20s%20.8e%20.8e%20.8e%20.8e%20.8e\n", _cells[i]->name().c_str(), _dL * (i + 0.5), _cells[i]->p(), _cells[i]->T(), _cells[i]->rho(), _cells[i]->h());

  fprintf(file, "\n%20s%20s%20s%20s%20s%20s\n", "Edge", "x", "v", "T", "mass_flow_rate", "rho_u_h_A");
  for (unsigned i = 0; i < _edges.size(); i++)
    fprintf(file, "%20s%20.8e%20.8e%20.8e%20.8e%20.8e\n", _edges[i]->name().c_str(),
                                                          _dL * i,
                                                          _edges[i]->v(),
                                                          _edges[i]->T_edge(),
                                                          _edges[i]->mass_flux() * _flow_area,
                                                          _edges[i]->enthalpy_flux() * _flow_area);
  fprintf(file, "\n");
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
