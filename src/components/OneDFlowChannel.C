#include "SystemBase.h" // for access solution
#include "OneDFlowChannel.h"
#include "SimpleHeatStructure.h"
#include "OneDFlowModel.h"
#include "DelPhiTypes.h"

registerMooseObject("delphiApp", OneDFlowChannel);

InputParameters
OneDFlowChannel::validParams()
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

  params.addParam<Real>("Tw", "A user-specified constant wall temperature");
  params.addParam<Real>("Hw", 0, "A user-specified constant wall heat transfer coefficient");
  params.addParam<Real>("HT_surface_area_density", 0, "A user-specified constant wall heat transfer area density");
  params.addParam<Real>("heat_source", 0, "volumetric heat source");

  return params;
}

OneDFlowChannel::OneDFlowChannel(const InputParameters & parameters)
  : OneDComponent(parameters),
    _order(getParam<unsigned>("order")),
    _length(getParam<Real>("length")),
    _n_elem(getParam<unsigned>("n_elems")),
    _dL(_length / _n_elem),
    _flow_area(getParam<Real>("A")),
    _dh(getParam<Real>("Dh")),
    _qv(getParam<Real>("heat_source")),
    _has_Tw(isParamValid("Tw")),
    _Tw(_has_Tw ? getParam<Real>("Tw") : 0.0),
    _hw(getParam<Real>("Hw")),
    _aw(getParam<Real>("HT_surface_area_density"))
{
  if (isParamValid("f") && isParamValid("WF_user_option"))
    mooseError("'f' and 'WF_user_option' cannot be both specified.");
}

OneDFlowChannel::~OneDFlowChannel()
{
  for(auto& cell : _cells)   delete cell;
  for(auto& edge : _edges)   delete edge;
}

void
OneDFlowChannel::buildMesh()
{
  const std::vector<Real> & pos = getParam<std::vector<Real>>("position");
  const std::vector<Real> & dir = getParam<std::vector<Real>>("orientation");

  Point position = Point(pos[0], pos[1], pos[2]);
  RealVectorValue direction = VectorValue<Real>(dir[0], dir[1], dir[2]);
  if (direction.norm() < 1e-16)
    mooseError("'orientation' cannot be a zero vector.");
  direction = direction.unit();

  _gL = direction * _sim.getParam<RealVectorValue>("gravity");

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
OneDFlowChannel::addPhysicalModel()
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
    pars.set<Real>("gL") = _gL;
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
OneDFlowChannel::addExternalVariables()
{
  _sim.addMooseAuxVar("p", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("T", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("rho", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("enthalpy", FEType(CONSTANT, MONOMIAL), {_subdomain_name});
  _sim.addMooseAuxVar("v", FEType(FIRST, LAGRANGE), {_subdomain_name});
  // _sim.addMooseAuxVar("T_edge", FEType(FIRST, LAGRANGE), {_subdomain_name});
}

void
OneDFlowChannel::setBoundaryEdge(DELPHI::EEndType end, EdgeBase* edge)
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
OneDFlowChannel::setExtendedNeighbors()
{
  for(auto& cell : _cells)  cell->setExtendedNeighborCells();
  for(auto& edge : _edges)  edge->setExtendedNeighborEdges();
}

void
OneDFlowChannel::setupIC(double * u)
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
OneDFlowChannel::updateSolution(double * u)
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
OneDFlowChannel::highOrderReconstruction()
{
  if (_order == 1) return;
  if (_n_elem < 2) return;

  boundaryEdge * inletEdge = dynamic_cast<boundaryEdge*>(_edges.front());
  Real T_in = 0.0;
  Real p_in = 0.0;
  if (inletEdge != NULL) // this is a boundary (connected to TDJ/TDV)
  {
    T_in = (inletEdge->v() > 0.0) ? inletEdge->T_bc() : _cells.front()->T();
    p_in = inletEdge->p_bc();
    _cells[0]->linearReconstruction(p_in, T_in, _cells[1]->p(), _cells[1]->T());
  }
  else
    _cells[0]->linearReconstructionIrregularMesh();

  boundaryEdge * outletEdge = dynamic_cast<boundaryEdge*>(_edges.back());
  Real T_out = 0.0;
  Real p_out = 0.0;
  if (outletEdge != NULL)
  {
    T_out = (outletEdge->v() > 0.0) ? _cells.back()->T() : outletEdge->T_bc();
    p_out = outletEdge->p_bc();
    _cells[_n_elem-1]->linearReconstruction(_cells[_n_elem-2]->p(), _cells[_n_elem-2]->T(), p_out, T_out);
  }
  else
    _cells.back()->linearReconstructionIrregularMesh();

  // remaining interior OneDCells
  for (unsigned i = 1; i < _n_elem - 1; i++)
    _cells[i]->linearReconstruction(_cells[i-1]->p(), _cells[i-1]->T(), _cells[i+1]->p(), _cells[i+1]->T());
}

void
OneDFlowChannel::computeHelperVariables()
{
  for(auto& edge : _edges)
    edge->computeFluxes();

  for(auto& cell : _cells)
    cell->computeHelperVariables();
}

void
OneDFlowChannel::computeTranRes(double * res)
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
OneDFlowChannel::computeSpatialRes(double * res)
{
  // Momentum equations RHS
  for(unsigned i = 0; i < _n_elem + 1; i++) // loop on edges
  {
    Real rho_edge = _edges[i]->rho_edge();
    Real v = _edges[i]->v();

    Real dv_dx = _edges[i]->dv_dx();
    Real dp_dx = _edges[i]->dp_dx();
    Real fric = _edges[i]->dp_dx_friction();
    Real gravity = _edges[i]->gravity();

    // assemble spatial terms
    res[3*i] = (rho_edge * v * dv_dx + dp_dx + fric - gravity) / _rho_ref;
  }

  _edges.front()->applyDirichletBC(res[0]);
  _edges.back()->applyDirichletBC(res[_n_DOFs-1]);

  // spatial terms for mass and energy equations
  for(unsigned i = 0; i < _n_elem; i++) //loop on elements
  {
    EdgeBase * w_edge = _cells[i]->wEdge();
    EdgeBase * e_edge = _cells[i]->eEdge();

    res[3*i+1] = (e_edge->mass_flux() - w_edge->mass_flux()) / _dL / _rho_ref;
    Real res_energy = (e_edge->enthalpy_flux() - w_edge->enthalpy_flux()) / _dL - _qv;

    if (_has_Tw)
      res_energy -= _hw * _aw * (_Tw - _cells[i]->T());

    if (_hs)
    {
      std::vector<std::vector<Real>> & Ts = _hs->Ts();
      Real Tw = (_hs_side == 0) ? Ts[i].front() : Ts[i].back();
      res_energy -= 1000.0 * 200.0 * (Tw - _cells[i]->T());
    }

    res[3*i+2] = res_energy / _rhoh_ref;
  }
}

void
OneDFlowChannel::onTimestepBegin()
{
}

void
OneDFlowChannel::onTimestepEnd()
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
OneDFlowChannel::writeTextOutput()
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
OneDFlowChannel::FillJacobianMatrixNonZeroEntry(MatrixNonZeroPattern * mnzp)
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

  if (_hs)
  {
    std::vector<std::vector<unsigned>> & Ts_DOFs = _hs->Ts_DOFs();
    for (unsigned j = 0; j < _cells.size(); j++)
    {
      unsigned Ts_DOF = (_hs_side == 0) ? Ts_DOFs[j].front() : Ts_DOFs[j].back();
      mnzp->addEntry(_cells[j]->TDOF(), Ts_DOF);
    }
  }
}
