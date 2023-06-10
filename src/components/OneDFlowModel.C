#include "OneDFlowModel.h"

CellBase::CellBase(const InputParameters & parameters) :
  _pars(parameters),
  _sim(*_pars.get<DelPhiSimulation *>("_sim")),
  _name(_pars.get<std::string>("name")),
  _dL_cell(_pars.get<Real>("dL")),
  _eos(_pars.get<const SinglePhaseFluidProperties *>("eos"))
{
}

void
CellBase::initialize(Real p, Real T)
{
  _p = p; _p_o = p; _p_oo = p;
  _T = T; _T_o = T; _T_oo = T;

  _rho = _eos->rho_from_p_T(p, T);
  _rho_o = _rho;
  _rho_oo = _rho;

  _h = _eos->h_from_p_T(p, T);
  _h_o = _h;
  _h_oo = _h;
}

void
CellBase::updateSolution(Real p, Real T)
{
  _p = p;
  _T = T;
  _p_w = _p_e = p;
  _T_w = _T_e = T;

  _rho = _eos->rho_from_p_T(p, T);
  _h = _eos->h_from_p_T(p, T);
  _rho_w = _rho_e = _rho;
  _h_w = _h_e = _h;
}

void
CellBase::saveOldSlns()
{
  _p_oo = _p_o;
  _T_oo = _T_o;
  _rho_oo = _rho_o;
  _h_oo = _h_o;

  _p_o = _p;
  _T_o = _T;
  _rho_o = _rho;
  _h_o = _h;
}

OneDCell::OneDCell(const InputParameters & parameters) :
  CellBase(parameters),
  _w_edge(NULL),
  _e_edge(NULL),
  _w_cell(NULL),
  _e_cell(NULL),
  _dh(_pars.get<Real>("Dh")),
  _wf_option(_pars.get<DELPHI::WallFrictionModel>("WF_option")),
  _f_const(_pars.isParamValid("f") ? _pars.get<Real>("f") : 0.0),
  _Welander_constant(_pars.isParamValid("Welander_constant") ? _pars.get<Real>("Welander_constant") : 0.0)
{
}

void
OneDCell::computeHelperVariables()
{
  Real v_avg = 0.5 * (_w_edge->v() + _e_edge->v());

  if (_wf_option == DELPHI::CONST_FRICTION)
    _dp_friction = 0.5 * _f_const * _dL_cell / _dh * _rho * v_avg * std::fabs(v_avg);
  else if (_wf_option == DELPHI::WELANDER)
    _dp_friction = _dL_cell * _Welander_constant * v_avg; // rho_0 * R * v
  else
    mooseError("Wall friction model does not exist.");
}

void
OneDCell::linearReconstruction(Real p_W, Real T_W, Real p_E, Real T_E)
{
  // Reference:
  // [1] Berger, M., Aftosmis, M.J., Murman, S.M., 2005. Analysis of slope limiters on irreg- ular grids.
  //     In: AIAA Paper 2005-0490, 43rd AIAA Aerospace Sciences Meeting, Jan. 10e13, Reno, NV, 2005.

  // pressure reconstruction
  if ((p_E - _p)*(_p - p_W) > 0.0)      // _p is in between p_W and p_E, i.e., not a local extremum
    if (std::fabs(p_E - p_W) > 1.e-10)  // The difference is large enough, so reconstruction makes sense, and that no divided-by-zero
    {
      Real f = (_p - p_W) / (p_E - p_W);
      Real wf = 2.0 * f * (1.0 - f) / ((1.0 - f) * (1.0 - f) + f * f);   // Eqn. (10) of Ref. [1]

      _p_w = _p - 0.25 * wf * (p_E - p_W);   // Eqn. (5) of Ref. [1]
      _p_e = _p + 0.25 * wf * (p_E - p_W);
    }

  // repeat for temperature
  if ((T_E - _T)*(_T - T_W) > 0.0)
    if (std::fabs(T_E - T_W) > 1.e-10)
    {
      Real f = (_T - T_W) / (T_E - T_W);
      Real wf = 2.0 * f * (1.0 - f) / ((1.0 - f) * (1.0 - f) + f * f);

      _T_w = _T - 0.25 * wf * (T_E - T_W);
      _T_e = _T + 0.25 * wf * (T_E - T_W);
    }

  _rho_w = _eos->rho_from_p_T(_p_w, _T_w);
  _rho_e = _eos->rho_from_p_T(_p_e, _T_e);
  _h_w = _eos->h_from_p_T(_p_w, _T_w);
  _h_e = _eos->h_from_p_T(_p_e, _T_e);
}

EdgeBase *
OneDCell::getOtherSideEdge(EdgeBase * edge)
{
  if (edge == _w_edge)        return _e_edge;
  else if (edge == _e_edge)   return _w_edge;
  else  mooseError("Cell '" + name() + "' is not connected to edge '" + edge->name() + "'.");
}

void
OneDCell::setExtendedNeighborCells()
{
  _w_cell = _w_edge->getOtherSideCell(this);
  _e_cell = _e_edge->getOtherSideCell(this);

  /*debug
  std::cerr << ((_w_cell == NULL) ? "NULL" : _w_cell->name()) << " - "
            << _w_edge->name() << " - "
            << name() << " - "
            << _e_edge->name() << " - "
            << ((_e_cell == NULL) ? "NULL" : _e_cell->name()) << std::endl;*/

  addConnectedDOFs(_pDOF);
  addConnectedDOFs(_TDOF);
  if (_w_cell)
  {
    addConnectedDOFs(_w_cell->pDOF());
    addConnectedDOFs(_w_cell->TDOF());
  }
  if (_e_cell)
  {
    addConnectedDOFs(_e_cell->pDOF());
    addConnectedDOFs(_e_cell->TDOF());
  }
  addConnectedDOFs(_w_edge->vDOF());
  addConnectedDOFs(_e_edge->vDOF());
}

BranchCell::BranchCell(const InputParameters & parameters) :
  CellBase(parameters)
{
}

void
BranchCell::addEdge(EdgeBase* edge, Real edge_out_norm)
{
  _connected_edges.push_back(edge);
  _edge_out_norms.push_back(edge_out_norm);
}

void
BranchCell::setExtendedNeighborCells()
{
  for (unsigned i = 0; i < _connected_edges.size(); i++)
    _connected_cells.push_back(_connected_edges[i]->getOtherSideCell(this));

  /* debug
  std::cerr << name() << std::endl;
  for (unsigned i = 0; i < _connected_edges.size(); i++)
  {
    std::cerr << "  edge = " << _connected_edges[i]->name() << std::endl;
    std::cerr << "  cell = " << _connected_cells[i]->name() << std::endl;
  }*/

  addConnectedDOFs(_pDOF);
  addConnectedDOFs(_TDOF);
  for (auto & edge : _connected_edges)
    if (edge)
      addConnectedDOFs(edge->vDOF());

  for (auto & cell : _connected_cells)
    if (cell)
    {
      addConnectedDOFs(cell->pDOF());
      addConnectedDOFs(cell->TDOF());
    }
}

EdgeBase::EdgeBase(const InputParameters & parameters) :
  _pars(parameters),
  _sim(*_pars.get<DelPhiSimulation *>("_sim")),
  _name(_pars.get<std::string>("name")),
  _dL_edge(0),
  _eos(_pars.get<const SinglePhaseFluidProperties *>("eos")),
  _w_cell(_pars.get<CellBase*>("west_cell")),
  _e_cell(_pars.get<CellBase*>("east_cell")),
  _w_edge(NULL),
  _e_edge(NULL)
{
}

CellBase *
EdgeBase::getOtherSideCell(CellBase * cell)
{
  if (cell == _w_cell)        return _e_cell;
  else if (cell == _e_cell)   return _w_cell;
  else  mooseError("Edge '" + name() + "' is not connected to cell '" + cell->name() + "'.");
}

void
EdgeBase::setExtendedNeighborEdges()
{
  if (_w_cell) _w_edge = _w_cell->getOtherSideEdge(this);
  if (_e_cell) _e_edge = _e_cell->getOtherSideEdge(this);

  addConnectedDOFs(_vDOF);
  if (_w_cell)
  {
    addConnectedDOFs(_w_cell->pDOF());
    addConnectedDOFs(_w_cell->TDOF());
  }
  if (_e_cell)
  {
    addConnectedDOFs(_e_cell->pDOF());
    addConnectedDOFs(_e_cell->TDOF());
  }
  if (_w_edge)
    addConnectedDOFs(_w_edge->vDOF());
  if (_e_edge)
    addConnectedDOFs(_e_edge->vDOF());
}


IntEdge::IntEdge(const InputParameters & parameters) :
  EdgeBase(parameters)
{
  // an IntEdge connects two OneDCell, so make sure this is happening
  if (_w_cell && dynamic_cast<OneDCell*>(_w_cell))
  {
    (dynamic_cast<OneDCell*>(_w_cell))->setEastEdge(this);
    _dL_edge += 0.5 * _w_cell->dL();
  }
  else
    mooseError("IntEdge: '" + name() + "' missing west_cell or west_cell is not a OneDCell.");

  if (_e_cell && dynamic_cast<OneDCell*>(_e_cell))
  {
    (dynamic_cast<OneDCell*>(_e_cell))->setWestEdge(this);
    _dL_edge += 0.5 * _e_cell->dL();
  }
  else
    mooseError("IntEdge: '" + name() + "' missing east_cell or east_cell is not a OneDCell.");
}

Real
IntEdge::dv_dx()
{
  if (_v > 0.0)
    return (_v - _w_edge->v()) / _w_cell->dL();
  else
    return (_e_edge->v() - _v) / _e_cell->dL();
}

vBCEdge::vBCEdge(const InputParameters & parameters) :
  boundaryEdge(parameters),
  _v_bc(_sim.getFunction(_pars.get<FunctionName>("v_bc"))),
  _T_bc(_sim.getFunction(_pars.get<FunctionName>("T_bc")))
{
}

vBCEdgeInlet::vBCEdgeInlet(const InputParameters & parameters) :
  vBCEdge(parameters)
{
  // expect no w_cell but e_cell which must be OneDCell
  if (_w_cell) mooseError("Not expecting west_cell in pBCEdgeInlet");

  if (_e_cell && dynamic_cast<OneDCell*>(_e_cell))
  {
    (dynamic_cast<OneDCell*>(_e_cell))->setWestEdge(this);
    _dL_edge += 0.5 * _e_cell->dL();
  }
  else
    mooseError("vBCEdgeInlet: '" + name() + "' missing east_cell or east_cell is not a OneDCell.");
}

vBCEdgeOutlet::vBCEdgeOutlet(const InputParameters & parameters) :
  vBCEdge(parameters)
{
  // expect no e_cell but w_cell which must be OneDCell
  if (_w_cell && dynamic_cast<OneDCell*>(_w_cell))
  {
    (dynamic_cast<OneDCell*>(_w_cell))->setEastEdge(this);
    _dL_edge += 0.5 * _w_cell->dL();
  }
  else
    mooseError("vBCEdgeOutlet: '" + name() + "' missing west_cell or west_cell is not a OneDCell.");

  if (_e_cell) mooseError("Not expecting east_cell in pBCEdgeOutlet");
}

pBCEdge::pBCEdge(const InputParameters & parameters) :
  boundaryEdge(parameters),
  _p_bc(_pars.get<Real>("p_bc")),
  _T_bc(_pars.get<Real>("T_bc"))
{
}

pBCEdgeInlet::pBCEdgeInlet(const InputParameters & parameters) :
  pBCEdge(parameters)
{
  // expect no w_cell but e_cell which must be OneDCell
  if (_w_cell) mooseError("Not expecting west_cell in pBCEdgeInlet");

  if (_e_cell && dynamic_cast<OneDCell*>(_e_cell))
  {
    (dynamic_cast<OneDCell*>(_e_cell))->setWestEdge(this);
    _dL_edge += 0.5 * _e_cell->dL();
  }
  else
    mooseError("pBCEdgeInlet: '" + name() + "' missing east_cell or east_cell is not a OneDCell.");
}

pBCEdgeOutlet::pBCEdgeOutlet(const InputParameters & parameters) :
  pBCEdge(parameters)
{
  // expect no e_cell but w_cell which must be OneDCell
  if (_w_cell && dynamic_cast<OneDCell*>(_w_cell))
  {
    (dynamic_cast<OneDCell*>(_w_cell))->setEastEdge(this);
    _dL_edge += 0.5 * _w_cell->dL();
  }
  else
    mooseError("pBCEdgeOutlet: '" + name() + "' missing west_cell or west_cell is not a OneDCell.");

  if (_e_cell) mooseError("Not expecting east_cell in pBCEdgeOutlet");
}

brvEdgeInlet::brvEdgeInlet(const InputParameters & parameters) :
  EdgeBase(parameters)
{
  // expect w_cell to be BranchCell
  if (_w_cell && dynamic_cast<BranchCell*>(_w_cell))
  {
    (dynamic_cast<BranchCell*>(_w_cell))->addEdge(this, /*out_norm*/ 1.0);
    _dL_edge += 0.5 * _w_cell->dL();
  }
  else
    mooseError("brvEdgeInlet: '" + name() + "' missing west_cell or west_cell is not a BranchCell.");

  // expect e_cell to be OneDCell
  if (_e_cell && dynamic_cast<OneDCell*>(_e_cell))
  {
    (dynamic_cast<OneDCell*>(_e_cell))->setWestEdge(this);
    _dL_edge += 0.5 * _e_cell->dL();
  }
  else
    mooseError("brvEdgeInlet: '" + name() + "' missing east_cell or east_cell is not a OneDCell.");
}

brvEdgeOutlet::brvEdgeOutlet(const InputParameters & parameters) :
  EdgeBase(parameters)
{
  // expect w_cell to be OneDCell
  if (_w_cell && dynamic_cast<OneDCell*>(_w_cell))
  {
    (dynamic_cast<OneDCell*>(_w_cell))->setEastEdge(this);
    _dL_edge += 0.5 * _w_cell->dL();
  }
  else
    mooseError("brvEdgeOutlet: '" + name() + "' missing west_cell or west_cell is not a BranchCell.");

  // expect e_cell to be BranchCell
  if (_e_cell && dynamic_cast<BranchCell*>(_e_cell))
  {
    (dynamic_cast<BranchCell*>(_e_cell))->addEdge(this, /*out_norm*/ -1.0);
    _dL_edge += 0.5 * _e_cell->dL();
  }
  else
    mooseError("brvEdgeOutlet: '" + name() + "' missing east_cell or east_cell is not a OneDCell.");
}
