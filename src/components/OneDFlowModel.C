#include "OneDFlowModel.h"

CellBase::CellBase(const InputParameters & parameters) :
  _pars(parameters),
  _name(_pars.get<std::string>("name")),
  _dL_cell(_pars.get<Real>("dL")),
  _eos(_pars.get<const SinglePhaseFluidProperties *>("eos")),
  _w_edge(NULL),
  _e_edge(NULL),
  _w_cell(NULL),
  _e_cell(NULL)
{
}

EdgeBase *
CellBase::getOtherSideEdge(EdgeBase * edge)
{
  if (edge == _w_edge)        return _e_edge;
  else if (edge == _e_edge)   return _w_edge;
  else  mooseError("Cell '" + name() + "' is not connected to edge '" + edge->name() + "'.");
}

void
CellBase::setExtendedNeighborCells()
{
  _w_cell = _w_edge->getOtherSideCell(this);
  _e_cell = _e_edge->getOtherSideCell(this);

  /* debug
  std::cerr << ((_w_cell == NULL) ? "NULL" : _w_cell->name()) << " - "
            << _w_edge->name() << " - "
            << name() << " - "
            << _e_edge->name() << " - "
            << ((_e_cell == NULL) ? "NULL" : _e_cell->name()) << std::endl;
  */

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

void
CellBase::initialize(Real p, Real T)
{
  _p = p; _p_o = p;
  _T = T; _T_o = T;
  _rho = _eos->rho_from_p_T(p, T);
  _rho_o = _rho;
  _h = _eos->h_from_p_T(p, T);
  _h_o = _h;
}

void
CellBase::updateSolution(Real p, Real T)
{
  _p = p;
  _T = T;
  _rho = _eos->rho_from_p_T(p, T);
  _h = _eos->h_from_p_T(p, T);
}

void
CellBase::saveOldSlns()
{
  _p_o = _p;
  _T_o = _T;
  _rho_o = _rho;
  _h_o = _h;
}


EdgeBase::EdgeBase(const InputParameters & parameters) :
  _pars(parameters),
  _name(_pars.get<std::string>("name")),
  _dL_edge(0),
  _eos(_pars.get<const SinglePhaseFluidProperties *>("eos")),
  _w_cell(_pars.get<CellBase*>("west_cell")),
  _e_cell(_pars.get<CellBase*>("east_cell")),
  _w_edge(NULL),
  _e_edge(NULL)
{
  if (_w_cell)
  {
    _w_cell->setEastEdge(this);
    _dL_edge += 0.5 * _w_cell->dL();
  }
  if (_e_cell)
  {
    _e_cell->setWestEdge(this);
    _dL_edge += 0.5 * _e_cell->dL();
  }
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
  EdgeBase(parameters),
  _v_bc(_pars.get<Real>("v_bc")),
  _T_bc(_pars.get<Real>("T_bc"))
{
}

Real
vBCEdge::dv_dx()
{
  if (_v > 0.0)
    return 0.0;
  else
    return (_e_edge->v() - _v) / _e_cell->dL();
}

pBCEdge::pBCEdge(const InputParameters & parameters) :
  EdgeBase(parameters),
  _p_bc(_pars.get<Real>("p_bc")),
  _T_bc(_pars.get<Real>("T_bc"))
{
}

Real
pBCEdge::dv_dx()
{
  if (_v > 0.0)
    return (_v - _w_edge->v()) / _w_cell->dL();
  else
    return 0.0;
}
