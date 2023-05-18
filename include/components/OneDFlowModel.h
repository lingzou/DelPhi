#pragma once

#include "SinglePhaseFluidProperties.h"

class CellBase;
class EdgeBase;

class CellBase
{
public:
  CellBase(const InputParameters & parameters);
  virtual ~CellBase() {}

  virtual std::string name() final { return _name; }
  virtual Real dL() final { return _dL_cell; }

  virtual Real p()    const final { return _p;    }
  virtual Real T()    const final { return _T;    }
  virtual Real rho()  const final { return _rho;  }
  virtual Real h()    const final { return _h;    }
  virtual Real rhoh()   const final { return _rho * _h; }
  virtual Real rho_o()  const final { return _rho_o;  }
  virtual Real h_o()    const final { return _h_o;    }
  virtual Real rhoh_o() const final { return _rho_o * _h_o; }

  // Residual related functions
  virtual void setDOF(unsigned pDOF, unsigned TDOF) final { _pDOF = pDOF; _TDOF = TDOF; }
  virtual unsigned pDOF() const final { return _pDOF; }
  virtual unsigned TDOF() const final { return _TDOF; }

  virtual void initialize(Real p, Real T) final;
  virtual void updateSolution(Real p, Real T) final;
  virtual void saveOldSlns() final;

  virtual EdgeBase * getOtherSideEdge(EdgeBase * edge) = 0;
  virtual void setExtendedNeighborCells() = 0;
  virtual void addConnectedDOFs(unsigned dof) final { _connected_DOFs.insert(dof); }
  virtual std::set<unsigned> & getConnectedDOFs() final { return _connected_DOFs; }

protected:
  const InputParameters & _pars;
  std::string _name;
  Real _dL_cell;

  Real _p, _p_o;
  Real _T, _T_o;
  Real _rho, _rho_o;
  Real _h, _h_o;

  unsigned _pDOF, _TDOF;

  const SinglePhaseFluidProperties * _eos;

  std::set<unsigned> _connected_DOFs;
};

class OneDCell : public CellBase
{
public:
  OneDCell(const InputParameters & parameters);
  virtual ~OneDCell() {}

  virtual void setWestEdge(EdgeBase * w_edge) { _w_edge = w_edge; }
  virtual void setEastEdge(EdgeBase * e_edge) { _e_edge = e_edge; }
  virtual EdgeBase * wEdge() const { return _w_edge; }
  virtual EdgeBase * eEdge() const { return _e_edge; }
  virtual EdgeBase * getOtherSideEdge(EdgeBase * edge) override;
  virtual void setExtendedNeighborCells() override;

protected:
  EdgeBase * _w_edge;
  EdgeBase * _e_edge;
  CellBase * _w_cell;
  CellBase * _e_cell;
};

class BranchCell : public CellBase
{
public:
  BranchCell(const InputParameters & parameters);
  virtual ~BranchCell() {}

  virtual void addEdge(EdgeBase* edge, Real edge_out_norm);

  // there is no 'other' side because there is unknown (could be many) connected edges
  virtual EdgeBase * getOtherSideEdge(EdgeBase * /*edge*/) override { return NULL; }
  virtual void setExtendedNeighborCells() override;

protected:
  std::vector<CellBase*> _connected_cells;
  std::vector<EdgeBase*> _connected_edges;
  std::vector<Real> _edge_out_norms;
};

class EdgeBase
{
public:
  EdgeBase(const InputParameters & parameters);
  virtual ~EdgeBase() {}

  virtual std::string name() final { return _name; }
  virtual Real dL() final { return _dL_edge; }

  virtual Real v()           const final { return _v; }
  virtual Real v_o()         const final { return _v_o; }
  virtual Real mass_flux() = 0;
  virtual Real enthalpy_flux() = 0;
  virtual Real rho_edge() = 0;
  virtual Real dv_dt(Real dt) = 0;
  virtual Real dv_dx() = 0;
  virtual Real dp_dx() = 0;
  // only for output purpose
  virtual Real T_edge() = 0;

  virtual void setDOF(unsigned vDOF) final { _vDOF = vDOF; }
  virtual unsigned vDOF() const final { return _vDOF; }
  virtual void initialize(Real v)     final { _v = v; _v_o = v; }
  virtual void updateSolution(Real v) final { _v = v; }
  virtual void saveOldSlns()          final { _v_o = _v; }

  virtual CellBase * getOtherSideCell(CellBase * cell);
  virtual void setExtendedNeighborEdges();

  virtual void addConnectedDOFs(unsigned dof) { _connected_DOFs.insert(dof); }
  virtual std::set<unsigned> & getConnectedDOFs() { return _connected_DOFs; }

protected:
  const InputParameters & _pars;
  std::string _name;
  Real _dL_edge;

  Real _v, _v_o;
  unsigned _vDOF;

  const SinglePhaseFluidProperties * _eos;

  CellBase * _w_cell;
  CellBase * _e_cell;
  EdgeBase * _w_edge;
  EdgeBase * _e_edge;

  std::set<unsigned> _connected_DOFs;
};

class IntEdge : public EdgeBase
{
public:
  IntEdge(const InputParameters & parameters);
  virtual ~IntEdge() {}

  virtual Real T_edge() override { return (_v > 0.0) ? _w_cell->T() : _e_cell->T(); }
  virtual Real mass_flux() override { return (_v > 0.0) ? _v * _w_cell->rho() : _v * _e_cell->rho(); }
  virtual Real enthalpy_flux() override { return (_v > 0.0) ? _v * _w_cell->rhoh() : _v * _e_cell->rhoh(); }
  virtual Real dv_dt(Real dt) override { return (_v - _v_o) / dt; }
  virtual Real dv_dx() override;
  virtual Real dp_dx() override { return (_e_cell->p() - _w_cell->p()) / _dL_edge; }

  // TODO: update to volume based average
  virtual Real rho_edge() final { return 0.5 * (_w_cell->rho() + _e_cell->rho()); }
};

// snjEdge is a type of IntEdge for connecting the ends of two pipes
// not like IntEdge is inside a single pipe
// it will need special treatment such as area change, irregular connection such as inlet-inlet situations
class snjEdge : public IntEdge
{
public:
  snjEdge(const InputParameters & parameters) : IntEdge(parameters) {}
  virtual ~snjEdge() {}
};

class vBCEdge : public EdgeBase
{
public:
  vBCEdge(const InputParameters & parameters);
  virtual ~vBCEdge() {}

  virtual Real dv_dt(Real /*dt*/) override { return 0.0; }

  virtual void updateGhostPressure(Real p_ghost) { _p_ghost = p_ghost; }
  virtual Real computeDirichletBCResidual() { return _v - _v_bc; }

protected:
  Real _v_bc, _T_bc;
  Real _p_ghost;
};

class vBCEdgeInlet : public vBCEdge
{
public:
  vBCEdgeInlet(const InputParameters & parameters);
  virtual ~vBCEdgeInlet() {}

  virtual Real T_edge() override { return (_v_bc > 0.0) ? _T_bc : _e_cell->T(); }
  virtual Real mass_flux() override { return (_v_bc > 0.0) ? _v_bc * _eos->rho_from_p_T(_p_ghost, _T_bc) : _v_bc * _e_cell->rho(); }
  virtual Real enthalpy_flux() override { return (_v_bc > 0.0) ? _v_bc * _eos->rho_from_p_T(_p_ghost, _T_bc) * _eos->h_from_p_T(_p_ghost, _T_bc) : _v_bc * _e_cell->rhoh(); }

  virtual Real dv_dx() override { return (_v_bc > 0.0) ? 0.0 : (_e_edge->v() - _v) / _e_cell->dL(); }
  virtual Real dp_dx() override { return (_e_cell->p() - _p_ghost) / _dL_edge; }

  virtual Real rho_edge() final { return _e_cell->rho(); }
};

class vBCEdgeOutlet : public vBCEdge
{
public:
  vBCEdgeOutlet(const InputParameters & parameters);
  virtual ~vBCEdgeOutlet() {}

  virtual Real T_edge() override { return (_v_bc > 0.0) ? _w_cell->T() : _T_bc; }
  virtual Real mass_flux() override { return (_v_bc > 0.0) ? _v_bc * _w_cell->rho() : _v_bc * _eos->rho_from_p_T(_p_ghost, _T_bc); }
  virtual Real enthalpy_flux() override { return (_v_bc > 0.0) ? _v_bc * _w_cell->rhoh() : _v_bc * _eos->rho_from_p_T(_p_ghost, _T_bc) * _eos->h_from_p_T(_p_ghost, _T_bc); }
  //virtual Real dv_dt(Real /*dt*/) override { return 0.0; }
  virtual Real dv_dx() override { return (_v_bc > 0.0) ? (_v_bc - _w_edge->v() / _w_cell->dL()) : 0.0; }
  virtual Real dp_dx() override { return (_p_ghost - _w_cell->p()) / _dL_edge; }

  virtual Real rho_edge() final { return _w_cell->rho(); }
};



class snjShadowEdge : public vBCEdge
{
public:
  snjShadowEdge(const InputParameters & parameters) : vBCEdge(parameters),
  _real_edge(_pars.get<snjEdge*>("real_edge"))
  {}
  virtual ~snjShadowEdge() {}

  virtual Real T_edge() override { return _real_edge->T_edge(); }
  virtual Real mass_flux() override { return _real_edge->mass_flux(); }
  virtual Real enthalpy_flux() override { return _real_edge->enthalpy_flux(); }
  virtual Real dv_dx() override { return _real_edge->dv_dx(); };
  virtual Real dp_dx() override { return _real_edge->dp_dx(); }
  virtual Real rho_edge() final { return _real_edge->rho_edge(); }

  virtual void updateGhostPressure(Real /*p_ghost*/) override { /* nothing to do */ }

  // Key implementation: show edge has the same velocity as the real edge
  virtual Real computeDirichletBCResidual() override { return _v - _real_edge->v(); }

  virtual void setExtendedNeighborEdges() override
  {
    // shadow edge has no connected cells but follow the real edge
    _connected_DOFs.insert(_vDOF);
    _connected_DOFs.insert(_real_edge->vDOF());
  }

protected:
  snjEdge * _real_edge;
};


class pBCEdge : public EdgeBase
{
public:
  pBCEdge(const InputParameters & parameters);
  virtual ~pBCEdge() {}

  virtual Real dv_dt(Real dt) override final { return (_v - _v_o) / dt; }

protected:
  Real _p_bc, _T_bc;
};

class pBCEdgeInlet : public pBCEdge
{
public:
  pBCEdgeInlet(const InputParameters & parameters);
  virtual ~pBCEdgeInlet() {}

  virtual Real T_edge() override { return (_v > 0.0) ? _T_bc : _e_cell->T(); }
  virtual Real mass_flux() override { return (_v > 0.0) ? _v * _eos->rho_from_p_T(_p_bc, _T_bc) : _v * _e_cell->rho(); }
  virtual Real enthalpy_flux() override { return (_v > 0.0) ? _v * _eos->rho_from_p_T(_p_bc, _T_bc) * _eos->h_from_p_T(_p_bc, _T_bc) : _v * _e_cell->rhoh(); }
  virtual Real dv_dx() override { return (_v > 0.0) ? 0.0 : (_e_edge->v() - _v) / _e_cell->dL(); }
  virtual Real dp_dx() override { return (_e_cell->p() - _p_bc) / _dL_edge; }

  virtual Real rho_edge() override final { return _e_cell->rho(); }
};

class pBCEdgeOutlet : public pBCEdge
{
public:
  pBCEdgeOutlet(const InputParameters & parameters);
  virtual ~pBCEdgeOutlet() {}

  virtual Real T_edge() override { return (_v > 0.0) ? _w_cell->T() : _T_bc; }
  virtual Real mass_flux() override { return (_v > 0.0) ? _v * _w_cell->rho() : _v * _eos->rho_from_p_T(_p_bc, _T_bc); }
  virtual Real enthalpy_flux() override { return (_v > 0.0) ? _v * _w_cell->rhoh() : _v * _eos->rho_from_p_T(_p_bc, _T_bc) * _eos->h_from_p_T(_p_bc, _T_bc); }
  virtual Real dv_dx() override { return (_v > 0.0) ? (_v - _w_edge->v()) / _w_cell->dL() : 0.0; }
  virtual Real dp_dx() override { return (_p_bc - _w_cell->p()) / _dL_edge; }

  virtual Real rho_edge() override final { return _w_cell->rho(); }
};

class brvEdgeInlet : public EdgeBase
{
public:
  brvEdgeInlet(const InputParameters & parameters);
  virtual ~brvEdgeInlet() {}

  virtual Real T_edge() override { return (_v > 0.0) ? _w_cell->T() : _e_cell->T(); }
  virtual Real mass_flux() override { return (_v > 0.0) ? _v * _w_cell->rho() : _v * _e_cell->rho(); }
  virtual Real enthalpy_flux() override { return (_v > 0.0) ? _v * _w_cell->rhoh() : _v * _e_cell->rhoh(); }
  virtual Real dv_dt(Real dt) override { return (_v - _v_o) / dt; }
  virtual Real dv_dx() override { return (_v > 0.0) ? 0.0 : (_e_edge->v() - _v) / _e_cell->dL(); }
  virtual Real dp_dx() override { return (_e_cell->p() - _w_cell->p()) / _dL_edge; }

  // TODO: update to volume based average
  virtual Real rho_edge() final { return 0.5 * (_w_cell->rho() + _e_cell->rho()); }
};

class brvEdgeOutlet : public EdgeBase
{
public:
  brvEdgeOutlet(const InputParameters & parameters);
  virtual ~brvEdgeOutlet() {}

  virtual Real T_edge() override { return (_v > 0.0) ? _w_cell->T() : _e_cell->T(); }
  virtual Real mass_flux() override { return (_v > 0.0) ? _v * _w_cell->rho() : _v * _e_cell->rho(); }
  virtual Real enthalpy_flux() override { return (_v > 0.0) ? _v * _w_cell->rhoh() : _v * _e_cell->rhoh(); }
  virtual Real dv_dt(Real dt) override { return (_v - _v_o) / dt; }
  virtual Real dv_dx() override { return (_v > 0.0) ? (_v - _w_edge->v()) / _w_cell->dL() : 0.0; }
  virtual Real dp_dx() override { return (_e_cell->p() - _w_cell->p()) / _dL_edge; }

  // TODO: update to volume based average
  virtual Real rho_edge() final { return 0.5 * (_w_cell->rho() + _e_cell->rho()); }
};
