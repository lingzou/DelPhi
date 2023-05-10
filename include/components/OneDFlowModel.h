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

  virtual void initialize(Real p, Real T);
  virtual void updateSolution(Real p, Real T);
  virtual void saveOldSlns();

  virtual void setWestEdge(EdgeBase * w_edge) { _w_edge = w_edge; }
  virtual void setEastEdge(EdgeBase * e_edge) { _e_edge = e_edge; }
  virtual EdgeBase * wEdge() const { return _w_edge; }
  virtual EdgeBase * eEdge() const { return _e_edge; }
  virtual EdgeBase * getOtherSideEdge(EdgeBase * edge);
  virtual void setExtendedNeighborCells();

  virtual void addConnectedDOFs(unsigned dof) { _connected_DOFs.insert(dof); }
  virtual std::set<unsigned> & getConnectedDOFs() { return _connected_DOFs; }

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

  EdgeBase * _w_edge;
  EdgeBase * _e_edge;
  CellBase * _w_cell;
  CellBase * _e_cell;

  std::set<unsigned> _connected_DOFs;
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
  //IntEdge(CellBase * w_cell, CellBase * e_cell) : _w_cell(w_cell), _e_cell(e_cell)
  virtual ~IntEdge() {}

  virtual Real mass_flux() override { return (_v > 0.0) ? _v * _w_cell->rho() : _v * _e_cell->rho(); }
  virtual Real enthalpy_flux() override { return (_v > 0.0) ? _v * _w_cell->rhoh() : _v * _e_cell->rhoh(); }
  virtual Real dv_dt(Real dt) override { return (_v - _v_o) / dt; }
  virtual Real dv_dx() override;
  virtual Real dp_dx() override { return (_e_cell->p() - _w_cell->p()) / _dL_edge; }

  // TODO: update to volume based average
  virtual Real rho_edge() final { return 0.5 * (_w_cell->rho() + _e_cell->rho()); }
};

class vBCEdge : public EdgeBase
{
public:
  vBCEdge(const InputParameters & parameters);
  virtual ~vBCEdge() {}

  virtual Real mass_flux() override { return (_v_bc > 0.0) ? _v_bc * _eos->rho_from_p_T(_e_cell->p(), _T_bc) : _v_bc * _e_cell->rho(); }
  virtual Real enthalpy_flux() override { return (_v_bc > 0.0) ? _v_bc * _eos->rho_from_p_T(_e_cell->p(), _T_bc) * _eos->h_from_p_T(_e_cell->p(), _T_bc) : _v_bc * _e_cell->rhoh(); }
  virtual Real dv_dt(Real /*dt*/) override { return 0.0; }
  virtual Real dv_dx() override;
  virtual Real dp_dx() override { return (_e_cell->p() - _p_ghost) / _dL_edge; }

  virtual Real rho_edge() final { return _e_cell->rho(); }

  virtual void updateGhostPressure(Real p_ghost) { _p_ghost = p_ghost; }
  virtual Real computeDirichletBCResidual() { return _v - _v_bc; }

protected:
  Real _v_bc, _T_bc;
  Real _p_ghost;
};

class pBCEdge : public EdgeBase
{
public:
  pBCEdge(const InputParameters & parameters);
  virtual ~pBCEdge() {}

  virtual Real mass_flux() override { return (_v > 0.0) ? _v * _w_cell->rho() : _v * _eos->rho_from_p_T(_p_bc, _T_bc); }
  virtual Real enthalpy_flux() override { return (_v > 0.0) ? _v * _w_cell->rhoh() : _v * _eos->rho_from_p_T(_p_bc, _T_bc) * _eos->h_from_p_T(_p_bc, _T_bc); }
  virtual Real dv_dt(Real dt) override { return (_v - _v_o) / dt; }
  virtual Real dv_dx() override;
  virtual Real dp_dx() override { return (_p_bc - _w_cell->p()) / _dL_edge; }

  virtual Real rho_edge() override final { return _w_cell->rho(); }

protected:
  Real _p_bc, _T_bc;
};
