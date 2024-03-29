#pragma once

// moose include
#include "SinglePhaseFluidProperties.h"
#include "Function.h"

#include "utils.h"

class CellBase;
class EdgeBase;

/**
 * Under the context of finite volume method (FVM), a cell (or control volume) is one of the most
 * fundamental discretization units where volume-averaged quantities are defined and
 * in general, conservation laws are observed, e.g., mass and energy conservations.
 *
 * As for stagerred-grid FVM in flow problem applications, a cell is where the scalar quantities are
 * defined, for example, pressure (p) and temperature (T) for single-phase flow, and dependent
 * variables such as density and enthalpy naturally follows. The vector variables such as velocity
 * instead are defined on edges (or surfaces), defined later.
 *
 * In this application, the cell can be a one-dimensional cell in a flow channel type of component,
 * which has length, direction, and expecting two edges on its two ends (we call them west and east
 * end).
 *
 *                   -------------------------
 *  west_edge (v)  ->|    OneDCell (p, T)    |-> east_edge (v)
 *                   -------------------------
 *
 * The cell can also be an arbitrary one that having no real shape, or we do not need and cannot
 * capture its shape, but rather more concerned on its averaged behavior, such as a well-mixing
 * large volume with multiple inlet/outlet connections.
 */

/**
 * CellBase is the base class for all cells (control volumes) where the most comman data
 * and methods are defined.
 */
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
  virtual Real rho_oo()  const final { return _rho_oo;  }
  virtual Real h_oo()    const final { return _h_oo;    }
  virtual Real rhoh_oo() const final { return _rho_oo * _h_oo; }

  virtual Real drho_dt() const final
  {
    if ((_sim.TS() == Moose::TI_IMPLICIT_EULER) || (_sim.timeStep() == 1)) // bdf1 or step 1 of bdf2
      return (_rho - _rho_o) / _sim.dt();
    else //bdf2 and step > 1
      return DELPHI::BDF2_ddt(_rho, _rho_o, _rho_oo, _sim.dt(), _sim.dtOld());
  }
  virtual Real drhoh_dt() const final
  {
    if ((_sim.TS() == Moose::TI_IMPLICIT_EULER) || (_sim.timeStep() == 1)) // bdf1 or step 1 of bdf2
      return (_rho * _h - _rho_o * _h_o) / _sim.dt();
    else
      return DELPHI::BDF2_ddt(_rho * _h, _rho_o * _h_o, _rho_oo * _h_oo, _sim.dt(), _sim.dtOld());
  }

  virtual Real dp_friction()
  { mooseError("dp_friction() is not expected"); return 0.0; }
  virtual Real dp_gravity()
  { mooseError("dp_gravity() is not expected"); return 0.0; }

  // high-order data access
  virtual Real T_w() const final { return _T_w; }
  virtual Real T_e() const final { return _T_e; }
  virtual Real rho_w() const final { return _rho_w; }
  virtual Real rho_e() const final { return _rho_e; }
  virtual Real h_w() const final { return _h_w; }
  virtual Real h_e() const final { return _h_e; }
  virtual Real rhoh_w() const final { return _rho_w * _h_w; }
  virtual Real rhoh_e() const final { return _rho_e * _h_e; }

  // Residual related functions
  virtual void setDOF(unsigned pDOF, unsigned TDOF) final { _pDOF = pDOF; _TDOF = TDOF; }
  virtual unsigned pDOF() const final { return _pDOF; }
  virtual unsigned TDOF() const final { return _TDOF; }

  virtual void initialize(Real p, Real T) final;
  virtual void updateSolution(Real p, Real T) final;
  virtual void saveOldSlns() final;

  virtual void computeHelperVariables() {}

  virtual EdgeBase * getOtherSideEdge(EdgeBase * edge) = 0;
  virtual void setExtendedNeighborCells() = 0;
  virtual void addConnectedDOFs(unsigned dof) final { _connected_DOFs.insert(dof); }
  virtual std::set<unsigned> & getConnectedDOFs() final { return _connected_DOFs; }

protected:
  const InputParameters & _pars;
  DelPhiSimulation & _sim;
  std::string _name;

  Real _dL_cell;

  // current value, old (o) value, old old (oo) value
  Real _p, _p_o, _p_oo;
  Real _T, _T_o, _T_oo;
  Real _rho, _rho_o, _rho_oo;
  Real _h, _h_o, _h_oo;
  // high-order data
  Real _p_e, _p_w;
  Real _T_e, _T_w;
  Real _rho_e, _rho_w;
  Real _h_e, _h_w;

  unsigned _pDOF, _TDOF;

  const SinglePhaseFluidProperties * _eos;

  std::set<unsigned> _connected_DOFs;
};

/**
 * OneDCell is a derived class from CellBase for the directional one-dimensional cell
 * for one-dimensional flow channels.
 * It is expected that it is connected with two edges, west and east edges.
 * However, the cells to its west and east sides are only conditional. For example, if the
 * OneDCell is on an inlet boundary, the west side cell does not exist.
 */
class OneDCell : public CellBase
{
public:
  OneDCell(const InputParameters & parameters);
  virtual ~OneDCell() {}

  virtual void setWestEdge(EdgeBase * w_edge) { _w_edge = w_edge; }
  virtual void setEastEdge(EdgeBase * e_edge) { _e_edge = e_edge; }
  virtual EdgeBase * wEdge() const { return _w_edge; }
  virtual EdgeBase * eEdge() const { return _e_edge; }
  virtual CellBase * wCell() const { return _w_cell; }
  virtual CellBase * eCell() const { return _e_cell; }
  virtual EdgeBase * getOtherSideEdge(EdgeBase * edge) override;
  virtual void setExtendedNeighborCells() override;

  virtual void computeHelperVariables() override;
  virtual void linearReconstruction(Real p_W, Real T_W, Real p_E, Real T_E);
  virtual void linearReconstructionIrregularMesh();

  virtual Real dp_friction() override { return _dp_friction; }
  virtual Real dp_gravity() override { return _dp_gravity; }

protected:
  EdgeBase * _w_edge;
  EdgeBase * _e_edge;
  CellBase * _w_cell;
  CellBase * _e_cell;

  Real _gL;

  Real _dh;
  DELPHI::WallFrictionModel _wf_option;
  Real _f_const;
  Real _Welander_constant;

  Real _dp_friction;
  Real _dp_gravity;
};

/**
 * BranchCell is a derived class from CellBase which represents an arbitrary volume
 * with possibly multiple connected flow channels.
 */
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

/**
 * EdgeBase is the base class for all edges where velocity is defined and the momentum
 * equation is solved.
 * Edges (or surfaces) always attach to cells (control volume) and compute
 * mass and energy flux such that conservations of mass and energy in the connected
 * cells can be computed.
 */
class EdgeBase
{
public:
  EdgeBase(const InputParameters & parameters);
  virtual ~EdgeBase() {}

  virtual std::string name() final { return _name; }
  virtual Real dL() final { return _dL_edge; }

  virtual Real v()           const final { return _v; }
  virtual Real v_o()         const final { return _v_o; }
  virtual Real v_oo()        const final { return _v_oo; }

  virtual void computeFluxes() = 0;
  virtual Real mass_flux() { return _mass_flux; }
  virtual Real enthalpy_flux() { return _enthalpy_flux; }
  virtual Real rho_edge() = 0;
  virtual Real dv_dt()
  {
    if ((_sim.TS() == Moose::TI_IMPLICIT_EULER) || (_sim.timeStep() == 1)) // bdf1 or step 1 of bdf2
      return (_v - _v_o) / _sim.dt();
    else //bdf2 and step > 1
      return DELPHI::BDF2_ddt(_v, _v_o, _v_oo, _sim.dt(), _sim.dtOld());
  }
  virtual Real dv_dx() = 0;
  virtual Real dp_dx() = 0;
  virtual Real dp_dx_friction() = 0;
  virtual Real gravity() = 0;
  // only for output purpose
  virtual Real T_edge() { return _T_edge; }

  virtual void setDOF(unsigned vDOF) final { _vDOF = vDOF; }
  virtual unsigned vDOF() const final { return _vDOF; }
  virtual void initialize(Real v)     final { _v = v; _v_o = v; _v_oo = v; }
  virtual void updateSolution(Real v) final { _v = v; }
  virtual void saveOldSlns()          final { _v_oo = _v_o; _v_o = _v; }
  virtual void applyDirichletBC(Real & /*res*/) { /*not all edges has DirichletBC*/ }

  virtual CellBase * getOtherSideCell(CellBase * cell);
  virtual void setExtendedNeighborEdges();

  virtual void addConnectedDOFs(unsigned dof) { _connected_DOFs.insert(dof); }
  virtual std::set<unsigned> & getConnectedDOFs() { return _connected_DOFs; }

protected:
  const InputParameters & _pars;
  DelPhiSimulation & _sim;
  std::string _name;

  Real _dL_edge;

  Real _v, _v_o, _v_oo;
  unsigned _vDOF;

  Real _T_edge; // the temperature used to compute enthalpy flux, also for output
  Real _mass_flux; // rho * u
  Real _enthalpy_flux; // rho * u * h

  const SinglePhaseFluidProperties * _eos;

  CellBase * _w_cell;
  CellBase * _e_cell;
  EdgeBase * _w_edge;
  EdgeBase * _e_edge;

  std::set<unsigned> _connected_DOFs;
};

/**
 * IntEdge is a derived class from EdgeBase which is in between two OneDCells.
 * For an IntEdge, the west cell and edge, east cell and edge are always expected.
 */
class IntEdge : public EdgeBase
{
public:
  IntEdge(const InputParameters & parameters);
  virtual ~IntEdge() {}

  virtual void computeFluxes() override
  {
    if (_v > 0.0)
    {
      _T_edge = _w_cell->T_e();
      _mass_flux = _v * _w_cell->rho_e();
      _enthalpy_flux = _v * _w_cell->rhoh_e();
    }
    else
    {
      _T_edge = _e_cell->T_w();
      _mass_flux = _v * _e_cell->rho_w();
      _enthalpy_flux = _v * _e_cell->rhoh_w();
    }
  }

  virtual Real dv_dx() override;
  virtual Real dp_dx() override { return (_e_cell->p() - _w_cell->p()) / _dL_edge; }
  virtual Real dp_dx_friction() override
  {
    return 0.5 * (_e_cell->dp_friction() + _w_cell->dp_friction()) / _dL_edge;
  }
  virtual Real gravity() override
  {
    return 0.5 * (_e_cell->dp_gravity() + _w_cell->dp_gravity()) / _dL_edge;
  }

  // TODO: update to volume based average
  virtual Real rho_edge() override { return 0.5 * (_w_cell->rho() + _e_cell->rho()); }
};

/**
 * snjEdge is a special type of IntEdge for connecting the ends of two pipes, not like IntEdge
 * inside a single pipe. For now, we limit the use of snjEdge to connect the outlet of a pipe and
 * inlet of the other pipe, while the order does not matter. For now, we don't allow for connection
 * like pipe-1(out) - pipe-2(out).
 */
class snjEdge : public IntEdge
{
public:
  snjEdge(const InputParameters & parameters) : IntEdge(parameters) {}
  virtual ~snjEdge() {}
};

class boundaryEdge : public EdgeBase
{
public:
  boundaryEdge(const InputParameters & parameters) : EdgeBase(parameters) {}
  virtual ~boundaryEdge() {}

  virtual Real p_bc() = 0; // ghost p for velocity BC, real p for p BC
  virtual Real T_bc() = 0; // real T BC
};

/**
 * vBCEdge derives from EdgeBase, and serves as the base class for velocity boundary type of
 * implementation. For velocity boundary, v_bc and T_bc are required, while pressure is a projected
 * value from interior domains.
 */
class vBCEdge : public boundaryEdge
{
public:
  vBCEdge(const InputParameters & parameters);
  virtual ~vBCEdge() {}

  virtual Real p_bc() override { return _p_ghost; }
  virtual Real T_bc() override { return _T_bc.value(_sim.time(), Point()); }
  virtual Real dv_dt() override { return 0.0; }

  virtual void updateGhostPressure(Real p_ghost) { _p_ghost = p_ghost; }
  virtual void applyDirichletBC(Real & res) override { res = _v - _v_bc.value(_sim.time(), Point()); }

protected:
  const Function & _v_bc;
  const Function & _T_bc;
  Real _p_ghost;
};

/**
 * vBCEdgeInlet is for setting velocity boundary at the 'inlet' side of a pipe.
 * Here, 'inlet' does not mean flow will defintely go into the flow channel, while it is
 * more a geometric definition.
 */
class vBCEdgeInlet : public vBCEdge
{
public:
  vBCEdgeInlet(const InputParameters & parameters);
  virtual ~vBCEdgeInlet() {}

  virtual void computeFluxes() override
  {
    Real v_bc = _v_bc.value(_sim.time(), Point());
    Real T_bc = vBCEdge::T_bc();
    Real rho_bc = _eos->rho_from_p_T(_p_ghost, T_bc);
    Real h_bc = _eos->h_from_p_T(_p_ghost, T_bc);

    if (v_bc > 0.0)
    {
      _T_edge = T_bc;
      _mass_flux = v_bc * rho_bc;
      _enthalpy_flux = v_bc * rho_bc * h_bc;
    }
    else
    {
      _T_edge = _e_cell->T_w();
      _mass_flux = v_bc * _e_cell->rho_w();
      _enthalpy_flux = v_bc * _e_cell->rhoh_w();
    }
  }

  virtual Real dv_dx() override
  {
    Real v_bc = _v_bc.value(_sim.time(), Point());
    return (v_bc > 0.0) ? 0.0 : (_e_edge->v() - _v) / _e_cell->dL();
  }
  virtual Real dp_dx() override { return (_e_cell->p() - _p_ghost) / _dL_edge; }
  virtual Real dp_dx_friction() override { return 0.5 * _e_cell->dp_friction() / _dL_edge; }
  virtual Real gravity() override { return 0.5 * _e_cell->dp_gravity() / _dL_edge; }

  virtual Real rho_edge() final { return _e_cell->rho(); }
};

/**
 * vBCEdgeOutlet is for setting velocity boundary at the 'outlet' side of a pipe.
 * Similarly, 'outlet' is a geometric concept.
 */
class vBCEdgeOutlet : public vBCEdge
{
public:
  vBCEdgeOutlet(const InputParameters & parameters);
  virtual ~vBCEdgeOutlet() {}

  virtual void computeFluxes() override
  {
    Real v_bc = _v_bc.value(_sim.time(), Point());
    Real T_bc = vBCEdge::T_bc();
    Real rho_bc = _eos->rho_from_p_T(_p_ghost, T_bc);
    Real h_bc = _eos->h_from_p_T(_p_ghost, T_bc);

    if (v_bc > 0.0)
    {
      _T_edge = _w_cell->T_e();
      _mass_flux = v_bc * _w_cell->rho_e();
      _enthalpy_flux = v_bc * _w_cell->rhoh_e();
    }
    else
    {
      _T_edge = T_bc;
      _mass_flux = v_bc * rho_bc;
      _enthalpy_flux = v_bc * rho_bc * h_bc;
    }
  }

  virtual Real dv_dx() override
  {
    Real v_bc = _v_bc.value(_sim.time(), Point());
    return (v_bc > 0.0) ? (v_bc - _w_edge->v() / _w_cell->dL()) : 0.0;
  }
  virtual Real dp_dx() override { return (_p_ghost - _w_cell->p()) / _dL_edge; }
  virtual Real dp_dx_friction() override { return 0.5 * _w_cell->dp_friction() / _dL_edge; }
  virtual Real gravity() override { return 0.5 * _w_cell->dp_gravity() / _dL_edge; }

  virtual Real rho_edge() final { return _w_cell->rho(); }
};

/**
 * snjShadowEdge is a very special type of Edge.
 * When two pipes are connected with a SingleJunction (snjEdge), there are two boundary edges
 * overlapping. One of the two overlapping edges is the snjEdge, which handles the real connections,
 * computations, etc. There is therefore a redundant edge, which you cannot and should not handle
 * the connections, otherwise the connected cells got confused who they are really talking to. Thus
 * the snjShadowEdge design, which just mimics what the 'real' edge does.
 */
class snjShadowEdge : public EdgeBase
{
public:
  snjShadowEdge(const InputParameters & parameters) : EdgeBase(parameters),
  _real_edge(_pars.get<snjEdge*>("real_edge"))
  {}
  virtual ~snjShadowEdge() {}

  // a shadow edge mimics what the 'real' edge does
  virtual void computeFluxes() override { /*_real_edge->computeFluxes();*/ }
  virtual Real T_edge() override final { return _real_edge->T_edge(); }
  virtual Real mass_flux() override final { return _real_edge->mass_flux(); }
  virtual Real enthalpy_flux() override final { return _real_edge->enthalpy_flux(); }
  virtual Real dv_dt() override final { return 0.0; }
  virtual Real dv_dx() override final { return _real_edge->dv_dx(); }
  virtual Real dp_dx() override final { return _real_edge->dp_dx(); }
  virtual Real dp_dx_friction() override { return _real_edge->dp_dx_friction(); }
  virtual Real gravity() override { return _real_edge->gravity(); }
  virtual Real rho_edge() override final { return _real_edge->rho_edge(); }

  // Key implementation: shadow edge has the same velocity as the real edge
  virtual void applyDirichletBC(Real & res) override final { res = _v - _real_edge->v(); }

  virtual void setExtendedNeighborEdges() override final
  {
    // shadow edge has no connected cells but follow the real edge
    _connected_DOFs.insert(_vDOF);
    _connected_DOFs.insert(_real_edge->vDOF());
  }

protected:
  snjEdge * _real_edge;
};

/**
 * pBCEdge derives from EdgeBase, and serves as the base class for pressure boundary type of
 * implementation. For pressure boundary, p_bc and T_bc are required, while velocity is part of the
 * results.
 */
class pBCEdge : public boundaryEdge
{
public:
  pBCEdge(const InputParameters & parameters);
  virtual ~pBCEdge() {}

  virtual Real p_bc() override { return _p_bc.value(_sim.time(), Point()); }
  virtual Real T_bc() override { return _T_bc.value(_sim.time(), Point()); }

protected:
  const Function & _p_bc;
  const Function & _T_bc;
};

/**
 * pBCEdgeInlet is for setting pressure boundary at the 'inlet' side of a pipe.
 */
class pBCEdgeInlet : public pBCEdge
{
public:
  pBCEdgeInlet(const InputParameters & parameters);
  virtual ~pBCEdgeInlet() {}

  virtual void computeFluxes() override
  {
    Real p_bc = pBCEdge::p_bc();
    Real T_bc = pBCEdge::T_bc();
    Real rho_bc = _eos->rho_from_p_T(p_bc, T_bc);
    Real h_bc = _eos->h_from_p_T(p_bc, T_bc);

    if (_v > 0.0)
    {
      _T_edge = T_bc;
      _mass_flux = _v * rho_bc;
      _enthalpy_flux = _v * rho_bc * h_bc;
    }
    else
    {
      _T_edge = _e_cell->T_w();
      _mass_flux = _v * _e_cell->rho_w();
      _enthalpy_flux = _v * _e_cell->rhoh_w();
    }
  }

  virtual Real dv_dx() override { return (_v > 0.0) ? 0.0 : (_e_edge->v() - _v) / _e_cell->dL(); }
  virtual Real dp_dx() override { return (_e_cell->p() - pBCEdge::p_bc()) / _dL_edge; }
  virtual Real dp_dx_friction() override { return 0.5 * _e_cell->dp_friction() / _dL_edge; }
  virtual Real gravity() override { return 0.5 * _e_cell->dp_gravity() / _dL_edge; }

  virtual Real rho_edge() override final { return _e_cell->rho(); }
};

/**
 * pBCEdgeOutlet is for setting pressure boundary at the 'outlet' side of a pipe.
 */
class pBCEdgeOutlet : public pBCEdge
{
public:
  pBCEdgeOutlet(const InputParameters & parameters);
  virtual ~pBCEdgeOutlet() {}

  virtual void computeFluxes() override
  {
    Real p_bc = pBCEdge::p_bc();
    Real T_bc = pBCEdge::T_bc();
    Real rho_bc = _eos->rho_from_p_T(p_bc, T_bc);
    Real h_bc = _eos->h_from_p_T(p_bc, T_bc);

    if (_v > 0.0)
    {
      _T_edge = _w_cell->T_e();
      _mass_flux = _v * _w_cell->rho_e();
      _enthalpy_flux = _v * _w_cell->rhoh_e();
    }
    else
    {
      _T_edge = T_bc;
      _mass_flux = _v * rho_bc;
      _enthalpy_flux = _v * rho_bc * h_bc;
    }
  }

  virtual Real dv_dx() override { return (_v > 0.0) ? (_v - _w_edge->v()) / _w_cell->dL() : 0.0; }
  virtual Real dp_dx() override { return (pBCEdge::p_bc() - _w_cell->p()) / _dL_edge; }
  virtual Real dp_dx_friction() override { return 0.5 * _w_cell->dp_friction() / _dL_edge; }
  virtual Real gravity() override { return 0.5 * _w_cell->dp_gravity() / _dL_edge; }

  virtual Real rho_edge() override final { return _w_cell->rho(); }
};

/**
 * brvEdgeInlet is for connecting the 'inlet' side of a pipe and a VolumeBranch (BranchCell),
 * so the positive flow is pointing from BranchCell (west side) to the OneDCell (east side).
 */
class brvEdgeInlet : public EdgeBase
{
public:
  brvEdgeInlet(const InputParameters & parameters);
  virtual ~brvEdgeInlet() {}

  virtual void computeFluxes() override
  {
    if (_v > 0.0)
    {
      _T_edge = _w_cell->T_e();
      _mass_flux = _v * _w_cell->rho_e();
      _enthalpy_flux = _v * _w_cell->rhoh_e();
    }
    else
    {
      _T_edge = _e_cell->T_w();
      _mass_flux = _v * _e_cell->rho_w();
      _enthalpy_flux = _v * _e_cell->rhoh_w();
    }
  }

  virtual Real dv_dx() override { return (_v > 0.0) ? 0.0 : (_e_edge->v() - _v) / _e_cell->dL(); }
  virtual Real dp_dx() override { return (_e_cell->p() - _w_cell->p()) / _dL_edge; }
  virtual Real dp_dx_friction() override { return 0.5 * _e_cell->dp_friction() / _dL_edge; }
  virtual Real gravity() override { return 0.5 * _e_cell->dp_gravity() / _dL_edge; }

  // TODO: update to volume based average
  virtual Real rho_edge() final { return 0.5 * (_w_cell->rho() + _e_cell->rho()); }
};

/**
 * brvEdgeOutlet is for connecting the 'outlet' side of a pipe and a VolumeBranch (BranchCell),
 * so the positive flow is pointing from OneDCell (west side) to the BranchCell (east side).
 */
class brvEdgeOutlet : public EdgeBase
{
public:
  brvEdgeOutlet(const InputParameters & parameters);
  virtual ~brvEdgeOutlet() {}

  virtual void computeFluxes() override
  {
    if (_v > 0.0)
    {
      _T_edge = _w_cell->T_e();
      _mass_flux = _v * _w_cell->rho_e();
      _enthalpy_flux = _v * _w_cell->rhoh_e();
    }
    else
    {
      _T_edge = _e_cell->T_w();
      _mass_flux = _v * _e_cell->rho_w();
      _enthalpy_flux = _v * _e_cell->rhoh_w();
    }
  }

  virtual Real dv_dx() override { return (_v > 0.0) ? (_v - _w_edge->v()) / _w_cell->dL() : 0.0; }
  virtual Real dp_dx() override { return (_e_cell->p() - _w_cell->p()) / _dL_edge; }
  virtual Real dp_dx_friction() override { return 0.5 * _w_cell->dp_friction() / _dL_edge; }
  virtual Real gravity() override { return 0.5 * _w_cell->dp_gravity() / _dL_edge; }

  // TODO: update to volume based average
  virtual Real rho_edge() final { return 0.5 * (_w_cell->rho() + _e_cell->rho()); }
};
