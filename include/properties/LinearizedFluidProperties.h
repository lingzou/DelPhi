#pragma once

#include "SinglePhaseFluidProperties.h"

class LinearizedFluidProperties : public SinglePhaseFluidProperties
{
public:
  static InputParameters validParams();

  LinearizedFluidProperties(const InputParameters & parameters);
  virtual ~LinearizedFluidProperties() {}

  virtual std::string fluidName() const override { return "linearized_fluid"; }

  virtual Real rho_from_p_T(Real p, Real T) const override;
  virtual Real h_from_p_T(Real p, Real T) const override;

protected:
  Real _rho_0;
  Real _p_0;
  Real _T_0;
  Real _drho_dp;
  Real _drho_dT;
  Real _cp;
};
