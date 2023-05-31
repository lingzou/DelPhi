#include "LinearizedFluidProperties.h"

registerMooseObject("delphiApp", LinearizedFluidProperties);

InputParameters
LinearizedFluidProperties::validParams()
{
  InputParameters params = SinglePhaseFluidProperties::validParams();
  params.addRequiredParam<Real>("rho_0", "Density at p_0 and T_0");
  params.addRequiredParam<Real>("p_0", "Reference pressure");
  params.addRequiredParam<Real>("T_0", "Reference temperature");
  params.addRequiredParam<Real>("drho_dp", "Constant value for density ");
  params.addRequiredParam<Real>("drho_dT", "Reference temperature");
  params.addRequiredParam<Real>("cp", "Constant value for specific heat at constant pressure");

  params.addClassDescription("Fluid properties for a simple linearized fluid with a constant drho_dp and drho_dT");
  return params;
}

LinearizedFluidProperties::LinearizedFluidProperties(const InputParameters & parameters)
  : SinglePhaseFluidProperties(parameters),
    _rho_0(getParam<Real>("rho_0")),
    _p_0(getParam<Real>("p_0")),
    _T_0(getParam<Real>("T_0")),
    _drho_dp(getParam<Real>("drho_dp")),
    _drho_dT(getParam<Real>("drho_dT")),
    _cp(getParam<Real>("cp"))
{
}

Real
LinearizedFluidProperties::rho_from_p_T(Real p, Real T) const
{
  return _rho_0 + (p - _p_0) * _drho_dp + (T - _T_0) * _drho_dT;
}

Real
LinearizedFluidProperties::h_from_p_T(Real /*p*/, Real T) const
{
  return _cp * T;
}
