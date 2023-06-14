#include "utils.h"

Real
DELPHI::BDF2_ddt(Real u, Real u_o, Real u_oo, Real dt, Real dt_o)
{
  Real rr = dt / dt_o;
  Real w0 = (2.0 * rr + 1.0) / ((rr + 1.0) * dt);
  Real w1 = -(rr + 1.0) / dt;
  Real w2 = rr * rr / ((rr + 1.0) * dt);

  return w0 * u + w1 * u_o + w2 * u_oo;
}

Real
DELPHI::SlopeLimiterIrregularMesh(Real u_W, Real u, Real u_E, Real dL_W, Real dL, Real dL_E)
{
  if (((u_E - u) * (u - u_W) <= 0.0) ||   // local extremum
      (std::fabs(u_E - u_W) < 1.e-10) )   // difference too small to have meaningful reconstruction
    return 0.0;
  else
  {
    Real a = dL_W / dL;
    Real b = dL_E / dL;
    Real fp = (1.0 + a) / (2.0 + a + b);
    Real J = u_E - u_W;
    Real f = (u - u_W) / J; // we know f must be between 0 and 1 now (the if statement)

    if (f < fp)
      return DELPHI::Berger_func1(a, f, fp, J, dL);
    else
      return DELPHI::Berger_func2(b, f, fp, J, dL);
  }
}

Real
DELPHI::Berger_func1(Real a, Real f, Real fp, Real J, Real h)
{
  if (a > 1.0e-6)
    return f * (1.0 - a / (1.0 + a) * std::pow(f / fp, 1.0 / a)) * 2.0 * J / h;
  else
    return 2.0 * f * J / h;
}

Real
DELPHI::Berger_func2(Real b, Real f, Real fp, Real J, Real h)
{
  if (b > 1.0e-6)
    return (1.0 - f) * (1.0 - b / (1.0 + b) * std::pow((1.0 - f) / (1.0 - fp), 1.0 / b)) * 2.0 * J / h;
  else
    return 2.0 * (1.0 - f) * J / h;
}
