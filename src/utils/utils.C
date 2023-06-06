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
