#pragma once

#include "Moose.h"

namespace DELPHI
{

Real BDF2_ddt(Real u, Real u_o, Real u_oo, Real dt, Real dt_o);

Real SlopeLimiterIrregularMesh(Real u_W, Real u, Real u_E, Real dL_W, Real dL, Real dL_E);
Real Berger_func1(Real a, Real f, Real fp, Real J, Real dx);
Real Berger_func2(Real b, Real f, Real fp, Real J, Real dx);
}
