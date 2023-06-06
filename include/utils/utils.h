#pragma once

#include "Moose.h"

namespace DELPHI
{

Real BDF2_ddt(Real u, Real u_o, Real u_oo, Real dt, Real dt_o);

}
