#pragma once

#include <string>

namespace DELPHI
{

enum WallFrictionModel
{
  CONST_FRICTION  = 0,
  WELANDER        = 1,
  INVALID_WFMODEL = 99
};

enum EEndType
{
  IN  = 0,
  OUT = 1,
  INVALID = 99
};

template <typename T>
T stringToEnum(const std::string & str);
template <>
EEndType stringToEnum<EEndType>(const std::string & s);
template <>
WallFrictionModel stringToEnum<WallFrictionModel>(const std::string & s);

void getConnection(const std::string & str, std::string & comp_name, EEndType & end_type);
}
