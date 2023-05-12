#pragma once

#include <string>

namespace DELPHI
{

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

void getConnection(const std::string & str, std::string & comp_name, EEndType & end_type);
}
