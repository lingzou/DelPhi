#include "DelPhiTypes.h"
#include "MooseError.h"

namespace DELPHI
{

template <>
EEndType
stringToEnum(const std::string & str)
{
  std::string str_upper(str);
  std::transform(str_upper.begin(), str_upper.end(), str_upper.begin(), ::toupper);

  if (str_upper == "IN")  return IN;
  else if (str_upper == "OUT")  return OUT;
  else mooseError("Unknown end type: " + str);

  // should not reach here
  return INVALID;
}

template <>
WallFrictionModel
stringToEnum(const std::string & str)
{
  std::string str_upper(str);
  std::transform(str_upper.begin(), str_upper.end(), str_upper.begin(), ::toupper);

  if (str_upper == "CONST_FRICTION")  return CONST_FRICTION;
  else if (str_upper == "WELANDER")  return WELANDER;
  else mooseError("Unknown end type: " + str);

  // should not reach here
  return INVALID_WFMODEL;
}

void
getConnection(const std::string & str, std::string & comp_name, EEndType & end_type)
{
  size_t bpos = str.find('(');
  if (bpos == std::string::npos)
    mooseError("Missing left parenthesis '('");

  size_t epos = str.find(')', bpos);
  if (epos == std::string::npos)
    mooseError("Missing right parenthesis '('");

  comp_name = str.substr(0, bpos);

  std::string end = str.substr(bpos + 1, epos - bpos - 1);

  end_type = stringToEnum<EEndType>(end);
}

}
