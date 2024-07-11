//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LBMDataAux.h"

registerMooseObject("MooseApp", LBMDataAux);

InputParameters
LBMDataAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>("lbm_uo", "The name of the LBM object to use");
  MooseEnum varType("vel_x vel_y vel_z speed density");
  params.addRequiredParam<MooseEnum>("var_type", varType, "The type of variable to be returned");
  params.addClassDescription("Get velocity from LBM user object.");

  return params;
}

LBMDataAux::LBMDataAux(const InputParameters & parameters)
  : AuxKernel(parameters), 
  _lbm_uo(getUserObject<LatticeBoltzmann>("lbm_uo")),
  _var_type(getParam<MooseEnum>("var_type").getEnum<VarType>())
{
  if (!isNodal())
    mooseError("LBMDataAux kernel must be defined as nodal");
}

Real
LBMDataAux::computeValue()
{
  switch (_var_type)
  {
    case VarType::vel_x:
      return _lbm_uo.getUx(_current_node->id());
    case VarType::vel_y:
      return _lbm_uo.getUy(_current_node->id());
    case VarType::vel_z:
      return _lbm_uo.getUz(_current_node->id());
    case VarType::speed:
      return _lbm_uo.getSpeed(_current_node->id());
    case VarType::density:
      return _lbm_uo.getDensity(_current_node->id());
    default:
      mooseError("Unknown variable type.");
  }
}
