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
  params.addRequiredParam<std::string>("var_type", "The type of variable to be returned");
  params.addClassDescription("Get velocity from LBM user object.");

  return params;
}

LBMDataAux::LBMDataAux(const InputParameters & parameters)
  : AuxKernel(parameters), 
  _lbm_uo(getUserObject<LatticeBoltzmann>("lbm_uo")),
  _var_type(getParam<std::string>("var_type"))
{
}

Real
LBMDataAux::computeValue()
{
  if (isNodal())
  {
    if (_var_type == "vel_x")
      return _lbm_uo.getUx(_current_node->id());
    else if (_var_type == "vel_y")
      return _lbm_uo.getUy(_current_node->id());
    else if (_var_type == "vel_z")
      return _lbm_uo.getUz(_current_node->id());
    else if (_var_type == "speed")
      return _lbm_uo.getSpeed(_current_node->id());
    else if (_var_type == "density")
      return _lbm_uo.getDensity(_current_node->id());
    else
      mooseError("Unknown variable type: " + _var_type);
  }
  else
  {
    mooseError("LBMDataAux kernel must be defined as nodal");
  }
}
