//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LBMPressureAux.h"

registerMooseObject("MooseApp", LBMPressureAux);

InputParameters
LBMPressureAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<UserObjectName>("lbm_uo", "The name of the LBM object to use");
  params.addClassDescription("Get velocity from LBM user object.");
  return params;
}

LBMPressureAux::LBMPressureAux(const InputParameters & parameters)
  : AuxKernel(parameters), _lbm_uo(getUserObject<LatticeBoltzmann>("lbm_uo"))
{
}

Real
LBMPressureAux::computeValue()
{
  if (isNodal())
    return _lbm_uo.getPressure(*_current_node);
  else
    return _lbm_uo.getPressure(_current_elem->centroid());
}
