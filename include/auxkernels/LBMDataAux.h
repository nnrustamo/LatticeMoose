//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
#include "LatticeBoltzmann.h"

/**
 * AuxKernel to project LBM velocities to the moose mesh
 */
class LBMDataAux : public AuxKernel
{
public:
  static InputParameters validParams();

  LBMDataAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// the coupled user object
  const LatticeBoltzmann & _lbm_uo;

  enum class VarType
  {
    vel_x,
    vel_y,
    vel_z,
    speed,
    density
  } _var_type;
};
