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
class LBMVelocityAux : public AuxKernel
{
public:
  static InputParameters validParams();

  LBMVelocityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  /// the coupled user object
  const LatticeBoltzmann & _lbm_uo;
};
