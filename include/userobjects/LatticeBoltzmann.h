//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "LatticeBoltzmannCore.h"
#include "LatticeStencil.h"
#include "libmesh/utility.h"

class LatticeBoltzmann : public GeneralUserObject
{
public:
  static InputParameters validParams();

  LatticeBoltzmann(const InputParameters & parameters);

  void initialize();
  void execute();
  void finalize();
  void logStep();

  Real getSpeed(Point p) const;
  Real getPressure(Point p) const;

protected:
  LatticeBoltzmannCore _simulation_object;
  long long _n_subcycles; // LB iterations per Moose iterations
  double _tolerance;
  bool _isConverged = false;
  long long _tsteps = 0;
};
