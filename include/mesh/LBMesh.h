//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneratedMesh.h"

/**
 * Generate structured mesh for Lattice Boltzmann simulations
 */

class LBMesh : public GeneratedMesh
{
public:
  static InputParameters validParams();

  LBMesh(const InputParameters & parameters);
  LBMesh(const LBMesh & /* other_mesh */) = default;

  // No copy
  LBMesh & operator=(const LBMesh & other_mesh) = delete;
  virtual void buildMesh() override;
  unsigned int getNx() const;
  unsigned int getNy() const;
  unsigned int getNz() const;
  MooseEnum getDim() const;

protected:   
    // Number of nodes in x, y, z direction
    unsigned int _lbnx, _lbny, _lbnz;
};

