//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LatticeBase.h"

// Forward declaration
class LatticeBoltzmann;
class LatticeBoltzmannCore;

/**
 * Lattice configurations heaeder file
 */
class StencilBase
{
  /**
   * Base to create lattice configureation objects
   */
public:
  StencilBase() = default;
  StencilBase & operator=(const StencilBase & other);
  StencilBase(const StencilBase & other);

protected:
  friend class LatticeBoltzmann;
  friend class LatticeBoltzmannCore;
  int64_t _q = 0;
  torch::Tensor _ex;
  torch::Tensor _ey;
  torch::Tensor _ez;
  torch::Tensor _weights;
  torch::Tensor _op;
  torch::Tensor _M;
  torch::Tensor _M_inv;
  torch::Tensor _S;
  // unknown distribution functions at the input and output in flow direction (x)
  torch::Tensor _input;  
  torch::Tensor _output;
  // distributions that are neither input nor output
  torch::Tensor _neutral;
  // traverse momentum correction tensor
  torch::Tensor _traverseM_y;
  torch::Tensor _traverseM_z;
};

class D2Q9 : public StencilBase
{
  /**
   * 2-dimensional 9 velocity lattice configuration
   */
public:
  D2Q9();
};

class D3Q19 : public StencilBase
{
  /**
   * 3-dimensional 19 velocity lattice configuration
   */
public:
  D3Q19();
};
