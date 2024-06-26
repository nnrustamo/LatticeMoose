//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#pragma once

#include "LatticeBase.hh"

/**
 * Lattice configurations heaeder file
 */
struct StencilBase
{
    /**
     * Base to create lattice configureation objects
     */
public:
    int64_t _q = 0;
    torch::Tensor _ex;
    torch::Tensor _ey;
    torch::Tensor _ez;
    torch::Tensor _weights;
    torch::Tensor _op;
    torch::Tensor _M;
    torch::Tensor _M_inv;
    torch::Tensor _S;

public:
    StencilBase() = default;   
    StencilBase& operator=(const StencilBase& other);
    StencilBase(const StencilBase& other);
    ~StencilBase() = default;
};

struct D2Q9 : public StencilBase
{
    /**
     * 2-dimensional 9 velocity lattice configuration
     */
public:
    D2Q9();
};

struct D3Q19 : public StencilBase
{
    /**
     * 3-dimensional 19 velocity lattice configuration
     */
public:
    D3Q19();
};
