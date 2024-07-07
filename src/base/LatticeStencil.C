//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LatticeStencil.h"

StencilBase &
StencilBase::operator=(const StencilBase & other)
{
  /**
   * Copy Assigment operator
   */
  if (this == &other)
  {
    return *this;
  }
  else if (other._q != 0)
  {
    _q = other._q;

    // Copy tensors
    _ex = other._ex.clone();
    _ey = other._ey.clone();
    _ez = other._ez.clone();
    _weights = other._weights.clone();
    _op = other._op.clone();
    _M = other._M.clone();
    _M_inv = other._M_inv.clone();
    _S = other._S.clone();
    _input = other._input.clone();
    _output = other._output.clone();
    _neutral = other._neutral.clone();
    _traverseM_y = other._traverseM_y.clone();
    _traverseM_z = other._traverseM_z.clone();
  }
  return *this;
}

StencilBase::StencilBase(const StencilBase & other)
{
  /**
   * Copy constructor
   */
  if (other._q != 0)
  {
    _q = other._q;

    // Copy tensors
    _ex = other._ex.clone();
    _ey = other._ey.clone();
    _ez = other._ez.clone();
    _weights = other._weights.clone();
    _op = other._op.clone();
    _M = other._M.clone();
    _M_inv = other._M_inv.clone();
    _S = other._S.clone();
    _input = other._input.clone();
    _output = other._output.clone();
    _traverseM_y = other._traverseM_y.clone();
    _traverseM_z = other._traverseM_z.clone();
  }
}

D2Q9::D2Q9() : StencilBase()
{
  _q = 9;
  _ex = torch::tensor({0, 1, 0, -1, 0, 1, -1, -1, 1}, torch::kInt64);
  _ey = torch::tensor({0, 0, 1, 0, -1, 1, 1, -1, -1}, torch::kInt64);
  _ez = torch::tensor({0, 0, 0, 0, 0, 0, 0, 0, 0}, torch::kInt64);
  _weights = torch::tensor({4.0 / 9.0,
                            1.0 / 9.0,
                            1.0 / 9.0,
                            1.0 / 9.0,
                            1.0 / 9.0,
                            1.0 / 36.0,
                            1.0 / 36.0,
                            1.0 / 36.0,
                            1.0 / 36.0},
                           torch::kFloat64);

  _op = torch::tensor({0, 3, 4, 1, 2, 7, 8, 5, 6}, torch::kInt64);
  _M = torch::tensor({{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                      {-4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0},
                      {4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0},
                      {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0},
                      {0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0},
                      {0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0},
                      {0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0},
                      {0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0},
                      {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0}},
                     torch::kFloat64);
  _M_inv = torch::linalg::inv(_M);

  // placeholders for                                            taud        tauq,       taud, tauq
  // taus         taus
  _S = torch::diag(torch::tensor({1.0 / 1.0,
                                  1.0 / 1.1,
                                  1.0 / 1.2,
                                  1.0 / 1.0000,
                                  1.0 / 1.0000,
                                  1.0 / 1.0000,
                                  1.0 / 1.0000,
                                  1.0 / 1.0000,
                                  1.0 / 1.0000},
                                 torch::kFloat64));
    /**
   * unknown distribution functions at the input and output in flow direction (x)
   * the order of the input and output is important
   * E.g. _input[0] -> _output[0]
   */   
  _input = torch::tensor({1, 5, 8}, torch::kInt64);
  _output = torch::tensor({3, 7, 6}, torch::kInt64);
  // distributions that are neither input nor output
  _neutral = torch::tensor({0, 2, 4}, torch::kInt64);
  // traverse momentum correction tensor
  _traverseM_y = torch::tensor({2, 4}, torch::kInt64);
  _traverseM_z = torch::tensor({0, 0}, torch::kInt64); 
}

D3Q19::D3Q19() : StencilBase()
{
  
  _q = 19;
  _ex = torch::tensor({0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1,  -1}, torch::kInt64);
  _ey = torch::tensor({0,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1,  -1, 1,  -1}, torch::kInt64);
  _ez = torch::tensor({0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,   0,  0,  0}, torch::kInt64);
  _weights = torch::tensor(
      {
          1.0 / 3.0,  1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0, 1.0 / 18.0,
          1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
          1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
      },
      torch::kFloat64);

  _op = torch::tensor({0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15},
                      torch::kInt64);

  _M = torch::tensor(
      {{1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
       {-30., -11., -11., -11., -11., -11., -11., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8.},
       {12., -4., -4., -4., -4., -4., -4., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
       {0., 1., -1., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0.},
       {0., -4., 4., 0., 0., 0., 0., 1., -1., 1., -1., 1., -1., 1., -1., 0., 0., 0., 0.},
       {0., 0., 0., 1., -1., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0., 1., -1., 1., -1.},
       {0., 0., 0., -4., 4., 0., 0., 1., 1., -1., -1., 0., 0., 0., 0., 1., -1., 1., -1.},
       {0., 0., 0., 0., 0., 1., -1., 0., 0., 0., 0., 1., 1., -1., -1., 1., 1., -1., -1.},
       {0., 0., 0., 0., 0., -4., 4., 0., 0., 0., 0., 1., 1., -1., -1., 1., 1., -1., -1.},
       {0., 2., 2., -1., -1., -1., -1., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2.},
       {0., -4., -4., 2., 2., 2., 2., 1., 1., 1., 1., 1., 1., 1., 1., -2., -2., -2., -2.},
       {0., 0., 0., 1., 1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0.},
       {0., 0., 0., -2., -2., 2., 2., 1., 1., 1., 1., -1., -1., -1., -1., 0., 0., 0., 0.},
       {0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0., 0., 0., 0., 0.},
       {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1.},
       {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -1., -1., 1., 0., 0., 0., 0.},
       {0., 0., 0., 0., 0., 0., 0., 1., -1., 1., -1., -1., 1., -1., 1., 0., 0., 0., 0.},
       {0., 0., 0., 0., 0., 0., 0., -1., -1., 1., 1., 0., 0., 0., 0., 1., -1., 1., -1.},
       {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1., -1., -1., 1., 1.}},
      torch::kFloat64);
  _M_inv = torch::linalg::inv(_M);

  _S = torch::diag(torch::tensor({1. / 1.,
                                  1. / 1.19,
                                  1. / 1.4,
                                  1. / 1.4,
                                  1. / 1.0000,
                                  1. / 1.,
                                  1. / 1.0000,
                                  1. / 1.,
                                  1. / 1.0000,
                                  1. / 1.0000,
                                  1. / 1.4,
                                  1. / 1.000,
                                  1. / 1.4,
                                  1. / 1.000,
                                  1. / 1.000,
                                  1. / 1.000,
                                  1. / 1.98,
                                  1. / 1.98,
                                  1. / 1.98},
                                 torch::kFloat64));

  /**
   * unknown distribution functions at the input and output in flow direction (x)
   * the order of the input and output is important
   * E.g. _input[0] -> _output[0]
   */                                 
  _input = torch::tensor({5, 11, 12, 15, 16}, torch::kInt64);
  _output = torch::tensor({6, 14, 13, 18, 17}, torch::kInt64);
  // distributions that are neither input nor output
  _neutral = torch::tensor({0, 1, 2, 3, 4, 7, 8, 9, 10}, torch::kInt64);
  // traverse momentum correction tensor
  _traverseM_y = torch::tensor({3, 4, 7, 8, 9, 10}, torch::kInt64); 
  _traverseM_z = torch::tensor({1, 2, 7, 8, 9, 10}, torch::kInt64);
};
