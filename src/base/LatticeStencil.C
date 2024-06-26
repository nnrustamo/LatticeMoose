//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LatticeStencil.h"
#include "Registry.h"

StencilBase& StencilBase::operator=(const StencilBase& other) 
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
    }
    return *this;
}

StencilBase::StencilBase(const StencilBase& other) 
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
    }
}

D2Q9::D2Q9() : StencilBase(){
    _q = 9;
    _ex = torch::tensor({0, 1, 0, -1, 0, 1, -1, -1, 1}, torch::kInt64);
    _ey = torch::tensor({0, 0, 1, 0, -1, 1, 1, -1, -1}, torch::kInt64);
    _ez = torch::tensor({0, 0, 0, 0, 0, 0, 0, 0, 0}, torch::kInt64);
    _weights = torch::tensor({4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
        1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0}, torch::kFloat64);

    _op = torch::tensor({0, 3, 4, 1, 2, 7, 8, 5, 6}, torch::kInt64);
    _M = torch::tensor({
                        {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
                        {-4.0, -1.0, -1.0, -1.0, -1.0, 2.0, 2.0, 2.0, 2.0},
                        {4.0, -2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0},
                        {0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0},
                        {0.0, -2.0, 0.0, 2.0, 0.0, 1.0, -1.0, -1.0, 1.0},
                        {0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0},
                        {0.0, 0.0, -2.0, 0.0, 2.0, 1.0, 1.0, -1.0, -1.0},
                        {0.0, 1.0, -1.0, 1.0, -1.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0}}, torch::kFloat64);
    _M_inv = torch::linalg::inv(_M);
    
    // placeholders for                                            taud        tauq,       taud,       tauq       taus         taus
    _S = torch::diag(torch::tensor({1.0/1.0, 1.0/1.1, 1.0/1.2, 1.0/1.0000, 1.0/1.0000, 1.0/1.0000, 1.0/1.0000, 1.0/1.0000, 1.0/1.0000}, torch::kFloat64));
}

D3Q19::D3Q19() : StencilBase() {
    _q = 19;
    _ex = torch::tensor({0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1,}, torch::kInt64);   
    _ey = torch::tensor({0, 0,  0, 1, -1, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1,  -1, -1, 1,}, torch::kInt64);
    _ez = torch::tensor({0, 0, 0, 0, 0, 1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 0, 0, 0, 0}, torch::kInt64);
    _weights = torch::tensor({1.0 / 3.0, 1.0/ 18.0, 1.0/ 18.0, 1.0/ 18.0, 1.0/ 18.0, 1.0/ 18.0, 1.0/ 18.0,
                        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,
                        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0,}, torch::kFloat64);
    

    _op = torch::tensor({0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13,
                            16, 15, 18, 17}, torch::kInt64);

    _M = torch::tensor({
                        {1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.},
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
                        {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., -1., -1., -1., -1., 1., 1.}}, torch::kFloat64);
    _M_inv = torch::linalg::inv(_M);

    // placeholders for                                                 tauq               tauq              tauq       taus               taus               taus      taus      taus
    _S = torch::diag(torch::tensor({1./1., 1./1.19, 1./ 1.4, 1./1.4, 1./1.0000, 1./1.,  1./1.0000, 1./1., 1./1.0000, 1./1.0000, 1./1.4, 1./1.000, 1./1.4,  1./1.000, 1./1.000, 1./1.000, 1./1.98, 1./1.98, 1./1.98}, torch::kFloat64));
};
