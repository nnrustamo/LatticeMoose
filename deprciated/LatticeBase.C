//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#include "LatticeBase.h"
#include "UserObject.h"

/**
 * Source file for base of lattice Boltzmann 
 * Stores all the tensor data
 */

InputParameters
LatticeBase::validParams()
{
  InputParameters params = UserObject::validParams();
  params.addRequiredParam<long long>("q", "The number of streaming directions");
  params.addRequiredParam<int64_t>("nx", "Mesh size in x directions");
  params.addRequiredParam<int64_t>("ny", "Mesh size in y directions");
  params.addParam<int64_t>("nz", 1, "Mesh size in z directions");
  params.addParam<double>("taus", 1, "Relaxation parameter");
  params.addParam<double>("initial_density", 1, "Initial lattice density");
  params.addClassDescription("Base lattice Boltzmann object that stores large simulation tensors");
  return params;
}

LatticeBase::LatticeBase(const InputParameters & parameters)
            : UserObject(parameters), 
            _q(getParam<std::size_t>("q")),
            _nx(getParam<std::size_t>("nx")),
            _ny(getParam<std::size_t>("ny")),
            _nz(getParam<std::size_t>("nz")),
            _taus(getParam<std::size_t>("taus")),
            _initial_density(getParam<std::size_t>("initial_density"))
{   
    if (_nz == 0)
    {
        _nz = 1;
    }

    _f_sizes = {_nz, _ny, _nx, _q};
    _u_sizes = {_nz, _ny, _nx};

    //
    _mesh = torch::ones(_u_sizes, torch::kInt16);
    _ux = torch::zeros(_u_sizes, torch::kFloat64);
    _uy = torch::zeros(_u_sizes, torch::kFloat64);
    _uz = torch::zeros(_u_sizes, torch::kFloat64);
    _u = torch::zeros(_u_sizes, torch::kFloat64);
    _u_old = torch::zeros(_u_sizes, torch::kFloat64);
    _rho = torch::zeros(_u_sizes, torch::kFloat64);
    _rho.fill_(_initial_density);

    _f = torch::zeros(_f_sizes, torch::kFloat64);
    _feq = torch::zeros(_f_sizes, torch::kFloat64);
}

LatticeBase& LatticeBase::operator=(const LatticeBase& other) 
{
    /**
     * Copy sssignment operator 
     */

    if (this == &other) {
        return *this;
    }

    // Copy data
    _nx = other._nx;
    _ny = other._ny;
    _nz = other._nz;
    _q = other._q;
    _taus = other._taus;
    _initial_density = other._initial_density;

    _f_sizes = other._f_sizes;
    _u_sizes = other._u_sizes;

    _mesh = other._mesh.clone();
    _f = other._f.clone();
    _feq = other._feq.clone();
    _ux = other._ux.clone();
    _uy = other._uy.clone();
    _uz = other._uz.clone();
    _u = other._u.clone();
    _u_old = other._u_old.clone();
    _rho = other._rho.clone();
    if (other._Fx.sizes()[0] != 0)
    {
        _Fx = other._Fx.clone();
        _Fy = other._Fy.clone();
        _Fz = other._Fz.clone();
    }
    return *this;
}

LatticeBase::LatticeBase(const LatticeBase& other) 
{
    /**
     * Copy Constructor
     */

    _nx = other._nx;
    _ny = other._ny;
    _nz = other._nz;
    _q = other._q;
    _taus = other._taus;
    _initial_density = other._initial_density;

    _f_sizes = other._f_sizes;
    _u_sizes = other._u_sizes;

    _mesh = other._mesh.clone();
    _f = other._f.clone();
    _feq = other._feq.clone();
    _ux = other._ux.clone();
    _uy = other._uy.clone();
    _uz = other._uz.clone();
    _u = other._u.clone();
    _u_old = other._u_old.clone();
    _rho = other._rho.clone();
    if (other._Fx.sizes()[0] != 0)
    {
        _Fx = other._Fx.clone();
        _Fy = other._Fy.clone();
        _Fz = other._Fz.clone();
    }
}


void 
LatticeBase::createTensor(/*const Params& options*/)
{
    /**
     * This function will be used for tensor settings
     */
}
