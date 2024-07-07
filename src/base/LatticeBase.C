//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LatticeBase.h"

/**
 * Source file for base of lattice Boltzmann
 * Stores all the tensor data
 */

LatticeBase::LatticeBase()
{
  /**
   * Default constructor
   */
}

void
LatticeBase::InitVars(const long long & nx,
                      const long long & ny,
                      const long long & nz,
                      long long q,
                      const double & taus,
                      const double & initial_density,
                      const double & inlet_density,
                      const double & outlet_density)
{
  /**
   * Initialize the variables
   */
  _nx = nx;
  _ny = ny;
  _nz = nz;
  _q = q;
  _taus = taus;
  _initial_density = initial_density;
  _inlet_density = inlet_density;
  _outlet_density = outlet_density;

  if (_nz == 0)
    _nz = 1;

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
  _f_copy = torch::zeros(_f_sizes, torch::kFloat64);
  _feq = torch::zeros(_f_sizes, torch::kFloat64);
}

LatticeBase &
LatticeBase::operator=(const LatticeBase & other)
{
  /**
   * Copy sssignment operator
   */

  if (this == &other)
  {
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
  _f_copy = other._f_copy.clone();
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

LatticeBase::LatticeBase(const LatticeBase & other)
{
  /**
   * Copy Constructor
   *
   */
  // the initial state of the lattice is not set
  if (other._mesh.size(0) != 0)
  {

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
    _f_copy = other._f_copy.clone();
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
}

void LatticeBase::createTensor(/*const Params& options*/)
{
  /**
   * This function will be used for tensor settings
   */
}
