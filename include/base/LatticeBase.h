//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <iostream>
#include <vector>
#include <torch/torch.h>

// Forward declaration
class LatticeBoltzmann;
class LatticeBoltzmannCore;

/**
 * Base object that contains all necessary tensors including probability functions and macroscopic
 * variables Every other core operator will contain this object
 */

class LatticeBase
{
public:
  LatticeBase();
  void InitVars(const long long &,
                const long long &,
                const long long &,
                long long,
                const double &,
                const double &,
                const double &,
                const double &);

  LatticeBase & operator=(const LatticeBase & other);
  LatticeBase(const LatticeBase & other);
  void createTensor(/*const Params &*/);

protected:
  friend class LatticeBoltzmann;
  friend class LatticeBoltzmannCore;

  // Simulation dimensions
  long long _nx;
  long long _ny;
  long long _nz; // in case of 2D simulation, _nz = 1

  std::vector<long long> _f_sizes;
  std::vector<long long> _u_sizes;

  long long _q;                        // number of stream directions (including stationary)
  double _taus;                        // BGK relaxation parameter related to viscosity
  double _tauq;                        // BGK relaxation parameter related to momentum
  double _taud;                        // BGK relaxation parameter related to energy
  double _initial_density;             // initial density
  double _inlet_density;
  double _outlet_density;             // outlet density
  const double _c_s = 1.0 / sqrt(3.0); // lattice speed of sound
  const double _c_s_2 = _c_s * _c_s;
  const double _c_s_4 = _c_s_2 * _c_s_2;

  // Main data tensors
  torch::Tensor _mesh;
  torch::Tensor _f;
  torch::Tensor _f_copy;
  torch::Tensor _feq;
  torch::Tensor _ux;
  torch::Tensor _uy;
  torch::Tensor _uz;
  torch::Tensor _u;
  torch::Tensor _u_old;
  torch::Tensor _rho;
  torch::Tensor _Fx;
  torch::Tensor _Fy;
  torch::Tensor _Fz;
};
