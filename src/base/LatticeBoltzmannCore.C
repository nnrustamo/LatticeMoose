//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LatticeBoltzmannCore.h"

/**
 * Source file for core Lattice Boltzmann Method
 */

using namespace torch::indexing;

LatticeBoltzmannCore::LatticeBoltzmannCore(const LatticeBase & lattice, const StencilBase & stencil)
  : _lattice(lattice), _stencil(stencil){};

void
LatticeBoltzmannCore::generateMesh()
{
  /**
   * Mesh tensor was already created in _lattice
   * This function modifies it as necessary
   */

  // Adding simple cricrle to the mesh
  /*
  int64_t Center_x = 19;
  int64_t Center_y = 19;
  int64_t Center_z = 19;
  int64_t Radius = 5;
  */

  // Create Meshgrid
  auto x = torch::arange(0, _lattice._nx, torch::kInt32);
  auto y = torch::arange(0, _lattice._ny, torch::kInt32);
  auto z = torch::arange(0, _lattice._nz, torch::kInt32);
  std::vector<torch::Tensor> grids = torch::meshgrid({z, y, x}, "ij");

  // Extract individual grids
  torch::Tensor z_indices = grids[0];
  torch::Tensor y_indices = grids[1];
  torch::Tensor x_indices = grids[2];

  // Create boundary mask
  torch::Tensor boundary_mask = torch::zeros_like(_lattice._mesh, torch::kBool);
  boundary_mask.fill_(false);
  if (_lattice._nz == 1)
  {
    boundary_mask = (y_indices == 0) | (y_indices == _lattice._ny - 1);// | 
                  /* adds circle 
                  (torch::sqrt((y_indices - Center_y) * (y_indices - Center_y) + \
                              (x_indices - Center_x) * (x_indices - Center_x)) <= Radius);*/
  }
  else
  {
    boundary_mask = (y_indices == 0) | (y_indices == _lattice._ny - 1); // | (z_indices == 0) | (z_indices == _lattice._nz - 1);// | 
                    /* adds circle 
                    (torch::sqrt((z_indices - Center_z) * (z_indices - Center_z) + \
                                (y_indices - Center_y) * (y_indices - Center_y) + \
                                (x_indices - Center_x) * (x_indices - Center_x)) <= Radius);*/
  }

  _lattice._mesh.masked_fill_(boundary_mask, 0);

  /**
   * Nodes that are adjacent to boundary will be set to 2, this will later be used in determining
   * the nodes for bounce-back
   * This will be achieved by rolling the mesh around in streaming directions and finding where
   * boundary hit happens
   */

  torch::Tensor new_mesh = _lattice._mesh.clone();
  for (int64_t ic = 1; ic < _stencil._q; ic++)
  {
    int64_t ex = _stencil._ex[ic].item<int64_t>();
    int64_t ey = _stencil._ey[ic].item<int64_t>();
    int64_t ez = _stencil._ez[ic].item<int64_t>();
    torch::Tensor shifted_mesh = torch::roll(_lattice._mesh, {ez, ey, ex}, {0, 1, 2});
    torch::Tensor adjacent_to_boundary = (shifted_mesh == 0) & (_lattice._mesh == 1);
    new_mesh.masked_fill_(adjacent_to_boundary, 2);
  }

  // Deep copy new mesh
  _lattice._mesh = new_mesh.clone();
}

void LatticeBoltzmannCore::loadMeshFromFile(/*std::string &filename*/)
{
  /**
   * This function loads mesh from file
   */
}

void
LatticeBoltzmannCore::stream()
{
  /**
   * Generalized streaming function for 2D and 3D
   * Torch implementation automatically enables periodic BC in al directions
   */

  _lattice._f_copy = _lattice._f.clone(); // saving pre-streming distribution
  for (/* do not use unsigned int */ int i = 1; i < _stencil._q; i++)
  {
    torch::Tensor rolled_tensor = torch::roll(
        _lattice._f_copy.index(
            {Slice(), Slice(), Slice(), i}),
        /* shifts = */
        {_stencil._ez[i].item<int64_t>(),
         _stencil._ey[i].item<int64_t>(),
         _stencil._ex[i].item<int64_t>()},
        /* dims = */
        {0, 1, 2});

    _lattice._f.index(
        {Slice(), Slice(), Slice(), i}) =
        rolled_tensor;
  }

  // torch.roll will stream into the solid nodes, set them to zero
  setfSolidtoZero();
}

void
LatticeBoltzmannCore::BGKCollision()
{
  /**
   * Single relaxation time collision operator based on BGK scheme
   * https://journals.aps.org/pr/abstract/10.1103/PhysRev.94.511
   */
  _lattice._f = _lattice._f - 1 / _lattice._taus * (_lattice._f - _lattice._feq);
}

void
LatticeBoltzmannCore::MRTCollision()
{
  /**
   * Multi-relaxation time (with and wihtout Hermite regularization) collision operator
   * https://doi.org/10.1016/j.jgsce.2023.205131
   */

  if (_enableRegularization)
  {
    regularize();
  }
  torch::Tensor f_noneq = (_lattice._f - _lattice._feq).view({-1, _stencil._q}).t();
  torch::Tensor pre_collision_moments = torch::matmul(_stencil._M, f_noneq);
  torch::Tensor post_collision_moments = torch::matmul(_stencil._S, pre_collision_moments);
  torch::Tensor post_collision_f = torch::matmul(_stencil._M_inv, post_collision_moments);
  _lattice._f = _lattice._f - post_collision_f.t().view({_lattice._f_sizes});

  // MRT will result in NAN at the solid nodes, set them to zero
  setfSolidtoZero();
}

void
LatticeBoltzmannCore::computeEquilibrium()
{
  /**
   * Calculates equilibrium distribution
   */
  const int64_t q = _stencil._q;
  const torch::Tensor ex = _stencil._ex.view({1, 1, 1, q});
  const torch::Tensor ey = _stencil._ey.view({1, 1, 1, q});
  const torch::Tensor ez = _stencil._ez.view({1, 1, 1, q});
  const torch::Tensor w = _stencil._weights.view({1, 1, 1, q});
  const torch::Tensor ux = _lattice._ux.unsqueeze(3);
  const torch::Tensor uy = _lattice._uy.unsqueeze(3);
  const torch::Tensor uz = _lattice._uz.unsqueeze(3);
  const torch::Tensor rho = _lattice._rho.unsqueeze(3);

  torch::Tensor comp1 = (ex * ux + ey * uy + ez * uz) / _lattice._c_s_2;
  torch::Tensor comp2 =
      0.5 * ((ex * ux + ey * uy + ez * uz) * (ex * ux + ey * uy + ez * uz)) / _lattice._c_s_4;
  torch::Tensor comp3 = -0.5 * (ux * ux + uy * uy + uz * uz) / _lattice._c_s_2;
  _lattice._feq = w * rho * (1.0 + comp1 + comp2 + comp3);

  // Unkown behavior may occur at the solid nodes, set them to zero
  setfeqSolidtoZero();
}

void
LatticeBoltzmannCore::computeObservables()
{
  /**
   * Calculate macroscopic properties and residual
   */

  // lattice density
  _lattice._rho = torch::sum(_lattice._f, /*dim = */ 3);

  // fluxes
  torch::Tensor flux_x = torch::sum(_lattice._f * _stencil._ex, 3);
  torch::Tensor flux_y = torch::sum(_lattice._f * _stencil._ey, 3);
  torch::Tensor flux_z = torch::sum(_lattice._f * _stencil._ez, 3);

  // set pressure boundary
  _lattice._rho.index_put_({Slice(), Slice(), 0}, _lattice._inlet_density);
  _lattice._rho.index_put_({Slice(), Slice(), _lattice._nx - 1}, _lattice._outlet_density);

  // lattice velocity
  _lattice._ux = flux_x / _lattice._rho + _fBody;
  _lattice._uy = flux_y / _lattice._rho;
  _lattice._uz = flux_z / _lattice._rho;

  setObservablestoZero();
}

void
LatticeBoltzmannCore::residual()
{
  /**
   * Compute residual
   */

  _lattice._u = torch::sqrt(_lattice._ux * _lattice._ux + _lattice._uy * _lattice._uy +
                            _lattice._uz * _lattice._uz);

  double sumUsqareMinusUsqareOld =
      torch::sum(torch::abs(_lattice._u - _lattice._u_old)).item<double>();
  double sumUsquare = torch::sum(_lattice._u).item<double>();

  _residual = (sumUsquare == 0) ? 1.0 : sumUsqareMinusUsqareOld / sumUsquare;
}

void
LatticeBoltzmannCore::wallBoundary()
{
  /**
   * Halfway bounce-back
   */

  torch::Tensor mesh_expanded = _lattice._mesh.unsqueeze(-1).expand_as(_lattice._f);
  torch::Tensor boundary_mask = (mesh_expanded == 2) & (_lattice._f == 0);
  boundary_mask = boundary_mask.to(torch::kBool);
  torch::Tensor f_bounce_back = torch::zeros_like(_lattice._f);

  for (/* do not use unsigned int */ int ic = 1; ic < _stencil._q; ic++)
  {
    int64_t index = _stencil._op[ic].item<int64_t>();
    auto lattice_slice = _lattice._f_copy.index(
        {Slice(), Slice(), Slice(), index});
    auto bounce_back_slice = f_bounce_back.index(
        {Slice(), Slice(), Slice(), ic});

    // f_bounce_back.index({Slice(), Slice(),
    // Slice(), ic}) =
    //     _lattice._f.index({Slice(), Slice(),
    //     Slice(), _stencil._op[ic].item<int64_t>()});
    //  f_bounce_back.select(2, ic).copy_(_lattice._f.select(2, _stencil._op[ic].item<int64_t>()));

    bounce_back_slice.copy_(lattice_slice);
  }
  _lattice._f.index_put_({boundary_mask}, f_bounce_back.index({boundary_mask}));
}

void
LatticeBoltzmannCore::openBoundary()
{

  /**
   * Pressure boundary condition based on
   * https://doi.org/10.1063/1.869307  (2D)
   * https://doi.org/10.1088/1742-5468/2010/01/P01018 (3D)
   * There may be a way to combine the 2D and 3D into one set of equations
   */

  int64_t inlet_layer = 0;
  int64_t outlet_layer = _lattice._nx - 1;
  torch::Tensor u_inlet, u_outlet, sum_neutral, sum_output, sum_input, sum_z, sum_e_z, Nx_z, sum_y, sum_e_y, Nx_y;
  switch (_lattice._q)
  {
    case 9: // 2D
      /**
       * inlet 
       */
      u_inlet = 1.0 - (1.0 / _lattice._inlet_density) * (_lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._neutral[0]}) +
                                                          _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._neutral[1]}) +
                                                          _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._neutral[2]}) +
                                                    2.0 * (_lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[0]}) +
                                                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[1]}) +
                                                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[2]})));

      _lattice._f.index_put_({Slice(), Slice(), inlet_layer, _stencil._input[0]}, 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[0]}) + 2.0 / 3.0 * _lattice._inlet_density * u_inlet);
      _lattice._f.index_put_({Slice(), Slice(), inlet_layer, _stencil._input[1]}, 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[1]}) + 1.0 / 6.0 * _lattice._inlet_density * u_inlet -
                            0.5 * (_lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._neutral[1]}) - 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._neutral[2]})));
      _lattice._f.index_put_({Slice(), Slice(), inlet_layer, _stencil._input[2]}, 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[2]}) + 1.0 / 6.0 * _lattice._inlet_density * u_inlet +
                            0.5 * (_lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._neutral[1]}) - 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._neutral[2]})));

      /**
       * outlet
       */
      u_outlet = (1.0 / _lattice._outlet_density) * (_lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._neutral[0]}) +
                                                          _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._neutral[1]}) +
                                                          _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._neutral[2]}) +
                                                    2.0 * (_lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[0]}) +
                                                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[1]}) +
                                                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[2]}))) - 1.0;
      _lattice._f.index_put_({Slice(), Slice(), outlet_layer, _stencil._output[0]}, 
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[0]}) - 2.0 / 3.0 * _lattice._outlet_density * u_outlet);
      _lattice._f.index_put_({Slice(), Slice(), outlet_layer, _stencil._output[1]}, 
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[1]}) - 1.0 / 6.0 * _lattice._outlet_density * u_outlet + 
                            0.5 * (_lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._neutral[1]}) - 
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._neutral[2]})));
      _lattice._f.index_put_({Slice(), Slice(), outlet_layer, _stencil._output[2]},
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[2]}) - 1.0 / 6.0 * _lattice._outlet_density * u_outlet - 
                            0.5 * (_lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._neutral[1]}) - 
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._neutral[2]})));
      break;

    case 19:  // 3D
      /**
       * inlet
       */
      // velocity
      sum_neutral = _lattice._f.index({Slice(),  Slice(), inlet_layer, _stencil._neutral});
      sum_output = _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output});
      u_inlet = 1.0 - (1.0 / _lattice._inlet_density) * (torch::sum(sum_neutral, 2) + 2.0 * (torch::sum(sum_output, 2)));
      // tangential momentum
      sum_z =  _lattice._f.index({Slice(),  Slice(), inlet_layer, _stencil._traverseM_z});                                                
      sum_e_z = _stencil._ez.index({_stencil._traverseM_z});
      Nx_z = 0.5 * torch::sum(sum_z * sum_e_z, 2) - _lattice._inlet_density * _lattice._uz.index({Slice(), Slice(), inlet_layer}) * (1.0 / 3.0);

      sum_y =  _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._traverseM_y});                                                
      sum_e_y = _stencil._ey.index({_stencil._traverseM_y});
      Nx_y = 0.5 * torch::sum(sum_y * sum_e_y, 2) - _lattice._inlet_density * _lattice._uy.index({Slice(), Slice(), inlet_layer}) * (1.0 / 3.0);

      // distribution functions
      _lattice._f.index_put_({Slice(), Slice(), inlet_layer, _stencil._input[0]}, 
                          _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[0]}) + 1.0 / 3.0 * _lattice._inlet_density * u_inlet); // f[5]
      _lattice._f.index_put_({Slice(), Slice(), inlet_layer, _stencil._input[1]}, 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[1]}) + 1.0 / 6.0 * _lattice._inlet_density * 
                            (u_inlet - _lattice._uz.index({Slice(), Slice(), inlet_layer})) + Nx_z); // f[12]
      _lattice._f.index_put_({Slice(), Slice(), inlet_layer, _stencil._input[2]}, 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[2]}) + 1.0 / 6.0 * _lattice._inlet_density * 
                            (u_inlet + _lattice._uz.index({Slice(), Slice(), inlet_layer})) - Nx_z); // f[11]
      _lattice._f.index_put_({Slice(), Slice(), inlet_layer, _stencil._input[3]}, 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[3]}) + 1.0 / 6.0 * _lattice._inlet_density * 
                            (u_inlet - _lattice._uy.index({Slice(), Slice(), inlet_layer})) + Nx_y);  // f[16]
      _lattice._f.index_put_({Slice(), Slice(), inlet_layer, _stencil._input[4]}, 
                            _lattice._f.index({Slice(), Slice(), inlet_layer, _stencil._output[4]}) + 1.0 / 6.0 * _lattice._inlet_density * 
                            (u_inlet + _lattice._uy.index({Slice(), Slice(), inlet_layer})) - Nx_y); // f[15]
      /**
       * outlet
       */
      // velocity
      sum_neutral = _lattice._f.index({Slice(),  Slice(), outlet_layer, _stencil._neutral});
      sum_input = _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input});
      u_outlet = (1.0 / _lattice._outlet_density) * (torch::sum(sum_neutral, 2) + 2.0 * torch::sum(sum_input, 2)) - 1.0;
      // tangential momentum
      sum_z =  _lattice._f.index({Slice(),  Slice(), outlet_layer, _stencil._traverseM_z});
      Nx_z = 0.5 * torch::sum(sum_z * sum_e_z, 2) - _lattice._outlet_density * _lattice._uz.index({Slice(), Slice(), outlet_layer}) * (1.0 / 3.0);

      sum_y =  _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._traverseM_y});
      Nx_y = 0.5 * torch::sum(sum_y * sum_e_y, 2) - _lattice._outlet_density * _lattice._uy.index({Slice(), Slice(), outlet_layer}) * (1.0 / 3.0);

      // distribution functions
      _lattice._f.index_put_({Slice(), Slice(), outlet_layer, _stencil._output[0]}, 
                          _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[0]}) - 1.0 / 3.0 * _lattice._outlet_density * u_outlet); // f[6]
      _lattice._f.index_put_({Slice(), Slice(), outlet_layer, _stencil._output[1]}, 
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[1]}) + 1.0 / 6.0 * _lattice._outlet_density * 
                            ( -u_outlet + _lattice._uz.index({Slice(), Slice(), outlet_layer})) - Nx_z); // f[13]
      _lattice._f.index_put_({Slice(), Slice(), outlet_layer, _stencil._output[2]},
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[2]}) + 1.0 / 6.0 * _lattice._outlet_density * 
                            (-u_outlet - _lattice._uz.index({Slice(), Slice(), outlet_layer})) + Nx_z); // f[14]
      _lattice._f.index_put_({Slice(), Slice(), outlet_layer, _stencil._output[3]},
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[3]}) + 1.0 / 6.0 * _lattice._outlet_density * 
                            (-u_outlet + _lattice._uy.index({Slice(), Slice(), outlet_layer})) - Nx_y); // f[17]
      _lattice._f.index_put_({Slice(), Slice(), outlet_layer, _stencil._output[4]},
                            _lattice._f.index({Slice(), Slice(), outlet_layer, _stencil._input[4]}) + 1.0 / 6.0 * _lattice._outlet_density * 
                            (-u_outlet - _lattice._uy.index({Slice(), Slice(), outlet_layer})) + Nx_y); // f[18]
      break;                                                                                                                                                
  }
  setfSolidtoZero();
}

void
LatticeBoltzmannCore::setfSolidtoZero()
{
  /**
   * Set distribution functions to zero @ solid nodes
   */
  auto solid_mask = (_lattice._mesh == 0)
                        .unsqueeze(-1)
                        .expand({_lattice._nz, _lattice._ny, _lattice._nx, _lattice._q});
  ;

  _lattice._f.masked_fill_(solid_mask, 0);
}

void
LatticeBoltzmannCore::setfeqSolidtoZero()
{
  /**
   * Set equilibrium distribution functions to zero @ solid nodes
   */
  auto solid_mask = (_lattice._mesh == 0)
                        .unsqueeze(-1)
                        .expand({_lattice._nz, _lattice._ny, _lattice._nx, _lattice._q});
  ;

  _lattice._feq.masked_fill_(solid_mask, 0);
}

void
LatticeBoltzmannCore::setObservablestoZero()
{
  /**
   * Set observables to zero @ solid nodes
   */

  auto solid_mask = _lattice._mesh == 0;

  _lattice._rho.masked_fill_(solid_mask, 0);
  _lattice._ux.masked_fill_(solid_mask, 0);
  _lattice._uy.masked_fill_(solid_mask, 0);
  _lattice._uz.masked_fill_(solid_mask, 0);
  _lattice._u.masked_fill_(solid_mask, 0);
  _lattice._u_old.masked_fill_(solid_mask, 0);
}

void
LatticeBoltzmannCore::setStencil(const StencilBase & stencil_new)
{
  _stencil = stencil_new;
}

void
LatticeBoltzmannCore::setLatticeBase(const LatticeBase & lattice_new)
{
  _lattice = lattice_new;
}

void
LatticeBoltzmannCore::regularize()
{
  /**
   * Project non-equilibrium distribuin onto second order Hermite space
   */
  // To be implemented
}
