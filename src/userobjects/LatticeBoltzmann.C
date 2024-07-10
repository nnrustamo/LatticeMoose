//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LatticeBoltzmann.h"
#include "Registry.h"
#include "LBMesh.h"

registerMooseObject("MooseApp", LatticeBoltzmann);

InputParameters
LatticeBoltzmann::validParams()
{
  InputParameters params = UserObject::validParams();
  params += MaterialPropertyInterface::validParams();
  params.addParam<double>("taus", 1.0, "Relaxation parameter");
  params.addParam<double>("initial_density", 1.0, "Initial lattice density");
  params.addParam<double>("inlet_density", 1.0, "Inlet lattice density");
  params.addParam<double>("outlet_density", 1.0, "Outlet lattice density");
  params.addParam<std::size_t>("n_subcycles", 1, "LBM iterations per timestep");
  params.addParam<double>("tolerance", 1.0e-6, "LBM convergence criteria");
  params.addParam<double>("fBody", 0.0, "LBM body force");

  return params;
}

LatticeBoltzmann::LatticeBoltzmann(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _simulation_object(LatticeBase(), StencilBase()),
    _n_subcycles(getParam<std::size_t>("n_subcycles")),
    _tolerance(getParam<double>("tolerance"))
{
  // get a reference to LBM mesh
  auto * mesh = dynamic_cast<LBMesh *>(&_subproblem.mesh());

  if (!mesh)
    mooseError("Must use LBMesh!");

  // Creating lattice
  long long nx = static_cast<long long>(mesh->getNx());
  long long ny = static_cast<long long>(mesh->getNy());
  long long nz = static_cast<long long>(mesh->getNz());
  MooseEnum dim = mesh->getDim();
  double taus = getParam<double>("taus");
  double initial_density = getParam<double>("initial_density");
  double inlet_density = getParam<double>("inlet_density");
  double outlet_density = getParam<double>("outlet_density");
  mapIndices();
  
  // Creating stencil
  StencilBase stencil;
  switch (dim)
  {
    case 2:
      stencil = D2Q9();
      _simulation_object._lattice.InitVars(nx, ny, nz, 9, taus, initial_density, inlet_density, outlet_density);
      break;
    case 3:
      stencil = D3Q19();
      _simulation_object._lattice.InitVars(nx, ny, nz, 19, taus, initial_density, inlet_density, outlet_density);
      break;
  }
  _simulation_object.setStencil(stencil);

  // intialize mesh
  if (mesh->load_mesh())
    _simulation_object._lattice._mesh = mesh->loadMeshFromFile();
  else
    _simulation_object.generateMesh();

  _simulation_object.processMesh();

  // printField(_simulation_object._lattice._mesh, 0);
  // initialize LBM simulation
  _simulation_object.computeEquilibrium();
  _simulation_object._lattice._f = _simulation_object._lattice._feq;
  _simulation_object.computeObservables();
}

void
LatticeBoltzmann::initialize()
{
  /**
   * Initialize timestep counter
   */
  _tsteps = 0;
}

void
LatticeBoltzmann::execute()
{
  /**
   * Main loop for the simulation
   * The order of steps are relatively flexible
   */ 

  while (/*_simulation_object._residual > _tolerance ||*/ _tsteps < _n_subcycles)
  { 
    if (_simulation_object._residual < _tolerance)
    {
        std::cout<< "Lattice Boltzmann simulation converged" << std::endl;
        break;
    }
    _simulation_object._lattice._u_old = _simulation_object._lattice._u; //.clone();
    _simulation_object.computeObservables();
    _simulation_object.computeEquilibrium();
    _simulation_object.MRTCollision();
    _simulation_object.stream();
    _simulation_object.wallBoundary();
    _simulation_object.openBoundary();
    _simulation_object.residual();
    logStep();
    _tsteps++;
  }
}

void
LatticeBoltzmann::finalize()
{
  /**
   * Finalize the simulation and write results
   */
}

void
LatticeBoltzmann::logStep()
{
  std::string msg = "Lattice Boltmann Timestep : " + std::to_string(_tsteps) +
                    ", Residual: " + std::to_string(_simulation_object._residual);
  std::cout << msg << std::endl;
}


void
LatticeBoltzmann::mapIndices()
{
  /**
   * Map the indices of the LBM mesh to the MOOSE mesh
   */
  auto * mesh = dynamic_cast<LBMesh *>(&_subproblem.mesh());
  if (!mesh)
    mooseError("Must use LBMesh!");

  // get the mesh bounds
  long long nx = static_cast<long long>(mesh->getNx());
  long long ny = static_cast<long long>(mesh->getNy());
  long long nz = static_cast<long long>(mesh->getNz());

  // map the indices
  for (unsigned long long id = 0; id < nx * ny * nz; id++)
  {
    std::vector<int> p(3);
    p[0] = id % nx;
    p[1] = (id / nx) % ny;
    p[2] = id / (nx * ny);
    _index_map[id] = p;
  }
}

Real
LatticeBoltzmann::getSpeed(unsigned long long p) const
{
  /**
   * Get the speed at a given point
   */
  return std::sqrt(
      Utility::pow<2>(_simulation_object._lattice._ux[_index_map.at(p)[2]][_index_map.at(p)[1]][_index_map.at(p)[0]].item<double>()) +
      Utility::pow<2>(_simulation_object._lattice._uy[_index_map.at(p)[2]][_index_map.at(p)[1]][_index_map.at(p)[0]].item<double>()) +
      Utility::pow<2>(_simulation_object._lattice._uz[_index_map.at(p)[2]][_index_map.at(p)[1]][_index_map.at(p)[0]].item<double>()));
}

Real
LatticeBoltzmann::getUx(unsigned long long p) const
{
  /**
   * Get the x-component of velocity at a given point
   */
  return _simulation_object._lattice._ux[_index_map.at(p)[2]][_index_map.at(p)[1]][_index_map.at(p)[0]].item<double>();
}

Real
LatticeBoltzmann::getUy(unsigned long long p) const
{
  /**
   * Get the y-component of velocity at a given point
   */
  return _simulation_object._lattice._uy[_index_map.at(p)[2]][_index_map.at(p)[1]][_index_map.at(p)[0]].item<double>();
}

Real
LatticeBoltzmann::getUz(unsigned long long p) const
{
  /**
   * Get the z-component of velocity at a given point
   */
  return _simulation_object._lattice._uz[_index_map.at(p)[2]][_index_map.at(p)[1]][_index_map.at(p)[0]].item<double>();
}

Real
LatticeBoltzmann::getDensity(unsigned long long p) const
{
  /**
   * Get pressure (density) at a given point
   */
  return _simulation_object._lattice._rho[_index_map.at(p)[2]][_index_map.at(p)[1]][_index_map.at(p)[0]].item<double>();
}


void LatticeBoltzmann::printField(const torch::Tensor & field, const unsigned int & precision)
{   
    /**
     * Printing the entire field for debugging
     */
    
    for(int64_t i = 0; i <_simulation_object._lattice._nz; i++)
    {
        for(int64_t j = 0; j < _simulation_object._lattice._ny; j++)
        {
            for(int64_t k = 0; k < _simulation_object._lattice._nx; k++)
            {
                std::cout<<std::setprecision(precision) << field[i][j][k].item<float>() << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

