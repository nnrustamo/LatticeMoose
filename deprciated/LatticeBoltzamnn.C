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

reisterMooseObject("MooseApp", LatticeBoltzmann);

InputParameters
LatticeBoltzmann::validParams()
{
    InputParameters params = UserObject::validParams();
    params += MaterialPropertyInterface::validParams();
    params.addParam<unsigned int>("n_subcycles", 1, "LBM iterations per timestep");
    params.addParam<double>("tolerance", 1.0e-6, "LBM convergence criteria");
    return params;
}

LatticeBoltzmann::LatticeBoltzmann(const InputParameters & parameters)
        : GeneralUserObject(parameters),
        _simulation_object(LatticeBase(parameters), StencilBase()),
        _n_subcycles(getParam<unsigned int>("n_subcycles"))
{

}

Real 
LatticeBoltzmann::getSpeed(Point p) const 
{
    return std::sqrt(Utility::pow<2>(_simulation_object._lattice._ux[p(0)][p(1)][p(2)].item<double>()) + 
    Utility::pow<2>(_simulation_object._lattice._uy[p(0)][p(1)][p(2)].item<double>()) + 
    Utility::pow<2>(_simulation_object._lattice._uz[p(0)][p(1)][p(2)].item<double>()) );
}
