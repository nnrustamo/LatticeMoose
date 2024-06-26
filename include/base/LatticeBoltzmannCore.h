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
#include "LatticeStencil.h"

/**
 * Core functions of Lattice Boltzmann Method
 */

class LatticeBoltzmannCore
{
public:
    LatticeBoltzmannCore(const LatticeBase &lattice, const StencilBase &stencil);

    void generateMesh();

    void loadMeshFromFile(/*std::string &filename*/);

    void stream();
    
    void BGKCollision();

    void MRTCollision();

    void regularize();

    void equilibrium();

    void observables();

    void residual();

    void wallBoundary();

    void openBoundary();

    void setfSolidtoZero();

    void setfeqSolidtoZero();

    void setObservablestoZero();

    void setStencil(const StencilBase &stencil_new);

    void setLatticeBase(const LatticeBase &lattice_new);

    ~LatticeBoltzmannCore() = default;

public:
    LatticeBase _lattice;
    StencilBase _stencil;
    double _residual = 1.0;
    const double _fBody = 1.0e-4;
    bool _enableRegularization = false;
};
